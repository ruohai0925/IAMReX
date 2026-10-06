//
// Created by aoe on 25-7-7.
//


#include <AMReX_Math.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <AMReX_Utility.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_FillPatchUtil.H>
#include <iamr_constants.H>

#include "DiffusedIB_Parallel.H"
#include <AMReX_MPMD.H>

#include <algorithm>
#include <filesystem>
#include <sstream>
namespace fs = std::filesystem;

#define GHOST_CELLS 2

using namespace amrex;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                     global variable                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define LOCAL_LEVEL 0

const Vector<std::string> direction_str{"X","Y","Z"};

namespace ParticleProperties{
    Vector<Real> _x{}, _y{}, _z{}, _rho{};
    Vector<Real> Vx{}, Vy{}, Vz{};
    Vector<Real> Ox{}, Oy{}, Oz{};
    Vector<Real> _radius;
    Vector<Real> _radius2;
    Vector<Real> _radius3;
    Vector<int> _geometry_type;
    Real rd{0.0};
    Vector<int> TLX{}, TLY{},TLZ{},RLX{},RLY{},RLZ{};
    int euler_finest_level{0};
    int euler_velocity_index{0};
    int euler_force_index{0};
    Real euler_fluid_rho{0.0};
    int verbose{0};
    int loop_ns{2};
    int loop_solid{1};
    int Uhlmann{0};

    Vector<Real> GLO, GHI;
    int start_step{-1};
    int collision_model{0};
    int delta_type{1};
    int reduce_method{0};  // 0 = original ReduceSum loop, 1 = ParticleReduce 6-tuple + packed all-reduce
    int pvf_method{1};     // 0 = reference per-particle whole-domain phi/pvf path, 1 = bounding-box fused kernels + packed MPI

    int write_freq{1};
    bool init_particle_from_file{false};

    Real fin_length{0.0};
    Real fin_span{0.0};
    Real fin_amplitude_deg{0.0};
    Real fin_frequency{0.0};
    Real fin_wavelength{0.0};
    Real fin_phase{0.0};
    Real fin_wave_amplitude{0.0};
    Real fin_wave_number{1.0};
    Vector<Real> fin_wave_phase_deg;
    int fin_n_chord{0};
    int fin_n_span{0};

    GpuArray<Real, 3> plo{0.0,0.0,0.0}, phi{0.0,0.0,0.0}, dx{0.0, 0.0, 0.0};

    int RKPM{0};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_prescribed_fin(int geometry_type) noexcept
{
    return geometry_type == 3 || geometry_type == 4;
}

struct fin_marker_state {
    Real x;
    Real y;
    Real z;
    Real u;
    Real v;
    Real w;
    Real quadrature_scale;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
fin_marker_state evaluate_fin_marker(const kernel_gpu& body, Long marker_id)
{
    const int local_id = static_cast<int>(marker_id - body.start_id - 1);
    const int chord_id = local_id / body.fin_n_span;
    const int span_id = local_id - chord_id * body.fin_n_span;
    const Real s = body.fin_length * chord_id / Real(body.fin_n_chord - 1);

    fin_marker_state state{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    state.x = body.location[0] + s;
    if (chord_id == 0 || chord_id == body.fin_n_chord - 1) {
        state.quadrature_scale *= 0.5;
    }
    if (span_id == 0 || span_id == body.fin_n_span - 1) {
        state.quadrature_scale *= 0.5;
    }

    if (body.geometry_type == 3) {
        const Real r = body.fin_span * span_id / Real(body.fin_n_span - 1);
        const Real argument = 2.0 * Math::pi<Real>()
                            * body.fin_frequency * body.boundary_time
                            - 2.0 * Math::pi<Real>() * s / body.fin_wavelength
                            + body.fin_phase;
        const Real theta = body.fin_amplitude * std::sin(argument);
        const Real theta_dot = 2.0 * Math::pi<Real>() * body.fin_frequency
                             * body.fin_amplitude * std::cos(argument);
        state.y = body.location[1] + r * std::cos(theta);
        state.z = body.location[2] + r * std::sin(theta);
        state.v = -r * std::sin(theta) * theta_dot;
        state.w = r * std::cos(theta) * theta_dot;
    } else {
        const Real r = -0.5 * body.fin_span
                     + body.fin_span * span_id / Real(body.fin_n_span - 1);
        const Real argument = 2.0 * Math::pi<Real>() * body.fin_wave_number
                            * s / body.fin_length
                            - 2.0 * Math::pi<Real>() * body.fin_frequency
                            * body.boundary_time + body.fin_phase;
        const Real slope = body.fin_wave_amplitude * 2.0 * Math::pi<Real>()
                         * body.fin_wave_number / body.fin_length
                         * std::cos(argument);
        state.y = body.location[1]
                + body.fin_wave_amplitude * std::sin(argument);
        state.z = body.location[2] + r;
        state.v = -2.0 * Math::pi<Real>() * body.fin_frequency
                * body.fin_wave_amplitude * std::cos(argument);
        state.quadrature_scale *= std::sqrt(1.0 + slope * slope);
    }
    return state;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                     other function                            */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
AMREX_INLINE AMREX_GPU_DEVICE
Real nodal_phi_to_heavi(Real phi)
{
    if (phi <= 0.0) {
        return 0.0;
    }
    return 1.0;
}

// Nodal level-set value of one particle at node (Xn,Yn,Zn). Same formulas as
// calculate_phi_nodal(): sphere -> signed distance normalised by the radius,
// ellipsoid -> first-order approximation left unnormalised, anything else -> 0.
AMREX_INLINE AMREX_GPU_DEVICE
Real nodal_phi_at (Real Xn, Real Yn, Real Zn,
                   Real Xp, Real Yp, Real Zp,
                   Real a, Real b, Real c, int geometry_type) noexcept
{
    if (geometry_type == 1) {
        Real phi = std::sqrt( (Xn - Xp)*(Xn - Xp)
                            + (Yn - Yp)*(Yn - Yp) + (Zn - Zp)*(Zn - Zp)) - a;
        return phi / a;
    } else if (geometry_type == 2) {
        Real xp = Xn - Xp;
        Real yp = Yn - Yp;
        Real zp = Zn - Zp;
        Real xp2 = xp * xp;
        Real yp2 = yp * yp;
        Real zp2 = zp * zp;
        Real a2 = a * a;
        Real b2 = b * b;
        Real c2 = c * c;
        Real numerator = (xp2 / a2 + yp2 / b2 + zp2 / c2) - 1.0;
        Real xp4 = xp2 * xp2;
        Real yp4 = yp2 * yp2;
        Real zp4 = zp2 * zp2;
        Real a4 = a2 * a2;
        Real b4 = b2 * b2;
        Real c4 = c2 * c2;
        Real denominator = 2.0 * std::sqrt(xp4 / a4 + yp4 / b4 + zp4 / c4);
        return numerator / (denominator + 1.e-12);
    }
    return 0.0;
}

// Volume fraction of cell (i,j,k) from the 8 analytic nodal values of one
// particle. Same arithmetic (and summation order) as calculate_phi_nodal()
// followed by nodal_phi_to_pvf(), without materialising phi_nodal.
AMREX_INLINE AMREX_GPU_DEVICE
Real cell_pvf_at (int i, int j, int k,
                  GpuArray<Real,3> const& dx, GpuArray<Real,3> const& plo,
                  Real Xp, Real Yp, Real Zp,
                  Real a, Real b, Real c, int geometry_type) noexcept
{
    Real num = 0.0;
    Real deo = 0.0;
    for (int kk = k; kk <= k+1; ++kk) {
        for (int jj = j; jj <= j+1; ++jj) {
            for (int ii = i; ii <= i+1; ++ii) {
                const Real ph = nodal_phi_at(ii * dx[0] + plo[0], jj * dx[1] + plo[1], kk * dx[2] + plo[2],
                                             Xp, Yp, Zp, a, b, c, geometry_type);
                num += (-ph) * nodal_phi_to_heavi(-ph);
                deo += std::abs(ph);
            }
        }
    }
    return num / (deo + 1.e-12);
}

// Index-space box containing every cell whose pvf can be non-zero for a body
// of extent R centred at loc: such a cell has a node strictly inside the body,
// i.e. within R of the centre. One extra cell of margin on each side.
Box particle_index_box (const RealVect& loc, Real R)
{
    const auto& dx  = ParticleProperties::dx;
    const auto& plo = ParticleProperties::plo;
    IntVect lo, hi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        lo[d] = static_cast<int>(std::floor((loc[d] - R - plo[d]) / dx[d])) - 1;
        hi[d] = static_cast<int>(std::ceil ((loc[d] + R - plo[d]) / dx[d])) + 1;
    }
    return Box(lo, hi);
}

void nodal_phi_to_pvf(MultiFab& pvf, const MultiFab& phi_nodal)
{
    BL_PROFILE("DiffusedIB::nodal_phi_to_pvf");
    // Print() << "In the nodal_phi_to_pvf\n";

    pvf.setVal(0.0);

    // Only set the valid cells of pvf
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(pvf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& pvffab   = pvf.array(mfi);
        auto const& pnfab = phi_nodal.array(mfi);
        ParallelFor(bx, [pvffab, pnfab]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real num = 0.0;
            for(int kk=k; kk<=k+1; kk++) {
                for(int jj=j; jj<=j+1; jj++) {
                    for(int ii=i; ii<=i+1; ii++) {
                        num += (-pnfab(ii,jj,kk)) * nodal_phi_to_heavi(-pnfab(ii,jj,kk));
                    }
                }
            }
            Real deo = 0.0;
            for(int kk=k; kk<=k+1; kk++) {
                for(int jj=j; jj<=j+1; jj++) {
                    for(int ii=i; ii<=i+1; ii++) {
                        deo += std::abs(pnfab(ii,jj,kk));
                    }
                }
            }
            pvffab(i,j,k) = num / (deo + 1.e-12);
        });
    }

}

// Fill phi_nodal with the nodal level set of a single particle: negative inside
// the body, positive outside. Handles geometry_type = 1 (sphere, signed distance
// normalised by the radius) and geometry_type = 2 (ellipsoid, a first-order
// approximation that is deliberately left unnormalised). Any other value aborts.
void calculate_phi_nodal(MultiFab& phi_nodal, kernel& current_kernel)
{
    BL_PROFILE("DiffusedIB::calculate_phi_nodal");
    phi_nodal.setVal(0.0);

    Real Xp = current_kernel.location[0];
    Real Yp = current_kernel.location[1];
    Real Zp = current_kernel.location[2];
    Real a = current_kernel.radius;
    int geometry_type = current_kernel.geometry_type;

    // Only set the valid cells of phi_nodal
    for (MFIter mfi(phi_nodal,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& pnfab = phi_nodal.array(mfi);
        auto dx = ParticleProperties::dx;
        auto plo = ParticleProperties::plo;

        if (geometry_type == 1) {
            // Sphere geometry
            ParallelFor(bx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    Real Xn = i * dx[0] + plo[0];
                    Real Yn = j * dx[1] + plo[1];
                    Real Zn = k * dx[2] + plo[2];

                    pnfab(i,j,k) = std::sqrt( (Xn - Xp)*(Xn - Xp)
                            + (Yn - Yp)*(Yn - Yp)  + (Zn - Zp)*(Zn - Zp)) - a;
                    pnfab(i,j,k) = pnfab(i,j,k) / a;

                }
            );
        } else if (geometry_type == 2) {
            // Ellipsoid geometry
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(current_kernel.radius2 > 0.0 && current_kernel.radius3 > 0.0,
                "geometry_type = 2 (ellipsoid) requires radius2 and radius3 to be set to positive values");
            Real b = current_kernel.radius2;  // semi-axis b
            Real c = current_kernel.radius3; // semi-axis c
            // Floor the denominator relative to the smallest semi-axis so that the
            // centre of the particle stays representable, including in FP32 builds.
            const Real denom_floor = Real(1.e-6) * amrex::min(a, amrex::min(b, c));
            ParallelFor(bx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    Real Xn = i * dx[0] + plo[0];
                    Real Yn = j * dx[1] + plo[1];
                    Real Zn = k * dx[2] + plo[2];

                    // Relative coordinates to ellipsoid center
                    Real xp = Xn - Xp;
                    Real yp = Yn - Yp;
                    Real zp = Zn - Zp;

                    // Compute ellipsoid level set function using the formula:
                    // d ≈ (x'^2/a^2 + y'^2/b^2 + z'^2/c^2 - 1) / (2 * sqrt(x'^4/a^4 + y'^4/b^4 + z'^4/c^4))
                    Real xp2 = xp * xp;
                    Real yp2 = yp * yp;
                    Real zp2 = zp * zp;

                    Real a2 = a * a;
                    Real b2 = b * b;
                    Real c2 = c * c;

                    Real numerator = (xp2 / a2 + yp2 / b2 + zp2 / c2) - 1.0;

                    Real xp4 = xp2 * xp2;
                    Real yp4 = yp2 * yp2;
                    Real zp4 = zp2 * zp2;

                    Real a4 = a2 * a2;
                    Real b4 = b2 * b2;
                    Real c4 = c2 * c2;

                    Real denominator = 2.0 * std::sqrt(xp4 / a4 + yp4 / b4 + zp4 / c4);

                    // Do not normalize here!
                    pnfab(i,j,k) = numerator / amrex::max(denominator, denom_floor);

                }
            );
        } else {
            Print() << "Particle (" << current_kernel.id << ") has unsupported geometry_type: " << geometry_type << "\n";
            Abort("Unsupported geometry type. Only geometry_type = 1 (sphere) and 2 (ellipsoid) are supported.");
        }
    }
}

// May use ParReduce later, https://amrex-codes.github.io/amrex/docs_html/GPU.html#multifab-reductions
void CalculateSumU_cir (RealVect& sum,
                        const MultiFab& E,
                        const MultiFab& pvf,
                        int EulerVelIndex)
{
    auto const& E_data = E.const_arrays();
    auto const& pvf_data = pvf.const_arrays();
    const Real d = Math::powi<3>(ParticleProperties::dx[0]);
    GpuTuple<Real, Real, Real> tmpSum = ParReduce(TypeList<ReduceOpSum,ReduceOpSum,ReduceOpSum>{}, TypeList<Real, Real, Real>{},E, IntVect{0},
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real, Real, Real>{
        auto E_ = E_data[box_no];
        auto pvf_ = pvf_data[box_no];
        return {
            E_(i, j, k, EulerVelIndex    ) * d * pvf_(i,j,k),
            E_(i, j, k, EulerVelIndex + 1) * d * pvf_(i,j,k),
            E_(i, j, k, EulerVelIndex + 2) * d * pvf_(i,j,k)
        };
    });
    sum[0] = get<0>(tmpSum);
    sum[1] = get<1>(tmpSum);
    sum[2] = get<2>(tmpSum);
}

void CalculateSumT_cir (RealVect& sum,
                        const MultiFab& E,
                        const MultiFab& pvf,
                        const RealVect pLoc,
                        int EulerVelIndex)
{
    auto plo = ParticleProperties::plo;
    auto dx = ParticleProperties::dx;

    auto const& E_data = E.const_arrays();
    auto const& pvf_data = pvf.const_arrays();
    const Real d = Math::powi<3>(ParticleProperties::dx[0]);
    GpuTuple<Real, Real, Real> tmpSum = ParReduce(TypeList<ReduceOpSum,ReduceOpSum,ReduceOpSum>{}, TypeList<Real, Real, Real>{},E, IntVect{0},
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real, Real, Real>{
        auto E_ = E_data[box_no];
        auto pvf_ = pvf_data[box_no];

        Real x = plo[0] + i*dx[0] + 0.5*dx[0];
        Real y = plo[1] + j*dx[1] + 0.5*dx[1];
        Real z = plo[2] + k*dx[2] + 0.5*dx[2];

        Real vx = E_(i, j, k, EulerVelIndex    );
        Real vy = E_(i, j, k, EulerVelIndex + 1);
        Real vz = E_(i, j, k, EulerVelIndex + 2);

        RealVect tmp = RealVect(x - pLoc[0], y - pLoc[1], z - pLoc[2]).crossProduct(RealVect(vx, vy, vz));

        return {
            tmp[0] * d * pvf_(i, j, k),
            tmp[1] * d * pvf_(i, j, k),
            tmp[2] * d * pvf_(i, j, k)
        };
    });
    sum[0] = get<0>(tmpSum);
    sum[1] = get<1>(tmpSum);
    sum[2] = get<2>(tmpSum);
}

[[nodiscard]] AMREX_FORCE_INLINE
Real cal_momentum(Real rho, Real radius, int geometry_type = 1, int idir = 0, Real radius2 = 0.0, Real radius3 = 0.0)
{
    if (geometry_type == 1) {
        // Sphere: I = (2/5) * m * r² = (8/15) * π * ρ * r⁵
        return 8.0 * Math::pi<Real>() * rho * Math::powi<5>(radius) / 15.0;
    } else if (geometry_type == 2) {
        // Ellipsoid: I = (1/5) * m * (sum of squares of perpendicular semi-axes)
        // m = (4/3) * π * a * b * c * ρ
        Real a = radius;
        Real b = (radius2 > 0.0) ? radius2 : radius;
        Real c = (radius3 > 0.0) ? radius3 : radius;
        Real m = 4.0 * Math::pi<Real>() * rho * a * b * c / 3.0;

        // Moment of inertia depends on rotation axis
        Real I;
        if (idir == 0) {
            // Rotation around x-axis (a-axis): I_x = (1/5) * m * (b² + c²)
            I = m * (b * b + c * c) / 5.0;
        } else if (idir == 1) {
            // Rotation around y-axis (b-axis): I_y = (1/5) * m * (a² + c²)
            I = m * (a * a + c * c) / 5.0;
        } else {
            // Rotation around z-axis (c-axis): I_z = (1/5) * m * (a² + b²)
            I = m * (a * a + b * b) / 5.0;
        }
        return I;
    } else {
        // Default to sphere for unsupported geometry types
        return 8.0 * Math::pi<Real>() * rho * Math::powi<5>(radius) / 15.0;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void deltaFunction(Real xf, Real xp, Real h, Real& value, int type)
{
    Real rr = Math::abs(( xf - xp ) / h);

    switch (type) {
    case DELTA_FUNCTION_TYPE::FOUR_POINT_IB:
        if(rr >= 0 && rr < 1 ){
            value = 1.0 / 8.0 * ( 3.0 - 2.0 * rr + std::sqrt( 1.0 + 4.0 * rr - 4 * Math::powi<2>(rr))) / h;
        }else if (rr >= 1 && rr < 2) {
            value = 1.0 / 8.0 * ( 5.0 - 2.0 * rr - std::sqrt( -7.0 + 12.0 * rr - 4 * Math::powi<2>(rr))) / h;
        }else {
            value = 0;
        }
        break;
    case DELTA_FUNCTION_TYPE::THREE_POINT_IB:
        if(rr >= 0.5 && rr < 1.5){
            value = 1.0 / 6.0 * ( 5.0 - 3.0 * rr - std::sqrt( - 3.0 * Math::powi<2>( 1 - rr) + 1.0 )) / h;
        }else if (rr >= 0 && rr < 0.5) {
            value = 1.0 / 3.0 * ( 1.0 + std::sqrt( 1.0 - 3 * Math::powi<2>(rr))) / h;
        }else {
            value = 0;
        }
        break;
    default:
        break;
    }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                    mParticle member function                  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//loop all particels
void mParticle::InteractWithEuler(MultiFab &EulerVel,
                                  MultiFab &EulerForce,
                                  Real dt)
{
    BL_PROFILE("mParticle::InteractWithEuler");
    BL_PROFILE_VAR("mParticle::InteractWithEuler()", mParticle_InteractWithEuler);
    if (verbose) Print() << "[Particle] mParticle::InteractWithEuler\n";
    // clear time , start record
    spend_time = 0;
    const auto InteractWithEulerStart = ParallelDescriptor::second();
    // The parallel IBM interface does not receive the simulation time. Advance
    // the prescribed boundary to the same new-time level as EulerVel. The
    // value is re-synchronised from UpdateParticles after every fluid step.
    ib_time += dt;

    //clean preStep's IB_properties
    for(auto& kernel : particle_kernels) {
        kernel.ib_force.scale(0.0);
        kernel.ib_moment.scale(0.0);
    }
    int type = ParticleProperties::delta_type;
    //for 1 -> Ns是
    int loop = ParticleProperties::loop_ns;

    // Sync kernel data to device memory for GPU ParallelFor
    syncKernelsToDevice();
    RefreshRKPMWeights();

    BL_ASSERT(loop > 0);
    while(loop > 0){
        if(verbose) Print() << "[Particle] Ns loop index : " << loop << "\n";

        EulerForce.setVal(0.0);
        // 清空拉格朗日点的所有的信息，并更新位置
        UpdateLagrangianMarker();
        // 速度插值
        VelocityInterpolation(EulerVel, type);
        // IB力
        ComputeLagrangianForce(dt);
        // 弥散
        ForceSpreading(EulerForce, type);
        // 速度修正
        VelocityCorrection(EulerVel, EulerForce, dt);
        loop--;
    }
    BL_PROFILE_VAR_STOP(mParticle_InteractWithEuler);
    spend_time += ParallelDescriptor::second() - InteractWithEulerStart;
}

void mParticle::RefreshRKPMWeights()
{
    if (!do_RKPM || !amrex::MPMD::Initialized() || m_mpmd_exchange_done) return;
    // m_mpmd_exchange_done = true;

    const int Nm = static_cast<int>(LargrangianMarker.size());
    const int server_root = amrex::MPMD::NProcs() - 1;
    const int nw = Nm * RKPM_STENCIL_SIZE;

    Vector<float> h_w(nw);
    if (ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()) {
        const Real t0 = ParallelDescriptor::second();
        // send number of markers
        MPI_Send(&Nm, 1, MPI_INT, server_root, 200, MPI_COMM_WORLD);
        // send marker positions
        Vector<double> h_pos(Nm * 3);
        for (int i = 0; i < Nm; ++i) {
            h_pos[3*i+0] = static_cast<double>(LargrangianMarker[i][0]);
            h_pos[3*i+1] = static_cast<double>(LargrangianMarker[i][1]);
            h_pos[3*i+2] = static_cast<double>(LargrangianMarker[i][2]);
        }
        MPI_Send(h_pos.data(), Nm*3, MPI_DOUBLE, server_root, 201, MPI_COMM_WORLD);
        // receive weights
        MPI_Recv(h_w.data(), nw, MPI_FLOAT, server_root, 202, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        const Real t1 = ParallelDescriptor::second();
        double s = 0.0;
        for (int i = 0; i < nw; ++i) s += h_w[i];
        Print() << "[RKPM-MPMD] received " << Nm << " markers, sum=" << s
                << ", recv-send time (ML inference included)=" << (t1 - t0) << " s\n";
    }

    // broadcast weights to all procs
    ParallelDescriptor::Bcast(h_w.data(), nw, ParallelDescriptor::IOProcessorNumber());

    // copy to device weight
    AMREX_ALWAYS_ASSERT(h_rkpm_flat.size() == static_cast<size_t>(nw));
    for (int i = 0; i < nw; ++i) {
        h_rkpm_flat[i].weight = static_cast<Real>(h_w[i]);
    }
    Gpu::copyAsync(Gpu::hostToDevice, h_rkpm_flat.begin(), h_rkpm_flat.end(), d_rkpm_flat.begin());
    Gpu::streamSynchronize();

    if (ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()) {
        Print() << "[RKPM-MPMD] wrote " << nw << " weights to d_rkpm_flat\n";
    }
}

void mParticle::InitParticles(const Vector<Real>& x,
                              const Vector<Real>& y,
                              const Vector<Real>& z,
                              const Vector<Real>& rho_s,
                              const Vector<Real>& Vx,
                              const Vector<Real>& Vy,
                              const Vector<Real>& Vz,
                              const Vector<Real>& Ox,
                              const Vector<Real>& Oy,
                              const Vector<Real>& Oz,
                              const Vector<int>& TLXt,
                              const Vector<int>& TLYt,
                              const Vector<int>& TLZt,
                              const Vector<int>& RLXt,
                              const Vector<int>& RLYt,
                              const Vector<int>& RLZt,
                              const Vector<Real>& radius,
                              const Vector<Real>& radius2,
                              const Vector<Real>& radius3,
                              const Vector<int>& geometry_type,
                              Real h,
                              Real gravity,
                              int _verbose)
{
    verbose = _verbose;
    if (verbose) Print() << "[Particle] mParticle::InitParticles\n";

    m_gravity[2] = gravity;

    //pre judge
    if(!((x.size() == y.size()) && (x.size() == z.size()))){
        Print() << "particle's position container are all different size";
        return;
    }

    const bool has_prescribed_fin = std::any_of(
        geometry_type.begin(), geometry_type.end(),
        [](int type) { return is_prescribed_fin(type); });
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !(has_prescribed_fin && ParticleProperties::RKPM != 0),
        "Prescribed fins do not support RKPM markers");

    if (ParticleProperties::RKPM != 0) {
        // RKPM only one particle
        do_RKPM = true;
        ResolveLagrangianMarker("rkpm_mappings.id");
        ResolveWithRPKM("rkpm_mappings.lag");
    }

    //all the particles have different radius
    for(int index = 0; index < x.size(); index++){
        int real_index;
        // if initial with input file, initialized by [0] data
        if(ParticleProperties::init_particle_from_file){
            real_index = 0;
        }else{
            real_index = index;
        }

        kernel mKernel;
        mKernel.id = index;
        mKernel.location[0] = x[index];
        mKernel.location[1] = y[index];
        mKernel.location[2] = z[index];
        mKernel.velocity[0] = Vx[real_index];
        mKernel.velocity[1] = Vy[real_index];
        mKernel.velocity[2] = Vz[real_index];
        mKernel.omega[0] = Ox[real_index];
        mKernel.omega[1] = Oy[real_index];
        mKernel.omega[2] = Oz[real_index];

        // use current state to initialize old state
        mKernel.location_old = mKernel.location;
        mKernel.velocity_old = mKernel.velocity;
        mKernel.omega_old = mKernel.omega;

        mKernel.TL[0] = TLXt[real_index];
        mKernel.TL[1] = TLYt[real_index];
        mKernel.TL[2] = TLZt[real_index];
        mKernel.RL[0] = RLXt[real_index];
        mKernel.RL[1] = RLYt[real_index];
        mKernel.RL[2] = RLZt[real_index];
        mKernel.rho = rho_s[real_index];
        // geometry type
        if (geometry_type.size() > 0 && real_index < geometry_type.size()) {
            mKernel.geometry_type = geometry_type[real_index];
        } else {
            mKernel.geometry_type = 1;  // default to sphere if not provided
        }
        // sphere particle need
        mKernel.radius = radius[real_index];
        // ellipsoid particle need
        // Check if radius2 is provided and has enough elements
        if (radius2.size() > 0 && real_index < radius2.size()) {
            mKernel.radius2 = radius2[real_index];
        } else {
            mKernel.radius2 = radius[real_index];  // default to radius if not provided
        }
        // Check if radius3 is provided and has enough elements
        if (radius3.size() > 0 && real_index < radius3.size()) {
            mKernel.radius3 = radius3[real_index];
        } else {
            mKernel.radius3 = radius[real_index];  // default to radius if not provided
        }
        if (is_prescribed_fin(mKernel.geometry_type)) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ParticleProperties::RKPM == 0,
                "Prescribed fins do not support RKPM markers");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ParticleProperties::fin_length > 0.0,
                "A prescribed fin requires fin_length > 0");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ParticleProperties::fin_span > 0.0,
                "A prescribed fin requires fin_span > 0");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ParticleProperties::fin_n_chord >= 2
                && ParticleProperties::fin_n_span >= 2,
                "A prescribed fin requires fin_n_chord and fin_n_span >= 2");

            mKernel.fin_length = ParticleProperties::fin_length;
            mKernel.fin_span = ParticleProperties::fin_span;
            mKernel.fin_frequency = ParticleProperties::fin_frequency;
            mKernel.fin_n_chord = ParticleProperties::fin_n_chord;
            mKernel.fin_n_span = ParticleProperties::fin_n_span;

            if (mKernel.geometry_type == 3) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ParticleProperties::fin_wavelength > 0.0,
                    "A twisting fin requires fin_wavelength > 0");
                mKernel.fin_amplitude = ParticleProperties::fin_amplitude_deg
                                      * Math::pi<Real>() / 180.0;
                mKernel.fin_wavelength = ParticleProperties::fin_wavelength;
                mKernel.fin_phase = ParticleProperties::fin_phase;
            } else {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    ParticleProperties::fin_wave_amplitude >= 0.0,
                    "A transverse-wave fin requires fin_wave_amplitude >= 0");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ParticleProperties::fin_wave_number > 0.0,
                    "A transverse-wave fin requires fin_wave_number > 0");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    ParticleProperties::fin_wave_phase_deg.empty()
                    || real_index < static_cast<int>(ParticleProperties::fin_wave_phase_deg.size()),
                    "fin_wave_phase_deg must provide one value per fin");
                mKernel.fin_wave_amplitude = ParticleProperties::fin_wave_amplitude;
                mKernel.fin_wave_number = ParticleProperties::fin_wave_number;
                mKernel.fin_phase = ParticleProperties::fin_wave_phase_deg.empty()
                                  ? ParticleProperties::fin_phase
                                  : ParticleProperties::fin_wave_phase_deg[real_index]
                                    * Math::pi<Real>() / 180.0;
            }

            mKernel.ml = mKernel.fin_n_chord * mKernel.fin_n_span;
            mKernel.dv = h
                       * mKernel.fin_length / Real(mKernel.fin_n_chord - 1)
                       * mKernel.fin_span / Real(mKernel.fin_n_span - 1);
            mKernel.Vp = 0.0;
            max_largrangian_num = std::max(max_largrangian_num, mKernel.ml);
            particle_kernels.emplace_back(mKernel);
            continue;
        }

        mKernel.Vp = Math::pi<Real>() * 4 / 3 * Math::powi<3>(radius[real_index]);

        if (ParticleProperties::RKPM == 0) {

            //int Ml = static_cast<int>( Math::pi<Real>() / 3 * (12 * Math::powi<2>(mKernel.radius / h)));
            //Real dv = Math::pi<Real>() * h / 3 / Ml * (12 * mKernel.radius * mKernel.radius + h * h);
            const int Ml = static_cast<int>((Math::powi<3>(mKernel.radius - (ParticleProperties::rd - 0.5) * h)
                   - Math::powi<3>(mKernel.radius - (ParticleProperties::rd + 0.5) * h))/(3.*h*h*h/4./Math::pi<Real>()));
            Real dv = (Math::powi<3>(mKernel.radius - (ParticleProperties::rd - 0.5) * h)
                   - Math::powi<3>(mKernel.radius - (ParticleProperties::rd + 0.5) * h))/(3.*Ml/4./Math::pi<Real>());
            mKernel.ml = Ml;
            mKernel.dv = dv;
            if( Ml > max_largrangian_num ) max_largrangian_num = Ml;

            Real phiK = 0;
            Gpu::HostVector<Real> h_phiK, h_thetaK;
            for(int marker_index = 0; marker_index < Ml; marker_index++){
                const Real Hk = -1.0 + 2.0 * (marker_index) / ( Ml - 1.0);
                Real thetaK = std::acos(Hk);
                if(marker_index == 0 || marker_index == Ml - 1){
                    phiK = 0;
                }else {
                    phiK = std::fmod( phiK + 3.809 / std::sqrt(Ml) / std::sqrt( 1 - Math::powi<2>(Hk)) , 2 * Math::pi<Real>());
                }
                h_phiK.push_back(phiK);
                h_thetaK.push_back(thetaK);
            }
            mKernel.phiK.resize(Ml);
            mKernel.thetaK.resize(Ml);
            Gpu::copyAsync(Gpu::hostToDevice, h_phiK.begin(), h_phiK.end(), mKernel.phiK.begin());
            Gpu::copyAsync(Gpu::hostToDevice, h_thetaK.begin(), h_thetaK.end(), mKernel.thetaK.begin());
            Gpu::streamSynchronize();

            if (verbose) Print() << "h: " << h << ", Ml: " << Ml << ", D: " << Math::powi<3>(h) << " gravity : " << gravity << "\n"
                                        << "Kernel : " << index << ": Location (" << x[index] << ", " << y[index] << ", " << z[index]
                                        << "), Velocity : (" << mKernel.velocity[0] << ", " << mKernel.velocity[1] << ", "<< mKernel.velocity[2]
                                        << "), Radius: " << mKernel.radius << ", Radius2: " << mKernel.radius2 << ", Radius3: " << mKernel.radius3 << ", Ml: " << Ml << ", dv: " << dv << ", Rho: " << mKernel.rho << "\n";
        }else {
            mKernel.ml = NumOfLagrangianMarker(index);
            mKernel.dv = RKPM_MAP.at(index)[0].eps;
        }

        particle_kernels.emplace_back(mKernel);
    }
    // Keep the original collision setup for rigid bodies. Prescribed fins are
    // open surfaces and do not participate in the rigid-body collision model.
    const auto rigid = std::find_if(particle_kernels.begin(), particle_kernels.end(),
        [](const kernel& body) { return !is_prescribed_fin(body.geometry_type); });
    if (rigid != particle_kernels.end()) {
        m_Collision.SetGeometry(RealVect(ParticleProperties::GLO),
                                RealVect(ParticleProperties::GHI),
                                rigid->radius, h);
    }
}

void mParticle::syncKernelsToDevice ()
{
    const int nk = static_cast<int>(particle_kernels.size());
    Gpu::HostVector<kernel_gpu> h_kg(nk);
    for (int i = 0; i < nk; ++i) {
        auto const& pk = particle_kernels[i];
        h_kg[i].location  = pk.location;
        h_kg[i].velocity  = pk.velocity;
        h_kg[i].omega     = pk.omega;
        h_kg[i].radius    = pk.radius;
        h_kg[i].dv        = pk.dv;
        h_kg[i].phiK      = pk.phiK.data();
        h_kg[i].thetaK    = pk.thetaK.data();
        h_kg[i].start_id  = pk.start_id;
        h_kg[i].geometry_type = pk.geometry_type;
        h_kg[i].fin_length = pk.fin_length;
        h_kg[i].fin_span = pk.fin_span;
        h_kg[i].fin_amplitude = pk.fin_amplitude;
        h_kg[i].fin_frequency = pk.fin_frequency;
        h_kg[i].fin_wavelength = pk.fin_wavelength;
        h_kg[i].fin_phase = pk.fin_phase;
        h_kg[i].fin_wave_amplitude = pk.fin_wave_amplitude;
        h_kg[i].fin_wave_number = pk.fin_wave_number;
        h_kg[i].fin_n_chord = pk.fin_n_chord;
        h_kg[i].fin_n_span = pk.fin_n_span;
        h_kg[i].boundary_time = ib_time;
    }
    d_kernels.resize(nk);
    Gpu::copyAsync(Gpu::hostToDevice, h_kg.begin(), h_kg.end(), d_kernels.begin());
    Gpu::streamSynchronize();
}

void mParticle::UpdateLagrangianMarker() {
    BL_PROFILE("mParticle::UpdateLagrangianMarker");
    if (verbose) Print() << "\tmParticle::UpdateLagrangianMarker\n";
    // update lagrangian marker attributions
    for (mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti) {
        const Box& box = pti.validbox();

        auto* particles = pti.GetArrayOfStructs().data();
        const Long np = pti.numParticles();
        auto* attri = pti.GetAttribs().data();
        const auto *const ids = pti.GetIDs().data();
        auto *const vUP_ptr = attri[P_ATTR_REAL::U_Marker].data();
        auto *const vVP_ptr = attri[P_ATTR_REAL::V_Marker].data();
        auto *const vWP_ptr = attri[P_ATTR_REAL::W_Marker].data();
        auto *const fxP_ptr = attri[P_ATTR_REAL::Fx_Marker].data();
        auto *const fyP_ptr = attri[P_ATTR_REAL::Fy_Marker].data();
        auto *const fzP_ptr = attri[P_ATTR_REAL::Fz_Marker].data();
        auto *const mxP_ptr = attri[P_ATTR_REAL::Mx_Marker].data();
        auto *const myP_ptr = attri[P_ATTR_REAL::My_Marker].data();
        auto *const mzP_ptr = attri[P_ATTR_REAL::Mz_Marker].data();

        const auto *const ps = d_kernels.data();
        const bool do_RKPM_l = do_RKPM;

        ParallelFor(np,
            [=] AMREX_GPU_DEVICE (const int i) noexcept {
                if (!do_RKPM_l) {
                    const auto id = ids[i];
                    const auto m_id = particles[i].id();
                    const auto location = ps[id].location;
                    const auto radius = ps[id].radius;
                    const auto *const phiK = ps[id].phiK;
                    const auto *const thetaK = ps[id].thetaK;
                    const auto start_id = ps[id].start_id;
                    const auto geometry_type = ps[id].geometry_type;

                    if (is_prescribed_fin(geometry_type)) {
                        const auto state = evaluate_fin_marker(ps[id], m_id);
                        particles[i].pos(0) = state.x;
                        particles[i].pos(1) = state.y;
                        particles[i].pos(2) = state.z;
                    } else {
                        const auto ia = m_id - start_id - 1;
                        particles[i].pos(0) = location[0] + radius * std::sin(thetaK[ia]) * std::cos(phiK[ia]);
                        particles[i].pos(1) = location[1] + radius * std::sin(thetaK[ia]) * std::sin(phiK[ia]);
                        particles[i].pos(2) = location[2] + radius * std::cos(thetaK[ia]);
                    }
                }
                // RKPM forbi
                vUP_ptr[i] = 0.0;
                vVP_ptr[i] = 0.0;
                vWP_ptr[i] = 0.0;
                fxP_ptr[i] = 0.0;
                fyP_ptr[i] = 0.0;
                fzP_ptr[i] = 0.0;
                mxP_ptr[i] = 0.0;
                myP_ptr[i] = 0.0;
                mzP_ptr[i] = 0.0;
            }
        );
    }
    // https://amrex-codes.github.io/amrex/docs_html/Particle.html#redistribute
    // Rigid markers need redistribution when a centre moves. Prescribed fins
    // also need it as their surfaces deform around a fixed origin.
    if (MarkersNeedRedistribute()) {
        mContainer->Redistribute();
    }

    if (verbose) {
        Print() << "[particle] : particle num :" << mContainer->TotalNumberOfParticles() << "\n";
        // mContainer->WriteAsciiFile(Concatenate("particle", 1));
    }
}

// Redistribute deforming fins, moved rigid bodies, and newly created markers.
// Remember the current centres for the rigid-body fast path.
bool mParticle::MarkersNeedRedistribute()
{
    const int nk = static_cast<int>(particle_kernels.size());
    // A prescribed fin deforms even when its origin stays fixed.
    bool moved = (static_cast<int>(m_redistributed_loc.size()) != nk)
              || std::any_of(particle_kernels.begin(), particle_kernels.end(),
                  [](const kernel& body) { return is_prescribed_fin(body.geometry_type); });
    for (int k = 0; !moved && k < nk; ++k) {
        moved = (particle_kernels[k].location != m_redistributed_loc[k]);
    }
    if (moved) {
        m_redistributed_loc.resize(nk);
        for (int k = 0; k < nk; ++k) m_redistributed_loc[k] = particle_kernels[k].location;
    }
    return moved;
}

template <typename P = Particle<num_Real, num_Int>>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void VelocityInterpolation_cir(int p_iter, P const& p, ParticleReal& Up, ParticleReal& Vp, ParticleReal& Wp,
                               Array4<Real const> const& E, int EulerVIndex,
                               const int *lo, const int *hi,
                               GpuArray<Real, AMREX_SPACEDIM> const& plo,
                               GpuArray<Real, AMREX_SPACEDIM> const& dx,
                               int type)
{
    const Real d = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

    const Real lx = (p.pos(0) - plo[0]) / dx[0]; // x
    const Real ly = (p.pos(1) - plo[1]) / dx[1]; // y
    const Real lz = (p.pos(2) - plo[2]) / dx[2]; // z

    const int i = static_cast<int>(Math::floor(lx)); // i
    const int j = static_cast<int>(Math::floor(ly)); // j
    const int k = static_cast<int>(Math::floor(lz)); // k

    Up = 0;
    Vp = 0;
    Wp = 0;
    //Euler to Lagrangian
    for(int ii = -2 + type; ii < 3 - type; ii++){
        for(int jj = -2 + type; jj < 3 - type; jj++){
            for(int kk = -2 + type; kk < 3 - type; kk ++){
                Real tU, tV, tW;
                const Real xi = plo[0] + (i + ii) * dx[0] + dx[0]/2;
                const Real yj = plo[1] + (j + jj) * dx[1] + dx[1]/2;
                const Real kz = plo[2] + (k + kk) * dx[2] + dx[2]/2;
                deltaFunction( p.pos(0), xi, dx[0], tU, type);
                deltaFunction( p.pos(1), yj, dx[1], tV, type);
                deltaFunction( p.pos(2), kz, dx[2], tW, type);
                const Real delta_value = tU * tV * tW;
                Up += delta_value * E(i + ii, j + jj, k + kk, EulerVIndex    ) * d;
                Vp += delta_value * E(i + ii, j + jj, k + kk, EulerVIndex + 1) * d;
                Wp += delta_value * E(i + ii, j + jj, k + kk, EulerVIndex + 2) * d;
            }
        }
    }
}

template<typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void VelocityInterpolationRKPM_cir(
    P p,
    ParticleReal& U,
    ParticleReal& V,
    ParticleReal& W,
    MAP_INFO const* rkpm_data,
    Array4<Real const> const& E,
    GpuArray<Real, AMREX_SPACEDIM> const& plo,
    GpuArray<Real, AMREX_SPACEDIM> const& dx,
    int EulerVIndex,
    int stencil_size)
{
    amrex::ignore_unused(p, plo, dx);

    U = 0;
    V = 0;
    W = 0;

    // Apply every stencil weight to the exact Euler cell it was generated for,
    // using each entry's own stored (i,j,k) index. This makes no assumption
    // about the stencil being an ordered 3x3x3 cube (or its center sitting at
    // entry 13) and is correct for any support shape / ordering produced by the
    // RKPM generator. Padding slots carry weight 0 with a valid index, so the
    // fixed-size loop is always memory-safe.
    for (int c = 0; c < stencil_size; ++c) {
        const auto& rkpm = rkpm_data[c];
        const int i = rkpm.index[0];
        const int j = rkpm.index[1];
        const int k = rkpm.index[2];
        U += rkpm.weight * rkpm.Vcell * E(i, j, k, EulerVIndex    );
        V += rkpm.weight * rkpm.Vcell * E(i, j, k, EulerVIndex + 1);
        W += rkpm.weight * rkpm.Vcell * E(i, j, k, EulerVIndex + 2);
    }
}

void mParticle::VelocityInterpolation(MultiFab &EulerVel,
                                      int type)//
{
    BL_PROFILE("mParticle::VelocityInterpolation");
    if (verbose) Print() << "\tmParticle::VelocityInterpolation\n";

    //Print() << "euler_finest_level " << euler_finest_level << std::endl;
    const auto& gm = mContainer->GetParGDB()->Geom(LOCAL_LEVEL);
    auto plo = gm.ProbLoArray();
    auto dx = gm.CellSizeArray();
    // attention
    // velocity ghost cells will be up-to-date
    EulerVel.FillBoundary(ParticleProperties::euler_velocity_index, 3, gm.periodicity());

    const int EulerVelocityIndex = ParticleProperties::euler_velocity_index;

    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){

        const Box& box = pti.validbox();

        auto& particles = pti.GetArrayOfStructs();
        const auto *p_ptr = particles.data();
        const Long np = pti.numParticles();
        const auto& ids = pti.GetIDs().data();

        auto& attri = pti.GetAttribs();
        auto* Up = attri[P_ATTR_REAL::U_Marker].data();
        auto* Vp = attri[P_ATTR_REAL::V_Marker].data();
        auto* Wp = attri[P_ATTR_REAL::W_Marker].data();
        const auto& E = EulerVel.array(pti);

        if (do_RKPM) {
            const MAP_INFO* rkpm_ptr = d_rkpm_flat.data();
            constexpr int STENCIL = mParticle::RKPM_STENCIL_SIZE;
            ParallelFor(np,
                [=] AMREX_GPU_DEVICE (const int i) {
                const auto id = p_ptr[i].id() - 1;
                VelocityInterpolationRKPM_cir(p_ptr[i], Up[i], Vp[i], Wp[i],
                                              rkpm_ptr + id * STENCIL, E, plo, dx, EulerVelocityIndex, STENCIL);
            });
        }else {
            ParallelFor(np,
                [=] AMREX_GPU_DEVICE (const int i) noexcept{
                VelocityInterpolation_cir(i, p_ptr[i], Up[i], Vp[i], Wp[i], E, EulerVelocityIndex, box.loVect(), box.hiVect(), plo, dx, type);
            });
        }
    }

    // if (verbose) mContainer->WriteAsciiFile(Concatenate("particle", 2));
}

void mParticle::ComputeLagrangianForce(Real dt)
{
    BL_PROFILE("mParticle::ComputeLagrangianForce");
    if (verbose) Print() << "\tmParticle::ComputeLagrangianForce\n";

    for(mParIter pti( *mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        auto& attri = pti.GetAttribs();
        auto const* p_ptr = pti.GetArrayOfStructs().data();

        const auto* Up = attri[P_ATTR_REAL::U_Marker].data();
        const auto* Vp = attri[P_ATTR_REAL::V_Marker].data();
        const auto* Wp = attri[P_ATTR_REAL::W_Marker].data();
        auto* FxP = attri[P_ATTR_REAL::Fx_Marker].data();
        auto* FyP = attri[P_ATTR_REAL::Fy_Marker].data();
        auto* FzP = attri[P_ATTR_REAL::Fz_Marker].data();

        const auto p_ids = pti.GetIDs().data();
        const auto ps = d_kernels.data();

        ParallelFor(np,
        [=] AMREX_GPU_DEVICE (const int i) noexcept{
            const auto p_id = p_ids[i];
            const auto& p = ps[p_id];
            const Real Ub = p.velocity[0];
            const Real Vb = p.velocity[1];
            const Real Wb = p.velocity[2];
            const Real Px = p.location[0];
            const Real Py = p.location[1];
            const Real Pz = p.location[2];
            Real target_u = Ub;
            Real target_v = Vb;
            Real target_w = Wb;

            if (is_prescribed_fin(p.geometry_type)) {
                const auto state = evaluate_fin_marker(p, p_ptr[i].id());
                target_u = state.u;
                target_v = state.v;
                target_w = state.w;
            } else {
                const auto rotational_velocity = p.omega.crossProduct(
                    RealVect(p_ptr[i].pos(0) - Px,
                             p_ptr[i].pos(1) - Py,
                             p_ptr[i].pos(2) - Pz));
                target_u += rotational_velocity[0];
                target_v += rotational_velocity[1];
                target_w += rotational_velocity[2];
            }

            FxP[i] = (target_u - Up[i]) / dt;
            FyP[i] = (target_v - Vp[i]) / dt;
            FzP[i] = (target_w - Wp[i]) / dt;
        });
    }
    // if (verbose) mContainer->WriteAsciiFile(Concatenate("particle", 3));
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ForceSpreading_cic (P const& p,
                         Real Px,
                         Real Py,
                         Real Pz,
                         ParticleReal& fxP,
                         ParticleReal& fyP,
                         ParticleReal& fzP,
                         ParticleReal& mxP,
                         ParticleReal& myP,
                         ParticleReal& mzP,
                         Array4<Real> const& E,
                         int EulerForceIndex,
                         Real dv,
                         GpuArray<Real,AMREX_SPACEDIM> const& plo,
                         GpuArray<Real,AMREX_SPACEDIM> const& dx,
                         int type)
{
    //const Real d = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    //plo to ii jj kk
    Real lx = (p.pos(0) - plo[0]) / dx[0];
    Real ly = (p.pos(1) - plo[1]) / dx[1];
    Real lz = (p.pos(2) - plo[2]) / dx[2];

    int i = static_cast<int>(Math::floor(lx));
    int j = static_cast<int>(Math::floor(ly));
    int k = static_cast<int>(Math::floor(lz));
    fxP *= dv;
    fyP *= dv;
    fzP *= dv;
    RealVect moment = RealVect((p.pos(0) - Px), (p.pos(1) - Py), (p.pos(2) - Pz)).crossProduct(RealVect(fxP, fyP, fzP));
    mxP = moment[0];
    myP = moment[1];
    mzP = moment[2];
    //lagrangian to Euler
    for(int ii = -2 + type; ii < +3 - type; ii++){
        for(int jj = -2 + type; jj < +3 - type; jj++){
            for(int kk = -2 + type; kk < +3 - type; kk ++){
                Real tU, tV, tW;
                const Real xi =plo[0] + (i + ii) * dx[0] + dx[0]/2;
                const Real yj =plo[1] + (j + jj) * dx[1] + dx[1]/2;
                const Real kz =plo[2] + (k + kk) * dx[2] + dx[2]/2;
                deltaFunction( p.pos(0), xi, dx[0], tU, type);
                deltaFunction( p.pos(1), yj, dx[1], tV, type);
                deltaFunction( p.pos(2), kz, dx[2], tW, type);
                Real delta_value = tU * tV * tW;
                HostDevice::Atomic::Add(&E(i + ii, j + jj, k + kk, EulerForceIndex  ), Real(delta_value * fxP));
                HostDevice::Atomic::Add(&E(i + ii, j + jj, k + kk, EulerForceIndex+1), Real(delta_value * fyP));
                HostDevice::Atomic::Add(&E(i + ii, j + jj, k + kk, EulerForceIndex+2), Real(delta_value * fzP));
            }
        }
    }
}

template<typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ForceSpreadingRKPM_cir(
    P p,
    Real Px,
    Real Py,
    Real Pz,
    ParticleReal& fxP,
    ParticleReal& fyP,
    ParticleReal& fzP,
    ParticleReal& mxP,
    ParticleReal& myP,
    ParticleReal& mzP,
    MAP_INFO const* rkpm_data,
    Real dv,
    Array4<Real> const &E,
    GpuArray<Real,AMREX_SPACEDIM> const& plo,
    GpuArray<Real,AMREX_SPACEDIM> const& dx,
    int EulerForceIndex,
    int stencil_size)
{
    amrex::ignore_unused(plo);
    const Real cellvol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

    const Real sx = fxP * dv;
    const Real sy = fyP * dv;
    const Real sz = fzP * dv;

    fxP = sx * cellvol;
    fyP = sy * cellvol;
    fzP = sz * cellvol;
    RealVect moment = RealVect(Real(p.pos(0) - Px), Real(p.pos(1) - Py), Real(p.pos(2) - Pz)).crossProduct(
                      RealVect(Real(fxP), Real(fyP), Real(fzP)));
    mxP = moment[0];
    myP = moment[1];
    mzP = moment[2];

    // Spread to the exact Euler cell each weight was generated for, using the
    // entry's own (i,j,k) index (mirror of VelocityInterpolationRKPM_cir, so the
    // interpolation and spreading operators stay exact transposes). Padding
    // slots carry weight 0 with a valid index, keeping the loop memory-safe.
    for (int c = 0; c < stencil_size; ++c) {
        const auto& rkpm = rkpm_data[c];
        const int i = rkpm.index[0];
        const int j = rkpm.index[1];
        const int k = rkpm.index[2];
        HostDevice::Atomic::Add(&E(i, j, k, EulerForceIndex    ), Real(rkpm.weight * rkpm.Vcell * sx));
        HostDevice::Atomic::Add(&E(i, j, k, EulerForceIndex + 1), Real(rkpm.weight * rkpm.Vcell * sy));
        HostDevice::Atomic::Add(&E(i, j, k, EulerForceIndex + 2), Real(rkpm.weight * rkpm.Vcell * sz));
    }
}


void mParticle::ForceSpreading(MultiFab & EulerForce,
                               int type)
{
    BL_PROFILE("mParticle::ForceSpreading");
    if (verbose) Print() << "\tmParticle::ForceSpreading\n";
    const auto& gm = mContainer->GetParGDB()->Geom(LOCAL_LEVEL);
    auto plo = gm.ProbLoArray();
    auto dxi = gm.CellSizeArray();
    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        const auto& particles = pti.GetArrayOfStructs();
        auto Uarray = EulerForce[pti].array();
        auto& attri = pti.GetAttribs();
        const auto& ids = pti.GetIDs().data();

        auto *const fxP_ptr = attri[P_ATTR_REAL::Fx_Marker].data();
        auto *const fyP_ptr = attri[P_ATTR_REAL::Fy_Marker].data();
        auto *const fzP_ptr = attri[P_ATTR_REAL::Fz_Marker].data();
        auto *const mxP_ptr = attri[P_ATTR_REAL::Mx_Marker].data();
        auto *const myP_ptr = attri[P_ATTR_REAL::My_Marker].data();
        auto *const mzP_ptr = attri[P_ATTR_REAL::Mz_Marker].data();
        const auto *const p_ptr = particles().data();

        auto force_index = ParticleProperties::euler_force_index;
        const auto ps = d_kernels.data();

        if (do_RKPM) {
            const MAP_INFO* rkpm_ptr = d_rkpm_flat.data();
            constexpr int STENCIL = mParticle::RKPM_STENCIL_SIZE;
            ParallelFor(np,
            [=] AMREX_GPU_DEVICE (const int i) noexcept{
                const auto p_id = p_ptr[i].id() - 1; // lagrangian marker's id
                const auto id = ids[i];  // particle's id
                auto loc_ptr = ps[id].location;
                auto dv = rkpm_ptr[p_id * STENCIL].eps;
                ForceSpreadingRKPM_cir(p_ptr[i], loc_ptr[0], loc_ptr[1], loc_ptr[2],
                                fxP_ptr[i], fyP_ptr[i], fzP_ptr[i],
                                mxP_ptr[i], myP_ptr[i], mzP_ptr[i],
                                rkpm_ptr + p_id * STENCIL, dv, Uarray, plo, dxi, force_index, STENCIL);
            });
        }else {
            ParallelFor(np,
            [=] AMREX_GPU_DEVICE (const int i) noexcept{
                const auto id = ids[i];
                auto loc_ptr = ps[id].location;
                auto dv = ps[id].dv;
                if (is_prescribed_fin(ps[id].geometry_type)) {
                    const auto state = evaluate_fin_marker(ps[id], p_ptr[i].id());
                    dv *= state.quadrature_scale;
                }
                ForceSpreading_cic(p_ptr[i], loc_ptr[0], loc_ptr[1], loc_ptr[2],
                                   fxP_ptr[i], fyP_ptr[i], fzP_ptr[i],
                                   mxP_ptr[i], myP_ptr[i], mzP_ptr[i],
                                   Uarray, force_index, dv, plo, dxi, type);
            });
        }
    }
    using pc = mParticleContainer::SuperParticleType;
    // particle id => thread id
    BL_PROFILE_VAR("ForceSpreading::reduce_per_particle", blp_fsreduce);

    if (ParticleProperties::reduce_method == 1) {
        static bool printed_once = false;
        if (!printed_once) {
            printed_once = true;
            amrex::Print() << "[Particle] : reduce_method = 1 : optimized ParticleReduce (6-tuple) + packed all-reduce in use\n";
        }
        // ===== optimized: one 6-tuple ParticleReduce per kernel (1 scan per
        // ===== particle instead of 6), then one packed MPI all-reduce (6N -> 1).
        // Note: ParticleReduce returns per-rank partial sums only, so the packed
        // all-reduce below is required. See https://github.com/AMReX-Codes/amrex/discussions/4593
        const int nkern = particle_kernels.size();
        Vector<Real> ib_fm(6 * nkern, 0.0);
        // ===== optimized v2: single pass over all local markers with atomicAdd
        // ===== into a 6*nkern device accumulator (one scan instead of nkern
        // ===== full ParticleReduce scans), then one packed MPI all-reduce.
        amrex::Gpu::DeviceVector<Real> d_ib_fm(6 * nkern, 0.0);
        auto* d_ptr = d_ib_fm.data();

        for (mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti) {
            // M_ID lives in the SoA int data (the AoS type here is Particle<0,0>
            // and has no idata), so read the kernel id through GetIDs() like the
            // spreading loop above does.
            const auto *const ids = pti.GetIDs().data();
            auto& attri = pti.GetAttribs();
            const Long np = pti.numParticles();
            const auto *const fxP_ptr = attri[P_ATTR_REAL::Fx_Marker].data();
            const auto *const fyP_ptr = attri[P_ATTR_REAL::Fy_Marker].data();
            const auto *const fzP_ptr = attri[P_ATTR_REAL::Fz_Marker].data();
            const auto *const mxP_ptr = attri[P_ATTR_REAL::Mx_Marker].data();
            const auto *const myP_ptr = attri[P_ATTR_REAL::My_Marker].data();
            const auto *const mzP_ptr = attri[P_ATTR_REAL::Mz_Marker].data();

            ParallelFor(np,
                [=] AMREX_GPU_DEVICE (const int i) noexcept {
                    // M_ID equals the kernel index (0..nkern-1); Particle_id[M_ID] = cur_p.id
                    const int pid = ids[i];
                    const Long base = 6 * pid;
                    amrex::Gpu::Atomic::Add(d_ptr + base + 0, fxP_ptr[i]);
                    amrex::Gpu::Atomic::Add(d_ptr + base + 1, fyP_ptr[i]);
                    amrex::Gpu::Atomic::Add(d_ptr + base + 2, fzP_ptr[i]);
                    amrex::Gpu::Atomic::Add(d_ptr + base + 3, mxP_ptr[i]);
                    amrex::Gpu::Atomic::Add(d_ptr + base + 4, myP_ptr[i]);
                    amrex::Gpu::Atomic::Add(d_ptr + base + 5, mzP_ptr[i]);
                }
            );
        }
        Gpu::copyAsync(Gpu::deviceToHost, d_ib_fm.begin(), d_ib_fm.end(), ib_fm.begin());
        Gpu::streamSynchronize();
        ParallelAllReduce::Sum(ib_fm.data(), ib_fm.size(), ParallelDescriptor::Communicator());
        for (int k = 0; k < nkern; ++k) {
            particle_kernels[k].ib_force  += {ib_fm[6*k+0], ib_fm[6*k+1], ib_fm[6*k+2]};
            particle_kernels[k].ib_moment += {ib_fm[6*k+3], ib_fm[6*k+4], ib_fm[6*k+5]};
        }
    }
    else
    {
        static bool printed_once = false;
        if (!printed_once) {
            printed_once = true;
            amrex::Print() << "[Particle] : reduce_method = 0 : original ReduceSum loop in use\n";
        }
        // ===== original: 6 ReduceSum + 6 ParallelAllReduce per kernel =====
        for (auto& cur_p : particle_kernels) { // gm position
            auto fx = ReduceSum(*mContainer, [=] AMREX_GPU_HOST_DEVICE(const pc& p) -> ParticleReal{
                if (p.idata(M_ID) == cur_p.id) {
                    return p.rdata(P_ATTR_REAL::Fx_Marker);
                }
                return 0.;
            });
            auto fy = ReduceSum(*mContainer, [=] AMREX_GPU_HOST_DEVICE(const pc& p) -> ParticleReal{
                if (p.idata(M_ID) == cur_p.id) {
                    return p.rdata(P_ATTR_REAL::Fy_Marker);
                }
                return 0.;
            });
            auto fz = ReduceSum(*mContainer, [=] AMREX_GPU_HOST_DEVICE(const pc& p) -> ParticleReal{
                if (p.idata(M_ID) == cur_p.id) {
                    return p.rdata(P_ATTR_REAL::Fz_Marker);
                }
                return 0.;
            });
            auto mx = ReduceSum(*mContainer, [=] AMREX_GPU_HOST_DEVICE(const pc& p) -> ParticleReal{
                if (p.idata(M_ID) == cur_p.id) {
                    return p.rdata(P_ATTR_REAL::Mx_Marker);
                }
                return 0.;
            });
            auto my = ReduceSum(*mContainer, [=] AMREX_GPU_HOST_DEVICE(const pc& p) -> ParticleReal{
                if (p.idata(M_ID) == cur_p.id) {
                    return p.rdata(P_ATTR_REAL::My_Marker);
                }
                return 0.;
            });
            auto mz = ReduceSum(*mContainer, [=] AMREX_GPU_HOST_DEVICE(const pc& p) -> ParticleReal{
                if (p.idata(M_ID) == cur_p.id) {
                    return p.rdata(P_ATTR_REAL::Mz_Marker);
                }
                return 0.;
            });

            ParallelAllReduce::Sum(fx, ParallelDescriptor::Communicator());
            ParallelAllReduce::Sum(fy, ParallelDescriptor::Communicator());
            ParallelAllReduce::Sum(fz, ParallelDescriptor::Communicator());
            ParallelAllReduce::Sum(mx, ParallelDescriptor::Communicator());
            ParallelAllReduce::Sum(my, ParallelDescriptor::Communicator());
            ParallelAllReduce::Sum(mz, ParallelDescriptor::Communicator());

            cur_p.ib_force += {Real(fx), Real(fy), Real(fz)};
            cur_p.ib_moment += {Real(mx), Real(my), Real(mz)};
        }
    }
    BL_PROFILE_VAR_STOP(blp_fsreduce);
    EulerForce.SumBoundary(ParticleProperties::euler_force_index, 3, gm.periodicity());
}

void mParticle::VelocityCorrection(MultiFab &Euler, MultiFab &EulerForce, Real dt) const
{
    BL_PROFILE("mParticle::VelocityCorrection");
    if(verbose) Print() << "\tmParticle::VelocityCorrection\n";
    MultiFab::Saxpy(Euler, dt, EulerForce, ParticleProperties::euler_force_index, ParticleProperties::euler_velocity_index, 3, 0); //VelocityCorrection
}

// pvf_method = 0 : reference implementation. For every particle the nodal level
// set and the pvf are evaluated on the whole domain, the momentum sums are
// reduced over the whole domain, and one set of MPI collectives is issued per
// particle. Cost per step ~ N_particle x domain.
void mParticle::UpdateParticlesReference(int iStep,
                                         const MultiFab& Euler_old,
                                         const MultiFab& Euler,
                                         MultiFab& phi_nodal,
                                         MultiFab& pvf,
                                         Real dt)
{
    MultiFab AllParticlePVF(pvf.boxArray(), pvf.DistributionMap(), pvf.nComp(), pvf.nGrow());
    AllParticlePVF.setVal(0.0);
    BL_PROFILE_VAR("UpdateParticles::phi_pvf_per_particle", blp_phipvf);
    //continue condition 6DOF
    for(auto& kernel : particle_kernels){

        if (is_prescribed_fin(kernel.geometry_type)) {
            // A prescribed fin is an open surface. It has no closed particle
            // volume and its analytical motion must not enter rigid-body 6DOF.
            continue;
        }

        calculate_phi_nodal(phi_nodal, kernel);
        nodal_phi_to_pvf(pvf, phi_nodal);

        int ncomp = pvf.nComp();
        int ngrow = pvf.nGrow();
        MultiFab pvf_old(pvf.boxArray(), pvf.DistributionMap(), ncomp, ngrow);
        MultiFab::Copy(pvf_old, pvf, 0, 0, ncomp, ngrow);

        int loop = ParticleProperties::loop_solid;

        while (loop > 0 && iStep > ParticleProperties::start_step) {

                kernel.sum_u_new.scale(0.0);
                kernel.sum_u_old.scale(0.0);
                // sum U
                CalculateSumU_cir(kernel.sum_u_new, Euler, pvf, ParticleProperties::euler_velocity_index);
                CalculateSumU_cir(kernel.sum_u_old, Euler_old, pvf_old, ParticleProperties::euler_velocity_index);
                ParallelAllReduce::Sum(kernel.sum_u_new.dataPtr(), 3, ParallelDescriptor::Communicator());
                ParallelAllReduce::Sum(kernel.sum_u_old.dataPtr(), 3, ParallelDescriptor::Communicator());

                kernel.sum_t_new.scale(0.0);
                kernel.sum_t_old.scale(0.0);
                // sum T
                CalculateSumT_cir(kernel.sum_t_new, Euler, pvf, kernel.location, ParticleProperties::euler_velocity_index);
                CalculateSumT_cir(kernel.sum_t_old, Euler_old, pvf_old, kernel.location, ParticleProperties::euler_velocity_index);
                ParallelAllReduce::Sum(kernel.sum_t_new.dataPtr(), 3, ParallelDescriptor::Communicator());
                ParallelAllReduce::Sum(kernel.sum_t_old.dataPtr(), 3, ParallelDescriptor::Communicator());

            // 6DOF (I/O rank only), then broadcast the new state
            SixDOFUpdate(kernel, dt);
            ParallelDescriptor::Bcast(&kernel.location[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.location_old[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.velocity[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.velocity_old[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.omega[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.omega_old[0],3,ParallelDescriptor::IOProcessorNumber());

            loop--;

            if (loop > 0) {
                calculate_phi_nodal(phi_nodal, kernel);
                nodal_phi_to_pvf(pvf, phi_nodal);
            }

        }

        RecordOldValue(kernel);
        MultiFab::Add(AllParticlePVF, pvf, 0, 0, 1, 0); // do not copy ghost cell values
    }
    BL_PROFILE_VAR_STOP(blp_phipvf);
    // calculate the pvf based on the information of all particles
    MultiFab::Copy(pvf, AllParticlePVF, 0, 0, 1, pvf.nGrow());
}

// 6DOF update of one particle from its momentum sums, IB force/moment and
// collision force. Only the I/O rank computes; callers broadcast afterwards.
void mParticle::SixDOFUpdate(kernel& kernel, Real dt)
{
    if (is_prescribed_fin(kernel.geometry_type)) return;
    if (ParallelDescriptor::MyProc() != ParallelDescriptor::IOProcessorNumber()) return;

    for(auto idir : {0,1,2})
    {
        //TL
        if (kernel.TL[idir] == 0) {
            kernel.velocity[idir] = 0.0;
        }
        else if (kernel.TL[idir] == 1) {
            kernel.location[idir] = kernel.location_old[idir] + (kernel.velocity[idir] + kernel.velocity_old[idir]) * dt * 0.5;
        }
        else if (kernel.TL[idir] == 2) {
            if(!ParticleProperties::Uhlmann){
                kernel.velocity[idir] = kernel.velocity_old[idir]
                                    + ((kernel.sum_u_new[idir] - kernel.sum_u_old[idir]) * ParticleProperties::euler_fluid_rho / dt
                                    - kernel.ib_force[idir] * ParticleProperties::euler_fluid_rho
                                    + m_gravity[idir] * (kernel.rho - ParticleProperties::euler_fluid_rho) * kernel.Vp
                                    + kernel.Fcp[idir]) * dt / kernel.rho / kernel.Vp ;
            }else{
                //Uhlmann
                kernel.velocity[idir] = kernel.velocity_old[idir]
                                    + (ParticleProperties::euler_fluid_rho / kernel.Vp /(ParticleProperties::euler_fluid_rho - kernel.rho)*kernel.ib_force[idir]
                                    + m_gravity[idir]) * dt;
            }
            kernel.location[idir] = kernel.location_old[idir] + (kernel.velocity[idir] + kernel.velocity_old[idir]) * dt * 0.5;
        }
        else {
            Print() << "Particle (" << kernel.id << ") has wrong TL"<< direction_str[idir] <<" value\n";
            Abort("Stop here!");
        }
        //RL
        if (kernel.RL[idir] == 0) {
            kernel.omega[idir] = 0.0;
        }
        else if (kernel.RL[idir] == 1) {
        }
        else if (kernel.RL[idir] == 2) {
            if(!ParticleProperties::Uhlmann){
                kernel.omega[idir] = kernel.omega_old[idir]
                                + ((kernel.sum_t_new[idir] - kernel.sum_t_old[idir]) * ParticleProperties::euler_fluid_rho / dt
                                - kernel.ib_moment[idir] * ParticleProperties::euler_fluid_rho
                                + kernel.Tcp[idir]) * dt / cal_momentum(kernel.rho, kernel.radius, kernel.geometry_type, idir, kernel.radius2, kernel.radius3);
            }else{
                //Uhlmann
                kernel.omega[idir] = kernel.omega_old[idir]
                                + ParticleProperties::euler_fluid_rho /(ParticleProperties::euler_fluid_rho - kernel.rho) * kernel.ib_moment[idir] * kernel.dv
                                / cal_momentum(kernel.rho, kernel.radius, kernel.geometry_type, idir, kernel.radius2, kernel.radius3) * kernel.rho * dt;
            }
        }
        else {
            Print() << "Particle (" << kernel.id << ") has wrong RL"<< direction_str[idir] <<" value\n";
            Abort("Stop here!");
        }
    }
}

// Broadcast location/velocity/omega (+ their *_old) of all particles from the
// I/O rank in one message instead of 6 broadcasts per particle.
void mParticle::BroadcastKernelState()
{
    const int nk = static_cast<int>(particle_kernels.size());
    constexpr int NV = 18;
    Vector<Real> buf(NV * nk);
    if (ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()) {
        for (int k = 0; k < nk; ++k) {
            auto const& p = particle_kernels[k];
            Real* b = buf.data() + NV * k;
            for (int d = 0; d < 3; ++d) {
                b[d]      = p.location[d];
                b[3 + d]  = p.location_old[d];
                b[6 + d]  = p.velocity[d];
                b[9 + d]  = p.velocity_old[d];
                b[12 + d] = p.omega[d];
                b[15 + d] = p.omega_old[d];
            }
        }
    }
    ParallelDescriptor::Bcast(buf.data(), buf.size(), ParallelDescriptor::IOProcessorNumber());
    for (int k = 0; k < nk; ++k) {
        auto& p = particle_kernels[k];
        const Real* b = buf.data() + NV * k;
        for (int d = 0; d < 3; ++d) {
            p.location[d]     = b[d];
            p.location_old[d] = b[3 + d];
            p.velocity[d]     = b[6 + d];
            p.velocity_old[d] = b[9 + d];
            p.omega[d]        = b[12 + d];
            p.omega_old[d]    = b[15 + d];
        }
    }
}

// (particle, local grid, intersection box) list for the pvf kernels: every
// cell of grid `grid` in `ibx` may carry a non-zero pvf of particle `kidx`
// evaluated at loc[kidx] (and loc_old[kidx] when that differs). Cached and
// reused as long as the centres and the grids are unchanged.
const Vector<mParticle::PvfWork>&
mParticle::PvfWorkList(const BoxArray& ba, const DistributionMapping& dm,
                       const Vector<RealVect>& loc, const Vector<RealVect>& loc_old)
{
    if (loc == m_pvf_work_loc && loc_old == m_pvf_work_loc_old && ba == m_pvf_work_ba) {
        return m_pvf_work;
    }
    m_pvf_work.clear();
    const int myproc = ParallelDescriptor::MyProc();
    const int nk = static_cast<int>(loc.size());
    for (int k = 0; k < nk; ++k) {
        auto const& kern = particle_kernels[k];
        // Open fin surfaces have no closed-body volume fraction.
        if (is_prescribed_fin(kern.geometry_type)) continue;
        const Real R = (kern.geometry_type == 2)
                     ? amrex::max(kern.radius, amrex::max(kern.radius2, kern.radius3)) : kern.radius;
        Box bbox = particle_index_box(loc[k], R);
        if (loc_old[k] != loc[k]) bbox.minBox(particle_index_box(loc_old[k], R));
        for (auto const& is : ba.intersections(bbox)) {
            if (dm[is.first] == myproc) m_pvf_work.push_back({k, is.first, is.second});
        }
    }
    m_pvf_work_loc = loc;
    m_pvf_work_loc_old = loc_old;
    m_pvf_work_ba = ba;
    return m_pvf_work;
}

// pvf_method = 1 : same numerics as UpdateParticlesReference(), but
//  * every particle only touches the cells of its own bounding box (pvf is
//    exactly zero elsewhere), phi_nodal / pvf_old are never materialised, and
//    the pvf evaluation and the momentum sums are fused into one small kernel
//    per (particle, box) pair;
//  * all MPI traffic of one 6DOF iteration is one packed all-reduce (12 sums
//    per particle) plus one packed broadcast of the particle state.
// The per-step cost is therefore ~ (local particles) x (cells per particle)
// instead of N_particle x domain, which is what weak scaling needs.
void mParticle::UpdateParticlesFused(int iStep,
                                     const MultiFab& Euler_old,
                                     const MultiFab& Euler,
                                     MultiFab& pvf,
                                     Real dt)
{
    const int nk = static_cast<int>(particle_kernels.size());
    const auto dx  = ParticleProperties::dx;
    const auto plo = ParticleProperties::plo;
    const Real dvol = Math::powi<3>(dx[0]);
    const int vidx = ParticleProperties::euler_velocity_index;

    for (auto const& kern : particle_kernels) {
        if (kern.geometry_type > 2 && !is_prescribed_fin(kern.geometry_type)) {
            Print() << "Particle (" << kern.id << ") has unsupported geometry_type: " << kern.geometry_type << "\n";
            Abort("Unsupported geometry type. Only geometry_type = 1 (sphere) and 2 (ellipsoid) are supported.");
        }
    }

    // The reference path copies pvf into pvf_old right after entering the
    // function, so the "old" sums always use the pvf at the entry location.
    // loc_pvf tracks where the most recent pvf was evaluated; the final
    // accumulation into pvf uses it (see the reference path).
    Vector<RealVect> loc_entry(nk), loc_pvf(nk);
    for (int k = 0; k < nk; ++k) {
        loc_entry[k] = particle_kernels[k].location;
        loc_pvf[k]   = particle_kernels[k].location;
    }

    // 12 partial sums per particle: sum_u_new, sum_u_old, sum_t_new, sum_t_old
    constexpr int NS = 12;
    Gpu::DeviceVector<Real> d_sums(NS * nk);
    Vector<Real> h_sums(NS * nk);
    Real* const d_sums_ptr = d_sums.data();

    BL_PROFILE_VAR("UpdateParticles::phi_pvf_per_particle", blp_phipvf);
    int loop = ParticleProperties::loop_solid;
    while (loop > 0 && iStep > ParticleProperties::start_step) {

        // ---- local sums on the GPU: one small kernel per (particle, box) ----
        BL_PROFILE_VAR("UpdateParticles::pvf_sums_local", blp_local);
        ParallelFor(NS * nk, [=] AMREX_GPU_DEVICE (int i) noexcept { d_sums_ptr[i] = 0.0; });

        for (int kidx = 0; kidx < nk; ++kidx) loc_pvf[kidx] = particle_kernels[kidx].location;
        // (particle, local grid, intersection box) triples of this rank only:
        // no MFIter per particle, and nothing at all for particles that live
        // on other ranks, so the host cost is ~ local particles, not N_total.
        auto const& work = PvfWorkList(Euler.boxArray(), Euler.DistributionMap(), loc_pvf, loc_entry);
        for (auto const& w : work) {
            const int kidx = w.kidx;
            auto const& kern = particle_kernels[kidx];
            const RealVect loc_c = kern.location;   // current pvf and torque arm
            const RealVect loc_e = loc_entry[kidx]; // "old" pvf
            const bool same_loc = (loc_c == loc_e);
            const Real a = kern.radius;
            const Real b = kern.radius2;
            const Real c = kern.radius3;
            const int  gtype = kern.geometry_type;
            Real* const sums = d_sums_ptr + NS * kidx;
            {
                const Box ibx = w.ibx;
                auto const& En = Euler.const_array(w.grid);
                auto const& Eo = Euler_old.const_array(w.grid);
                ParallelFor(ibx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const Real pvf_c = cell_pvf_at(i, j, k, dx, plo, loc_c[0], loc_c[1], loc_c[2], a, b, c, gtype);
                    const Real pvf_e = same_loc ? pvf_c
                                     : cell_pvf_at(i, j, k, dx, plo, loc_e[0], loc_e[1], loc_e[2], a, b, c, gtype);
                    if (pvf_c == 0.0 && pvf_e == 0.0) return;

                    // same expressions as CalculateSumU_cir / CalculateSumT_cir
                    Real x = plo[0] + i*dx[0] + 0.5*dx[0];
                    Real y = plo[1] + j*dx[1] + 0.5*dx[1];
                    Real z = plo[2] + k*dx[2] + 0.5*dx[2];
                    const RealVect r(x - loc_c[0], y - loc_c[1], z - loc_c[2]);
                    const RealVect vn(En(i,j,k,vidx), En(i,j,k,vidx+1), En(i,j,k,vidx+2));
                    const RealVect vo(Eo(i,j,k,vidx), Eo(i,j,k,vidx+1), Eo(i,j,k,vidx+2));
                    const RealVect tn = r.crossProduct(vn);
                    const RealVect to = r.crossProduct(vo);

                    Gpu::Atomic::AddNoRet(sums + 0,  vn[0] * dvol * pvf_c);
                    Gpu::Atomic::AddNoRet(sums + 1,  vn[1] * dvol * pvf_c);
                    Gpu::Atomic::AddNoRet(sums + 2,  vn[2] * dvol * pvf_c);
                    Gpu::Atomic::AddNoRet(sums + 3,  vo[0] * dvol * pvf_e);
                    Gpu::Atomic::AddNoRet(sums + 4,  vo[1] * dvol * pvf_e);
                    Gpu::Atomic::AddNoRet(sums + 5,  vo[2] * dvol * pvf_e);
                    Gpu::Atomic::AddNoRet(sums + 6,  tn[0] * dvol * pvf_c);
                    Gpu::Atomic::AddNoRet(sums + 7,  tn[1] * dvol * pvf_c);
                    Gpu::Atomic::AddNoRet(sums + 8,  tn[2] * dvol * pvf_c);
                    Gpu::Atomic::AddNoRet(sums + 9,  to[0] * dvol * pvf_e);
                    Gpu::Atomic::AddNoRet(sums + 10, to[1] * dvol * pvf_e);
                    Gpu::Atomic::AddNoRet(sums + 11, to[2] * dvol * pvf_e);
                });
            }
        }
        BL_PROFILE_VAR_STOP(blp_local);

        // ---- one packed all-reduce, 6DOF on the I/O rank, one packed broadcast ----
        BL_PROFILE_VAR("UpdateParticles::6dof_comm", blp_comm);
        Gpu::copyAsync(Gpu::deviceToHost, d_sums.begin(), d_sums.end(), h_sums.begin());
        Gpu::streamSynchronize();
        ParallelAllReduce::Sum(h_sums.data(), NS * nk, ParallelDescriptor::Communicator());
        for (int kidx = 0; kidx < nk; ++kidx) {
            auto& kern = particle_kernels[kidx];
            if (is_prescribed_fin(kern.geometry_type)) continue;
            const Real* s = h_sums.data() + NS * kidx;
            kern.sum_u_new = RealVect(s[0], s[1],  s[2]);
            kern.sum_u_old = RealVect(s[3], s[4],  s[5]);
            kern.sum_t_new = RealVect(s[6], s[7],  s[8]);
            kern.sum_t_old = RealVect(s[9], s[10], s[11]);
            SixDOFUpdate(kern, dt);
        }
        BroadcastKernelState();
        BL_PROFILE_VAR_STOP(blp_comm);

        loop--;
    }

    for (auto& kern : particle_kernels) {
        if (!is_prescribed_fin(kern.geometry_type)) RecordOldValue(kern);
    }

    // ---- pvf of all particles: accumulate in particle order, valid cells only,
    //      ghost cells stay zero (as in the reference path) ----
    BL_PROFILE_VAR("UpdateParticles::pvf_accumulate", blp_acc);
    pvf.setVal(0.0);
    auto const& work_acc = PvfWorkList(pvf.boxArray(), pvf.DistributionMap(), loc_pvf, loc_pvf);
    for (auto const& w : work_acc) {
        auto const& kern = particle_kernels[w.kidx];
        const RealVect loc = loc_pvf[w.kidx];
        const Real a = kern.radius;
        const Real b = kern.radius2;
        const Real c = kern.radius3;
        const int  gtype = kern.geometry_type;
        auto const& pv = pvf.array(w.grid);
        ParallelFor(w.ibx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            pv(i,j,k) += cell_pvf_at(i, j, k, dx, plo, loc[0], loc[1], loc[2], a, b, c, gtype);
        });
    }
    BL_PROFILE_VAR_STOP(blp_acc);
    BL_PROFILE_VAR_STOP(blp_phipvf);
}

void mParticle::UpdateParticles(int iStep,
                                Real time,
                                const MultiFab& Euler_old,
                                const MultiFab& Euler,
                                MultiFab& phi_nodal,
                                MultiFab& pvf,
                                Real dt)
{
    BL_PROFILE("mParticle::UpdateParticles");
    if (verbose) Print() << "mParticle::UpdateParticles\n";
    ib_time = time + dt;
    // start record
    auto UpdateParticlesStart = ParallelDescriptor::second();

    //Particle Collision calculation
    BL_PROFILE_VAR("UpdateParticles::collision", blp_coll);
    DoParticleCollision(ParticleProperties::collision_model);
    BL_PROFILE_VAR_STOP(blp_coll);

    static bool printed_once = false;
    if (!printed_once) {
        printed_once = true;
        amrex::Print() << "[Particle] : pvf_method = " << ParticleProperties::pvf_method
                       << (ParticleProperties::pvf_method == 1
                           ? " : bounding-box fused pvf/sums + packed MPI in use\n"
                           : " : reference per-particle whole-domain pvf in use\n");
    }
    if (ParticleProperties::pvf_method == 1) {
        UpdateParticlesFused(iStep, Euler_old, Euler, pvf, dt);
    } else {
        UpdateParticlesReference(iStep, Euler_old, Euler, phi_nodal, pvf, dt);
    }

    BL_PROFILE_VAR("UpdateParticles::print", blp_print);
    spend_time += ParallelDescriptor::second() - UpdateParticlesStart;
    ParallelDescriptor::ReduceRealMax(spend_time);
    Print() << "[DIBM] IB and update particle, step : "<< iStep <<", time : " << spend_time << "\n";
    BL_PROFILE_VAR_STOP(blp_print);

    int particle_write_freq = ParticleProperties::write_freq;
    if (iStep % particle_write_freq == 0) {
        BL_PROFILE_VAR("UpdateParticles::write_csv", blp_csv);
        // Every rank holds the complete particle state (all of it is either
        // all-reduced or broadcast above), so spread the per-particle CSV
        // files over the ranks round-robin. Creating a file on the parallel
        // file system costs ~10-20 ms; doing all N of them on the I/O rank
        // serialises N x that time and makes every other rank wait at the
        // next collective.
        const int nprocs = ParallelDescriptor::NProcs();
        const int myproc = ParallelDescriptor::MyProc();
        const int nk = static_cast<int>(particle_kernels.size());
        for (int k = myproc; k < nk; k += nprocs) {
            const Real output_time = is_prescribed_fin(particle_kernels[k].geometry_type)
                                   ? ib_time : time;
            WriteIBForceAndMoment(iStep, output_time, dt, particle_kernels[k]);
        }
        BL_PROFILE_VAR_STOP(blp_csv);
    }

    // if (verbose) mContainer->WriteAsciiFile(Concatenate("particle", 4));
}

void mParticle::DoParticleCollision(int model)
{
    Vector<kernel*> rigid_bodies;
    for (auto& body : particle_kernels) {
        if (!is_prescribed_fin(body.geometry_type)) {
            rigid_bodies.push_back(&body);
        }
    }
    if(rigid_bodies.size() < 2 ) return ;

    if (verbose) Print() << "\tmParticle::DoParticleCollision\n";

    if(ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()){
        for(const auto* body : rigid_bodies){
            m_Collision.InsertParticle(body->location, body->velocity,
                                       body->radius, body->rho);
        }

        m_Collision.takeModel(model);

        for(auto* body : rigid_bodies){
            body->Fcp = m_Collision.Particles.front().preForece
                      * body->Vp * body->rho * m_gravity.vectorLength();
            m_Collision.Particles.pop_front();
        }
    }
    // one packed broadcast instead of one per particle
    const int nk = static_cast<int>(particle_kernels.size());
    Vector<Real> buf(3 * nk);
    for (int k = 0; k < nk; ++k) {
        for (int d = 0; d < 3; ++d) buf[3*k + d] = particle_kernels[k].Fcp[d];
    }
    ParallelDescriptor::Bcast(buf.data(), buf.size(), ParallelDescriptor::IOProcessorNumber());
    for (int k = 0; k < nk; ++k) {
        for (int d = 0; d < 3; ++d) particle_kernels[k].Fcp[d] = buf[3*k + d];
    }
}

void mParticle::RecordOldValue(kernel& kernel)
{
    kernel.location_old = kernel.location;
    kernel.velocity_old = kernel.velocity;
    kernel.omega_old = kernel.omega;
}

void mParticle::ResolveLagrangianMarker(std::string marker_file) {
    if (verbose) Print() << "\tmParticle::ResolveLagrangianMarker\n";
    // resolve point dict to LargrangianMarker
    /** point dict
    {
        0: (1.0565757615805609, 0.8544864844710995, 1.475),
        ...
    }
    {
        10: (0.7637951673013351, 1.1165430264314247, 1.425),
        ...
    }
    **/

    // txt file resolve
    std::ifstream marker(marker_file);
    std::string line;
    size_t number;
    size_t begin{0};
    while (std::getline(marker, line)) {
        // skip
        if (line.empty() || line.find('{') != std::string::npos ) {
            StartOfMarker.push_back(begin);
            number = 0;
            continue;
        }
        if (line.find('}') != std::string::npos) {
            NumOfMarker.push_back(number);
            continue;
        }
        // clear blank char
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);
        // sub location
        size_t colon_pos = line.find(':');
        size_t start = line.find('(');
        size_t end = line.find(')');
        std::string locations = line.substr(start + 1, end - start - 1);
        // location
        std::vector<Real> location;
        std::stringstream vs(locations);
        std::string v;
        while (getline(vs, v, ',')) {
            location.push_back(std::stod(v));
        }
        number++;
        begin++;
        LargrangianMarker.push_back(RealVect{location});
    }
    Print() << "\tLargrangianMarker size : " << LargrangianMarker.size()  << "\n";
}

void mParticle::ResolveWithRPKM(std::string RKPM_file) {
    if (verbose) Print() << "\tmParticle::ResolveWithRPKM\n";
    // resolve RKPM file
    /**
    0: [
        {"i": 4, "j": 3, "k": 6, "w": -0.309, "Vcell": 0.008, "eps":1.017},
        {"i": 4, "j": 3, "k": 7, "w": -1.365, "Vcell": 0.008, "eps":1.017},
        {"i": 4, "j": 3, "k": 8, "w": -1.327, "Vcell": 0.008, "eps":1.017},
        ... # 平均27个欧拉点
    ],
    1: [...],
    ...
    **/
    // txt file
    std::ifstream RKPM(RKPM_file);
    std::string line, dict_arr;
    int key_id;
    bool in_dict;
    while (getline(RKPM, line)) {
        // clear blank char
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        // clear empty and #
        if (line.empty() || line.find("#") == 0) continue;
        size_t colonPos = line.find(':');
        // dict start
        if (colonPos != std::string::npos && line.find('[') != std::string::npos) {
            key_id = std::stoi(line.substr(0, colonPos));
            in_dict = true;
            dict_arr = "";

            continue;
        }
        // get all {}
        if (in_dict) {
            dict_arr += line;
        }
        // dict end
        if (in_dict && line.find(']') != std::string::npos) {
            in_dict = false;
            Vector<MAP_INFO> r;
            int start = 0;
            while ((start = dict_arr.find('{', start)) != std::string::npos) {
                size_t end = dict_arr.find('}', start);
                std::string str = dict_arr.substr(start, end - start + 1);
                // resolve {}
                str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());
                int s = 0;
                MAP_INFO t;
                while ((s = str.find('"', s)) != std::string::npos) {
                    size_t e = str.find('"', s + 1);
                    std::string key = str.substr(s + 1, e - s - 1);

                    s = str.find(':', e) + 1;
                    size_t valueEnd = str.find_first_of(",}", s);
                    std::string valueStr = str.substr(s, valueEnd - s);

                    if (key == "i") t.index[0] = stoi(valueStr);
                    else if (key == "j") t.index[1] = stoi(valueStr);
                    else if (key == "k") t.index[2] = stoi(valueStr);
                    else if (key == "w") t.weight = stod(valueStr);
                    else if (key == "Vcell") t.Vcell = stod(valueStr);
                    else if (key == "eps") t.eps = stod(valueStr);

                    s = valueEnd + 1;
                }
                r.push_back(t);
                start = end + 1;
            }
            RKPM_MAP[key_id] = r;
            continue;
        }
    }
    Print() << "\tRKPM mapping size : " << RKPM_MAP.size() << "\n";

    // Flatten RKPM_MAP into a contiguous device-accessible array for GPU ParallelFor
    int max_key = 0;
    for (const auto& [key, vec] : RKPM_MAP) {
        max_key = std::max(max_key, key);
        AMREX_ALWAYS_ASSERT(int(vec.size()) <= RKPM_STENCIL_SIZE);
    }
    h_rkpm_flat.resize((max_key + 1) * RKPM_STENCIL_SIZE);
    for (const auto& [key, vec] : RKPM_MAP) {
        const int n = int(vec.size());
        for (int j = 0; j < n; ++j) {
            h_rkpm_flat[key * RKPM_STENCIL_SIZE + j] = vec[j];
        }
        // Pad unused slots with a zero-weight entry that carries a VALID cell
        // index (copied from the first real entry). The interpolation/spreading
        // loops run over the fixed RKPM_STENCIL_SIZE and dereference each entry's
        // index; without a valid padding index they would read cell (0,0,0),
        // which is out of range for the particle's box.
        MAP_INFO pad = (n > 0) ? vec[0] : MAP_INFO{};
        pad.weight = 0.0;
        pad.Vcell  = 0.0;
        for (int j = n; j < RKPM_STENCIL_SIZE; ++j) {
            h_rkpm_flat[key * RKPM_STENCIL_SIZE + j] = pad;
        }
    }
    d_rkpm_flat.resize(h_rkpm_flat.size());
    Gpu::copyAsync(Gpu::hostToDevice, h_rkpm_flat.begin(), h_rkpm_flat.end(),
                   d_rkpm_flat.begin());
    Gpu::streamSynchronize();
}

int mParticle::StartOfLagrangianMarker(size_t index) {
    return StartOfMarker.at(index);
}

int mParticle::NumOfLagrangianMarker(size_t index) {
    return NumOfMarker.at(index);
}

RealVect mParticle::GetPositionOfMarker(size_t index) {
    return LargrangianMarker.at(index);
}

void mParticle::WriteParticleFile(int index)
{
    mContainer->WriteAsciiFile(Concatenate("particle", index));
}

// Writes one particle's CSV line from the calling rank. The caller decides
// which rank writes which particle (see UpdateParticles).
void mParticle::WriteIBForceAndMoment(int step, Real time, Real dt, kernel& current_kernel)
{
    if (is_prescribed_fin(current_kernel.geometry_type)) {
        const std::string file = "IB_Fin_" + std::to_string(current_kernel.id)
                               + ".csv";
        const bool new_file = !fs::exists(file);
        std::ofstream output(file, std::ios::app);
        if (!output.is_open()) {
            Print() << "[DIBM] cannot open fin force file " << file << "\n";
            return;
        }
        if (new_file) {
            output << "step,time,phase,Fx,Fy,Fz,Mx,My,Mz\n";
        }

        const Real phase_direction = current_kernel.geometry_type == 4
                                   ? -1.0 : 1.0;
        const Real phase = phase_direction * 2.0 * Math::pi<Real>()
                         * current_kernel.fin_frequency * time
                         + current_kernel.fin_phase;
        const Real rho = ParticleProperties::euler_fluid_rho;
        output << step << ',' << time << ',' << phase << ','
               << -rho * current_kernel.ib_force[0] << ','
               << -rho * current_kernel.ib_force[1] << ','
               << -rho * current_kernel.ib_force[2] << ','
               << -rho * current_kernel.ib_moment[0] << ','
               << -rho * current_kernel.ib_moment[1] << ','
               << -rho * current_kernel.ib_moment[2] << '\n';
        return;
    }

    std::string file("IB_Particle_" + std::to_string(current_kernel.id) + ".csv");
    std::ofstream out_ib_force;

    std::string head;
    if(!fs::exists(file)){
        head = "iStep,time,X,Y,Z,Vx,Vy,Vz,Rx,Ry,Rz,Fx,Fy,Fz,Mx,My,Mz,Fcpx,Fcpy,Fcpz,Tcpx,Tcpy,Tcpz,SumUx,SumUy,SumUz,SumTx,SumTy,SumTz\n";
    }else{
        head = "";
    }

    out_ib_force.open(file, std::ios::app);
    if(!out_ib_force.is_open()){
        Print() << "[Particle] write particle file error , step: " << step;
    }else{
        out_ib_force << head << step << "," << time << ","
                     << current_kernel.location[0] << "," << current_kernel.location[1] << "," << current_kernel.location[2] << ","
                     << current_kernel.velocity[0] << "," << current_kernel.velocity[1] << "," << current_kernel.velocity[2] << ","
                     << current_kernel.omega[0] << "," << current_kernel.omega[1] << "," << current_kernel.omega[2] << ","
                     << current_kernel.ib_force[0] << "," << current_kernel.ib_force[1] << "," << current_kernel.ib_force[2] << ","
                     << current_kernel.ib_moment[0] << "," << current_kernel.ib_moment[1] << "," << current_kernel.ib_moment[2] << ","
                     << current_kernel.Fcp[0] << "," << current_kernel.Fcp[1] << "," << current_kernel.Fcp[2] << ","
                     << current_kernel.Tcp[0] << "," << current_kernel.Tcp[1] << "," << current_kernel.Tcp[2] << ","
                     << (current_kernel.sum_u_new[0] - current_kernel.sum_u_old[0])/dt << ","
                     << (current_kernel.sum_u_new[1] - current_kernel.sum_u_old[1])/dt << ","
                     << (current_kernel.sum_u_new[2] - current_kernel.sum_u_old[2])/dt << ","
                     << (current_kernel.sum_t_new[0] - current_kernel.sum_t_old[0])/dt << ","
                     << (current_kernel.sum_t_new[1] - current_kernel.sum_t_old[1])/dt << ","
                     << (current_kernel.sum_t_new[2] - current_kernel.sum_t_old[2])/dt << "\n";
    }
    out_ib_force.close();
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                    Particles member function                  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void Particles::create_particles(const Geometry &gm,
                                 const DistributionMapping & dm,
                                 const BoxArray & ba)
{
    BL_PROFILE("Particles::create_particles");
    Print() << "[Particle] : create Particle Container\n";
    if(particle->mContainer != nullptr){
        delete particle->mContainer;
        particle->mContainer = nullptr;
    }
    particle->mContainer = new mParticleContainer(gm, dm, ba);

    //get particle tile
    std::pair<int, int> key{0,0};
    auto& particleTileTmp = particle->mContainer->GetParticles(0)[key];
    //insert particle's markers
    BL_PROFILE_VAR("create_particles::push_back_markers", blp_push);
    // The container's tiles live in device memory: every push_back there is a
    // one-element device fill + stream synchronisation (~10 us), i.e. tens of
    // seconds for 1e5 markers x 11 attributes. Build the markers in a pinned
    // host tile instead and move them to the device with one copyParticles().
    using HostTile = amrex::ParticleTile<mParticleContainer::ParticleType, num_Real, num_Int, amrex::PinnedArenaAllocator>;
    HostTile host_tile;
    host_tile.define(0, 0);
    // Every rank walks the full particle list (so ids and start_id are the same
    // everywhere) but only inserts the markers of the particles whose centre
    // lies in one of its own grids. This keeps the initial Redistribute from
    // having to ship every marker out of rank 0 (~1 s per 1e5 markers).
    const auto dxinv = gm.InvCellSizeArray();
    const auto plo   = gm.ProbLoArray();
    const int myproc = ParallelDescriptor::MyProc();
    auto owns_centre = [&](const RealVect& loc) -> bool {
        IntVect iv;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            iv[d] = static_cast<int>(std::floor((loc[d] - plo[d]) * dxinv[d]));
        }
        const auto isects = ba.intersections(Box(iv, iv), /*first_only=*/true, /*ng=*/0);
        if (isects.empty()) {
            // centre outside the domain: fall back to the I/O rank
            return myproc == ParallelDescriptor::IOProcessorNumber();
        }
        return dm[isects[0].first] == myproc;
    };
    int marker_index = 0;
    {
        for (auto& cur_p: particle->particle_kernels){
            if (particle->do_RKPM) {
                cur_p.start_id = particle->StartOfLagrangianMarker(cur_p.id);
            }else {
                cur_p.start_id = marker_index;
            }
            if (!owns_centre(cur_p.location)) {
                marker_index += cur_p.ml;   // keep the global id numbering identical on all ranks
                continue;
            }
            for(int i = 0; i < cur_p.ml; i++){
                //insert code
                mParticleContainer::ParticleType markerP;
                markerP.id() = ++marker_index;
                markerP.cpu() = ParallelDescriptor::MyProc();
                markerP.pos(0) = cur_p.location[0];
                markerP.pos(1) = cur_p.location[1];
                markerP.pos(2) = cur_p.location[2];
                if (particle->do_RKPM) {
                    auto pos = particle->GetPositionOfMarker(i + cur_p.start_id);
                    markerP.pos(0) = pos[0];
                    markerP.pos(1) = pos[1];
                    markerP.pos(2) = pos[2];
                }

                std::array<ParticleReal, num_Real> Marker_attr{};
                Marker_attr[U_Marker] = 0.0;
                Marker_attr[V_Marker] = 0.0;
                Marker_attr[W_Marker] = 0.0;
                Marker_attr[Fx_Marker] = 0.0;
                Marker_attr[Fy_Marker] = 0.0;
                Marker_attr[Fz_Marker] = 0.0;

                std::array<int, num_Int> Particle_id;
                Particle_id[M_ID] = cur_p.id;

                host_tile.push_back(markerP);
                host_tile.push_back_real(Marker_attr);
                host_tile.push_back_int(Particle_id);
            }
        }
        const auto np_host = host_tile.numParticles();
        if (np_host > 0) {
            particleTileTmp.resize(np_host);
            amrex::copyParticles(particleTileTmp, host_tile);
        }
    }
    // start_id was computed identically on every rank above; no broadcast needed.
    BL_PROFILE_VAR_STOP(blp_push);
    BL_PROFILE_VAR("create_particles::Redistribute", blp_redist);
    particle->mContainer->Redistribute(); // Still needs to redistribute here!
    BL_PROFILE_VAR_STOP(blp_redist);

    ParticleProperties::plo = gm.ProbLoArray();
    ParticleProperties::phi = gm.ProbHiArray();
    ParticleProperties::dx = gm.CellSizeArray();
}

mParticle* Particles::get_particles()
{
    return particle;
}

void Particles::init_particle(Real gravity, Real h)
{
    BL_PROFILE("Particles::init_particle");
    Print() << "[Particle] : create Particle's kernel\n";
    particle = new mParticle;
    if(particle != nullptr){
        isInitial = true;
        particle->InitParticles(
            ParticleProperties::_x,
            ParticleProperties::_y,
            ParticleProperties::_z,
            ParticleProperties::_rho,
            ParticleProperties::Vx,
            ParticleProperties::Vy,
            ParticleProperties::Vz,
            ParticleProperties::Ox,
            ParticleProperties::Oy,
            ParticleProperties::Oz,
            ParticleProperties::TLX,
            ParticleProperties::TLY,
            ParticleProperties::TLZ,
            ParticleProperties::RLX,
            ParticleProperties::RLY,
            ParticleProperties::RLZ,
            ParticleProperties::_radius,
            ParticleProperties::_radius2,
            ParticleProperties::_radius3,
            ParticleProperties::_geometry_type,
            h,
            gravity,
            ParticleProperties::verbose);
    }

}

void Particles::Restart(Real gravity, Real h, int iStep)
{
    Print() << "[Particle] : restart Particle's kernel, step :" << iStep << "\n"
                   << "\tstart read particle csv file , default name is IB_Particle_x.csv\n"
                   << "\tdo not delete those file before \"restart\"\n\n";
    delete particle;
    particle = new mParticle;
            particle->InitParticles(
            ParticleProperties::_x,
            ParticleProperties::_y,
            ParticleProperties::_z,
            ParticleProperties::_rho,
            ParticleProperties::Vx,
            ParticleProperties::Vy,
            ParticleProperties::Vz,
            ParticleProperties::Ox,
            ParticleProperties::Oy,
            ParticleProperties::Oz,
            ParticleProperties::TLX,
            ParticleProperties::TLY,
            ParticleProperties::TLZ,
            ParticleProperties::RLX,
            ParticleProperties::RLY,
            ParticleProperties::RLZ,
            ParticleProperties::_radius,
            ParticleProperties::_radius2,
            ParticleProperties::_radius3,
            ParticleProperties::_geometry_type,
            h,
            gravity,
            ParticleProperties::verbose);
    //deal in IO processor
    //start read csv file
    for(auto& kernel : particle->particle_kernels){
        if (is_prescribed_fin(kernel.geometry_type)) {
            // Prescribed geometry is reconstructed analytically and does not
            // carry rigid-body state in the particle CSV file.
            continue;
        }
        //filename
        if(ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()){
            std::string fileName = "IB_Particle_" + std::to_string(kernel.id) + ".csv";
            std::string tmpfile = "tmp" + fileName;
            //file stream
            std::ifstream particle_data(fileName);
            std::ofstream particle_file(tmpfile);
            // open state
            if(!particle_data.is_open() || !particle_file.is_open()){
                Abort("\tCan not open particle file : " + fileName);
            }
            std::string lineData;
            int line{0};
            while(std::getline(particle_data, lineData)){
                line++;
                particle_file << lineData << "\n";
                if(line <= iStep) {
                    continue;
                }
                //old location
                //iStep,time,X,Y,Z,Vx,Vy,Vz,Rx,Ry,Rz,Fx,Fy,Fz,Mx,My,Mz,Fcpx,Fcpy,Fcpz,Tcpx,Tcpy,Tcpz
                if(line == iStep + 1){
                    std::stringstream ss(lineData);
                    std::string data;
                    std::vector<Real> dataStruct;
                    while(std::getline(ss, data, ',')){
                        dataStruct.emplace_back(std::stod(data));
                    }
                    kernel.location[0] = dataStruct[2];
                    kernel.location[1] = dataStruct[3];
                    kernel.location[2] = dataStruct[4];
                    kernel.velocity[0] = dataStruct[5];
                    kernel.velocity[1] = dataStruct[6];
                    kernel.velocity[2] = dataStruct[7];
                    kernel.omega[0] = dataStruct[8];
                    kernel.omega[1] = dataStruct[9];
                    kernel.omega[2] = dataStruct[10];

                    kernel.location_old = kernel.location;
                    kernel.velocity_old = kernel.velocity;
                    kernel.omega_old    = kernel.omega;
                    break;
                }
                else
                    break;
            }
            particle_data.close();
            particle_file.close();
            std::remove(fileName.c_str());
            std::rename(tmpfile.c_str(), fileName.c_str());
        }
        ParallelDescriptor::Bcast(&kernel.location[0], 3, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.velocity[0], 3,ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.omega[0], 3,ParallelDescriptor::IOProcessorNumber());

        ParallelDescriptor::Bcast(&kernel.location_old[0], 3, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.velocity_old[0], 3,ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.omega_old[0], 3,ParallelDescriptor::IOProcessorNumber());
    }

    isInitial = true;
}

void Particles::Initialize()
{
    ParmParse pp("particle");

    std::string particle_inputfile;
    std::string particle_init_file;
    pp.get("input",particle_inputfile);

    if(!particle_inputfile.empty()){
        ParmParse p_file(particle_inputfile);
        p_file.query("init", particle_init_file);
        p_file.getarr("x",          ParticleProperties::_x);
        p_file.getarr("y",          ParticleProperties::_y);
        p_file.getarr("z",          ParticleProperties::_z);
        p_file.getarr("rho",        ParticleProperties::_rho);
        p_file.getarr("velocity_x", ParticleProperties::Vx);
        p_file.getarr("velocity_y", ParticleProperties::Vy);
        p_file.getarr("velocity_z", ParticleProperties::Vz);
        p_file.getarr("omega_x",    ParticleProperties::Ox);
        p_file.getarr("omega_y",    ParticleProperties::Oy);
        p_file.getarr("omega_z",    ParticleProperties::Oz);
        p_file.getarr("TLX",        ParticleProperties::TLX);
        p_file.getarr("TLY",        ParticleProperties::TLY);
        p_file.getarr("TLZ",        ParticleProperties::TLZ);
        p_file.getarr("RLX",        ParticleProperties::RLX);
        p_file.getarr("RLY",        ParticleProperties::RLY);
        p_file.getarr("RLZ",        ParticleProperties::RLZ);
        p_file.getarr("radius",     ParticleProperties::_radius);
        p_file.queryarr("radius2",   ParticleProperties::_radius2);
        p_file.queryarr("radius3",   ParticleProperties::_radius3);
        p_file.queryarr("geometry_type", ParticleProperties::_geometry_type);
        p_file.query("fin_length", ParticleProperties::fin_length);
        p_file.query("fin_span", ParticleProperties::fin_span);
        p_file.query("fin_amplitude_deg", ParticleProperties::fin_amplitude_deg);
        p_file.query("fin_frequency", ParticleProperties::fin_frequency);
        p_file.query("fin_wavelength", ParticleProperties::fin_wavelength);
        p_file.query("fin_phase", ParticleProperties::fin_phase);
        p_file.query("fin_wave_amplitude", ParticleProperties::fin_wave_amplitude);
        p_file.query("fin_wave_number", ParticleProperties::fin_wave_number);
        p_file.queryarr("fin_wave_phase_deg",
                        ParticleProperties::fin_wave_phase_deg);
        p_file.query("fin_n_chord", ParticleProperties::fin_n_chord);
        p_file.query("fin_n_span", ParticleProperties::fin_n_span);
        p_file.query("RD",          ParticleProperties::rd);
        p_file.query("LOOP_NS",     ParticleProperties::loop_ns);
        p_file.query("LOOP_SOLID",  ParticleProperties::loop_solid);
        p_file.query("verbose",     ParticleProperties::verbose);
        p_file.query("start_step",  ParticleProperties::start_step);
        p_file.query("Uhlmann",     ParticleProperties::Uhlmann);
        p_file.query("collision_model", ParticleProperties::collision_model);
        p_file.query("write_freq",  ParticleProperties::write_freq);
        p_file.query("delta_type", ParticleProperties::delta_type);
        p_file.query("reduce_method", ParticleProperties::reduce_method);
        p_file.query("pvf_method", ParticleProperties::pvf_method);
        // update with RKPM method
        p_file.query("RKPM", ParticleProperties::RKPM);

        ParmParse ns("ns");
        ns.get("fluid_rho",      ParticleProperties::euler_fluid_rho);

        ParmParse level_parse("amr");
        level_parse.get("max_level", ParticleProperties::euler_finest_level);

        ParmParse geometry_parse("geometry");
        geometry_parse.getarr("prob_lo", ParticleProperties::GLO);
        geometry_parse.getarr("prob_hi", ParticleProperties::GHI);
        Print() << "[Particle] : Reading partilces cfg file : " << particle_inputfile << "\n"
                       << "             Particle's level : " << ParticleProperties::euler_finest_level << "\n";

        if(!particle_init_file.empty()){
            ParticleProperties::init_particle_from_file = true;
            //clear particle position container
            ParticleProperties::_x.clear();
            ParticleProperties::_y.clear();
            ParticleProperties::_z.clear();
            // parse particle's location data
            std::ifstream init_particle(particle_init_file);
            std::string line_data;
            while(std::getline(init_particle, line_data)){
                // id x_location y_location z_location
                std::istringstream line(line_data);
                std::vector<std::string> str_tokne;
                std::string token;
                while(line >> token){
                    str_tokne.push_back(token);
                }

                ParticleProperties::_x.push_back(std::stod(str_tokne[0]));
                ParticleProperties::_y.push_back(std::stod(str_tokne[1]));
                ParticleProperties::_z.push_back(std::stod(str_tokne[2]));
            }
            ParticleProperties::_x.shrink_to_fit();
            ParticleProperties::_y.shrink_to_fit();
            ParticleProperties::_z.shrink_to_fit();
            Print() << "             initial Particle by file : " << particle_init_file
                    << "             particle's size : " << ParticleProperties::_x.size() << "\n";
        }

    }else {
        Abort("[Particle] : can't read particles settings, pls check your config file \"particle.input\"");
    }
}

int Particles::ParticleFinestLevel()
{
    return ParticleProperties::euler_finest_level;
}
