/*
 * SPDX-FileCopyrightText: 1997 - 2023 Berkeley Lab; 2023 - 2026 Yadong Zeng<zdsjtu@gmail.com> & ZhuXu Li<1246206018@qq.com>
 *
 * SPDX-License-Identifier: LicenseRef-OpenSource
 * Modified from IAMR, originally developed at Lawrence Berkeley National Lab.
 * Original source: https://github.com/AMReX-Fluids/IAMR
 */

#ifndef IAMR_NavierStokes_H_
#define IAMR_NavierStokes_H_

#include <NavierStokesBase.H>
#include <AMReX_LevelBld.H>
#include "Utilities.H"

class NSBld
    :
    public amrex::LevelBld
{
    void variableSetUp () override;
    void variableCleanUp () override;
    amrex::AmrLevel *operator() () override;
    amrex::AmrLevel *operator() (amrex::Amr&            papa,
                                 int             lev,
                                 const amrex::Geometry& level_geom,
                                 const amrex::BoxArray& ba,
                                 const amrex::DistributionMapping& dm,
                                 amrex::Real            time) override;
};

/**
 * @defgroup NavierStokes fluid parameters and function
 *
 * introduction for Navier stokes
 */

/// @{

class NavierStokes
    :
    public NavierStokesBase
{
public:

    friend class Projection;
    friend class MacProj;
    friend class Diffusion;

    NavierStokes ();

    NavierStokes (amrex::Amr&            papa,
                  int             lev,
                  const amrex::Geometry& level_geom,
                  const amrex::BoxArray& bl,
                  const amrex::DistributionMapping& dm,
                  amrex::Real            time);

    ~NavierStokes () override = default;

    NavierStokes (NavierStokes const&) = delete;
    NavierStokes (NavierStokes &&) = delete;
    NavierStokes& operator= (NavierStokes const&) = delete;
    NavierStokes& operator= (NavierStokes &&) = delete;

    ////////////////////////////////////////////////////////////////////////////
    //    AmrLevel virtual functions                                          //
    ////////////////////////////////////////////////////////////////////////////
    //
    // Init grid data at problem start-up.
    //
    /**
     * Initialize grid data at problem start-up.
     * This function sets up initial conditions for the Navier-Stokes equations
     * including velocity, pressure, and scalar fields.
     */
    void initData () override;
    //
    // Problem initialization functions
    //
#include "prob_init.H"
    //
    // Build/fill any additional data structures after restart.
    //
    void post_restart () override;
    //
    // Do pre-plotfile work
    //
    void writePlotFilePre (const std::string& dir,
                           std::ostream&      os) override;
    //
    // Do post-plotfile work
    //
    void writePlotFilePost (const std::string& dir,
                            std::ostream&      os) override;
    //
    // Modify list of variables to be plotted
    //
    void setPlotVariables() override;
    //
    // Returns a amrex::MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the amrex::MultiFab is built on the appropriately
    // grown BoxArray.
    //
    std::unique_ptr<amrex::MultiFab>
    derive (const std::string& name, amrex::Real time, int ngrow) override;
    //
    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    //
    void derive (const std::string& name,
                 amrex::Real        time,
                 amrex::MultiFab&   mf,
                 int                dcomp) override;
    //
    // Insure state, and pressure are consistent.
    //
    void post_init (amrex::Real stop_time) override;
    //
    // Advance grids at this level in time.
    //
    /**
     * Advance grids at this level in time using the Navier-Stokes solver.
     * @param time Current simulation time
     * @param dt Time step size
     * @param iteration Current iteration number
     * @param ncycle Total number of cycles for this level
     * @return Estimated optimal timestep for next iteration
     */
    amrex::Real advance (amrex::Real time,
                         amrex::Real dt,
                         int  iteration,
                         int  ncycle) override;

    // ls related
    /**
     * Advance grids using semi-staggered two-phase level set method.
     * This function implements the level set method for two-phase flow simulations.
     * @param time Current simulation time
     * @param dt Time step size
     * @param iteration Current iteration number
     * @param ncycle Total number of cycles for this level
     * @return Estimated optimal timestep for next iteration
     */
    amrex::Real advance_semistaggered_twophase_ls (amrex::Real time,
                         amrex::Real dt,
                         int  iteration,
                         int  ncycle);
#ifdef AMREX_PARTICLES
    /**
     * Advance grids using semi-staggered fluid-structure interaction with diffused immersed boundary method.
     * This function implements FSI simulations using the diffused immersed boundary method.
     * @param time Current simulation time
     * @param dt Time step size
     * @param iteration Current iteration number
     * @param ncycle Total number of cycles for this level
     * @return Estimated optimal timestep for next iteration
     */
    amrex::Real advance_semistaggered_fsi_diffusedib (amrex::Real time,
                         amrex::Real dt,
                         int  iteration,
                         int  ncycle);
#endif

    /**
     * Advance grids using semi-staggered two-phase phase field method.
     * This function implements the phase field method for two-phase flow simulations.
     * @param time Current simulation time
     * @param dt Time step size
     * @param iteration Current iteration number
     * @param ncycle Total number of cycles for this level
     * @return Estimated optimal timestep for next iteration
     */
    amrex::Real advance_semistaggered_twophase_phasefield (amrex::Real time,
                         amrex::Real dt,
                         int  iteration,
                         int  ncycle);

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase public virtual functions                           //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Compute the mac sync correction.
    //
    void mac_sync () override;
    //
    // Reflux function.
    //
    void reflux () override;
    //
    // Calculate divU, which sets them to zero by default.
    //
    void calc_divu (amrex::Real      time,
                    amrex::Real      dt,
                    amrex::MultiFab& divu) override;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokes public static functions                                //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Define data descriptors.
    //
    static void variableSetUp ();

    //
    // ls related
    //
    // Max possible number of scalars is 5: Density, 2 Tracers, Temperature, ls
    static constexpr int NUM_STATE_MAX = AMREX_SPACEDIM+5;
    static amrex::GpuArray<amrex::GpuArray<amrex::Real, NUM_STATE_MAX>, AMREX_SPACEDIM*2>
    //static auto
    get_bc_values () { return m_bc_values; }

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokes public functions                                       //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Advect scalars.
    //
    /**
     * Advect scalar quantities using the chosen advection scheme.
     * @param dt Time step size
     * @param fscalar First scalar component index
     * @param lscalar Last scalar component index
     */
    void scalar_advection (amrex::Real dt,
                           int  fscalar,
                           int  lscalar);
    //
    /**
     * Update scalar quantities with diffusion terms.
     * @param dt Time step size
     * @param first_scalar First scalar component index
     * @param last_scalar Last scalar component index
     */
    void scalar_diffusion_update (amrex::Real dt,
                                  int  first_scalar,
                                  int  last_scalar);

protected:

    ////////////////////////////////////////////////////////////////////////////
    //    Overriding Virtual Functions in NavierStokesBase                    //
    ////////////////////////////////////////////////////////////////////////////

    /**
     * Average down all state types from fine to coarse level.
     */
    void avgDown () override;
    //
    // Calculate nonuniform viscosity and diffusivity
    //
    /**
     * Calculate nonuniform viscosity and diffusivity.
     * @param time Current simulation time
     * @param dt Time step size
     * @param iteration Current iteration number
     * @param ncycle Total number of cycles
     */
    void calcViscosity (amrex::Real time,
                        amrex::Real dt,
                        int  iteration,
                        int  ncycle) override;
    //
    /**
     * Calculate diffusivity coefficients.
     * @param time Current simulation time
     */
    void calcDiffusivity (amrex::Real time) override;
    //
    /**
     * Get viscosity coefficients for each spatial direction.
     * @param viscosity Array of MultiFab pointers for viscosity in each direction
     * @param time Current simulation time
     */
    void getViscosity (amrex::MultiFab*  viscosity[AMREX_SPACEDIM],
                       amrex::Real time) override;
    //
    // Compute viscous terms.
    //
    /**
     * Compute viscous terms for momentum equations.
     * @param visc_terms MultiFab containing viscous terms (output)
     * @param src_comp Source component index
     * @param ncomp Number of components
     * @param time Current simulation time
     */
    void getViscTerms (amrex::MultiFab& visc_terms,
                       int       src_comp,
                       int       ncomp,
                       amrex::Real      time) override;
    //
    /**
     * Sum integrated quantities for diagnostics and statistics.
     */
    void sum_integrated_quantities () override;
    //
    /**
     * Update velocities with diffusion terms.
     * @param dt Time step size
     */
    void velocity_diffusion_update (amrex::Real dt) override;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokes protected static functions                             //
    ////////////////////////////////////////////////////////////////////////////

    /**
     * Estimate error for adaptive mesh refinement.
     * @param tags TagBoxArray for marking cells for refinement
     * @param clearval Value to clear tags
     * @param tagval Value to set tags
     * @param time Current simulation time
     * @param n_error_buf Number of error buffer cells
     * @param ngrow Number of ghost cells
     */
    void errorEst (amrex::TagBoxArray& tags,
                   int                 clearval,
                   int                 tagval,
                   amrex::Real         time,
                   int                 n_error_buf,
                   int                 ngrow) override;

    /**
     * Initialize NavierStokes solver by reading input file.
     */
    static void Initialize ();
    /**
     * Initialize boundary conditions.
     */
    static void Initialize_bcs ();
    /**
     * Initialize diffusivity coefficients.
     */
    static void Initialize_diffusivities ();
    /**
     * Finalize NavierStokes solver and clean up resources.
     */
    static void Finalize ();

private:

    /**
     * Get diffusivity coefficients for scalar transport.
     * @param diffusivity Array of MultiFab pointers for diffusivity in each direction
     * @param time Current simulation time
     * @param state_comp State component index
     * @param dst_comp Destination component index
     * @param ncomp Number of components
     */
    static void getDiffusivity (amrex::MultiFab*  diffusivity[AMREX_SPACEDIM],
                                amrex::Real time,
                                int  state_comp,
                                int  dst_comp,
                                int  ncomp);

    NavierStokes& getLevel (int lev) {
    return dynamic_cast<NavierStokes&> ( parent->getLevel(lev) );
    }

    //
    // Initialize the pressure by iterating the initial timestep.
    //
    /**
     * Initialize pressure by iterating the initial timestep.
     * @param dt_init Initial time step size (output)
     * @param nc_save Vector to save cycle numbers
     * @param dt_save Vector to save time step sizes
     */
    void post_init_press (amrex::Real&        dt_init,
                          amrex::Vector<int>&  nc_save,
                          amrex::Vector<amrex::Real>& dt_save);
    //
    // Update scalars, (viscous solve in scalar_update).
    //
    /**
     * Update scalar quantities with viscous solve.
     * @param dt Time step size
     * @param first_scalar First scalar component index
     * @param last_scalar Last scalar component index
     */
    void scalar_update (amrex::Real dt,
                        int  first_scalar,
                        int  last_scalar);
    //
    // Define error estimators
    //
    /**
     * Set up error estimators for adaptive mesh refinement.
     */
    static void error_setup ();

    ///////////////////
    // Private data  //
    ///////////////////

    //
    // Runtime parameters
    //
    // For EB, set all covered cells to zero (1, default), or show whatever is
    // in covered cells in plotfile (0).
    static int set_plot_coveredCell_val;

    //
    // Holds info for dynamically generated error estimators
    //
    static amrex::Vector<amrex::AMRErrorTag> errtags;

    //
    // Hold external Dirichlet BC values
    //
    //
    // ls related
    //
    // Max possible number of scalars is 5: Density, 2 Tracers, Temperature, ls
    static amrex::GpuArray<amrex::GpuArray<amrex::Real, NUM_STATE_MAX>,
                   AMREX_SPACEDIM*2> m_bc_values;
};

/// @}
#endif /*_NavierStokes_H_*/
