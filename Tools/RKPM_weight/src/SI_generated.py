import numpy as np
from src import visual
from src.grid_index import containing_cell_indices

# Read Lagrangian marker coordinates (xp, yp, zp) from a whitespace-delimited
# point-cloud file. Lines beginning with # are comments.
#   # X Y Z
#   -0.0087703851 0.0371783379 0.0369941932
#   ...
def read_geometry_file(file_path):
    points = np.atleast_2d(np.loadtxt(file_path, comments='#'))
    if points.size == 0:
        raise ValueError(f"Geometry file is empty or cannot be parsed: {file_path}")
    if points.shape[1] != 3:
        raise ValueError(
            f"Geometry file must have 3 columns (X Y Z), got "
            f"{points.shape[1]}: {file_path}")
    return points

# Approximate ellipsoid surface area with the Knud Thomsen formula.
def ellipsoid_area_approx(a, b, c):
    """
    Approximate ellipsoid surface area with the Knud Thomsen formula.
    """
    p = 1.6075
    term = ((a * b)**p + (a * c)**p + (b * c)**p) / 3.0
    return 4 * np.pi * (term**(1 / p))

def cylinder_area(radius, height):
    """
    Compute total cylinder surface area, including both end caps.

    Args:
        radius (float): Cylinder radius.
        height (float): Cylinder height.

    Returns:
        float: Total area = 2 * pi * radius * (radius + height).
    """
    return 2 * np.pi * radius * (radius + height)

# Generate affected Eulerian cells and Lagrangian marker support domains.
#
# Design: operate directly in the global coordinate system of the solver's
# finest grid without selecting a local subregion or applying a local-to-global
# index offset. The old int(sx/dx) truncation was ambiguous.
#   1. Read prob_lo/prob_hi/n_cell/max_level and compute
#      dx_finest = (prob_hi - prob_lo) / (n_cell * 2**max_level).
#   2. Expand the point-cloud bounds by two finest cells and clamp the result to
#      the computational domain.
#   3. Generate global cell centers as (index + 0.5) * dx + prob_lo. RKPM weights
#      depend only on relative coordinates and are independent of the origin.
#
# Inputs: domain bounds, finest-cell size, point-cloud path, optional body-frame
# center, z-axis rotation angle in degrees, and marker-surface ellipsoid semi-axes.
# Outputs include affected Eulerian points, transformed Lagrangian points,
# nearest grid points, support parameters and all per-marker support domains.
def generate_grid(
    prob_lo, prob_hi, dx_finest, geometry_file=None, center=None, angle=0.0, *, axes=None
):

    prob_lo = np.asarray(prob_lo, dtype=float)
    prob_hi = np.asarray(prob_hi, dtype=float)
    dx_finest = np.asarray(dx_finest, dtype=float)
    dx, dy, dz = dx_finest
    # Number of finest-grid cells per direction, used for domain clamping.
    n_fine = np.round((prob_hi - prob_lo) / dx_finest).astype(int)

    # Lagrangian point cloud.
    if geometry_file is not None:
        lagrangian_points = read_geometry_file(geometry_file)
    else:
        raise NameError("A Lagrangian marker-coordinate file is required")
        # lagrangian_points = load_lagrangian_from_id_file('rkpm_mappings.id')

    # Estimate missing semi-axes from the input cloud before applying --angle.
    if axes is None:
        axes = np.ptp(lagrangian_points, axis=0) / 2.0
        axis_source = "input point-cloud bounds"
    else:
        axes = np.asarray(axes, dtype=float)
        if axes.shape != (3,) or not np.all(np.isfinite(axes) & (axes > 0.0)):
            raise ValueError("axes must contain three finite, positive ellipsoid semi-axes")
        axis_source = "specified"
    print(f"[geometry] ellipsoid semi-axes ({axis_source}): {axes.tolist()}")

    # Translate body-frame coordinates to world coordinates.
    if center is not None:
        center = np.asarray(center, dtype=float)
        lagrangian_points = lagrangian_points + center
    # Rotate around z after translation; use center as the body-frame pivot.
    if angle:
        pivot = center if center is not None else np.zeros(3)
        theta = np.radians(angle)
        c, s = np.cos(theta), np.sin(theta)
        Rz = np.array([[ c, -s, 0.],
                       [ s,  c, 0.],
                       [0., 0., 1.]])
        lagrangian_points = (lagrangian_points - pivot) @ Rz.T + pivot
    # visual.PointCloud(lagrangian_points)

    # Bounding region: expand point-cloud extrema by two finest-grid cells.
    pmin = lagrangian_points.min(axis=0)
    pmax = lagrangian_points.max(axis=0)
    # The point cloud must remain inside the computational domain.
    if np.any(pmin < prob_lo) or np.any(pmax > prob_hi):
        raise ValueError(
            f"Point cloud extends outside [prob_lo, prob_hi]:\n"
            f"  prob_lo = {prob_lo.tolist()}\n  prob_hi = {prob_hi.tolist()}\n"
            f"  point-cloud min = {pmin.tolist()}\n"
            f"  point-cloud max = {pmax.tolist()}\n"
            "Check the --geometry coordinate system or use --body-frame to "
            "translate the point cloud into the domain."
        )
    # Expand finest-grid global indices and clamp them to [0, n_fine-1].
    min_ijk = containing_cell_indices(pmin, prob_lo, dx_finest)
    max_ijk = containing_cell_indices(pmax, prob_lo, dx_finest)
    i_lo, j_lo, k_lo = (int(value) - 2 for value in min_ijk)
    i_hi, j_hi, k_hi = (int(value) + 2 for value in max_ijk)
    i_lo, j_lo, k_lo = max(0, i_lo), max(0, j_lo), max(0, k_lo)
    i_hi = min(int(n_fine[0]) - 1, i_hi)
    j_hi = min(int(n_fine[1]) - 1, j_hi)
    k_hi = min(int(n_fine[2]) - 1, k_hi)
    print(f"[grid] finest index range: i[{i_lo},{i_hi}] j[{j_lo},{j_hi}] k[{k_lo},{k_hi}] "
          f"-> cells ({i_hi-i_lo+1}, {j_hi-j_lo+1}, {k_hi-k_lo+1})")

    cell_volume = float(np.prod(dx_finest))

    # Estimate marker area and thickness.
    Ne = len(lagrangian_points)
    area = ellipsoid_area_approx(*axes) / Ne
    thickness = min(dx, dy, dz)
    V_lag = area * thickness
    print(f"Vl: {V_lag}, area:{area}, frac: {V_lag / cell_volume}")

    # Find each marker's containing global cell.
    indices_ijk = containing_cell_indices(
        lagrangian_points, prob_lo, dx_finest
    )
    nearest_grid_points = prob_lo + (indices_ijk + 0.5) * dx_finest

    # Compute support-domain scales.
    delta_I = np.full(Ne, dx + (1 / 1000) * dx)
    eta_I = np.full(Ne, dy + (1 / 1000) * dy)
    theta_I = np.full(Ne, dz + (1 / 1000) * dz)

    # Construct every 3x3x3 support by broadcasting integer offsets. This is
    # equivalent to searching the Eulerian point cloud around each marker, but
    # scales as O(27*Nmarkers) instead of O(Ngrid*Nmarkers).
    offset_grid = np.stack(
        np.meshgrid(
            np.arange(-1, 2),
            np.arange(-1, 2),
            np.arange(-1, 2),
            indexing="ij",
        ),
        axis=-1,
    ).reshape(-1, 3)
    support_indices = indices_ijk[:, None, :] + offset_grid[None, :, :]
    valid_support = np.all(
        (support_indices >= 0) & (support_indices < n_fine), axis=2
    )
    support_points = prob_lo + (support_indices + 0.5) * dx_finest

    all_S_I = [
        np.column_stack(
            (
                support_points[index, valid_support[index]],
                np.full(np.count_nonzero(valid_support[index]), cell_volume),
            )
        )
        for index in range(Ne)
    ]

    # Post-solve conservation checks need only cells touched by a stencil, not
    # every cell in the point-cloud bounding box.
    affected_indices = np.unique(support_indices[valid_support], axis=0)
    eulerian_points = prob_lo + (affected_indices + 0.5) * dx_finest

    return eulerian_points, Ne, lagrangian_points, nearest_grid_points, delta_I, eta_I, theta_I, all_S_I, V_lag
