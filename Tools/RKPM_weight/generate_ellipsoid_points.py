#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 IAMReX contributors
# SPDX-License-Identifier: BSD-3-Clause
"""Generate approximately area-uniform markers on an axis-aligned ellipsoid.

The output matches Tutorials/RKPM/test_cylinder_centers.txt: '# X Y Z',
then exactly N rows of Cartesian coordinates, with no point-count column.
Default placement minimizes E3=sum(1/r_ij**3) on the ellipsoid.
The command line and Python API interpret axes as x/y/z semi-axis lengths
in the supplied order.
Requires NumPy and SciPy. Coordinates use the same units as the supplied axes.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.special import elliprg
from scipy.stats import qmc


def _surface_samples(count, axes, seed):
    """Map sphere samples using the ellipsoid surface-area Jacobian.

    Under x=diag(a,b,c) u, dS=abc*||u/axes|| dOmega. Accepting a
    uniform sphere direction with probability min(axes)*||u/axes||
    therefore samples ellipsoid area, rather than solid angle.
    Scrambled Sobol samples reduce quadrature noise and are reproducible.
    """
    engine = qmc.Sobol(d=3, scramble=True, seed=seed)
    batch = 1 << max(10, (2 * count - 1).bit_length())
    accepted = []
    remaining = count
    while remaining:
        draws = engine.random(batch)
        z = 2.0 * draws[:, 0] - 1.0
        phi = 2.0 * np.pi * draws[:, 1]
        radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
        directions = np.column_stack((radius * np.cos(phi), radius * np.sin(phi), z))
        probability = np.linalg.norm(directions * (axes.min() / axes), axis=1)
        selected = directions[draws[:, 2] < probability][:remaining]
        accepted.append(selected * axes)
        remaining -= len(selected)
        # Doubling the total Sobol draw count preserves its balance property.
        batch = engine.num_generated
    return np.concatenate(accepted)


def _riesz3_energy_gradient(angles, axes, block_size=256):
    """Mean pair E3 energy and analytic angular gradient on a scaled ellipsoid.

    Blocking bounds temporary storage to O(block_size*N); all pairs remain
    included. Mean-pair normalization does not change the minimizers.
    """
    theta, phi = np.asarray(angles).reshape(-1, 2).T
    sine, cosine = np.sin(theta), np.cos(theta)
    sin_phi, cos_phi = np.sin(phi), np.cos(phi)
    points = np.column_stack((sine*cos_phi, sine*sin_phi, cosine)) * axes
    count = len(points)
    pair_count = count * (count - 1) / 2
    cartesian_gradient = np.empty_like(points)
    energy = 0.0
    for start in range(0, count, block_size):
        stop = min(start + block_size, count)
        squared = cdist(points[start:stop], points, "sqeuclidean")
        squared[np.arange(stop-start), np.arange(start, stop)] = np.inf
        if np.any(squared == 0):
            return np.inf, np.zeros(2*count)
        inverse = squared**-1.5
        energy += float(inverse.sum())
        weights = inverse / squared
        # optimize=False avoids threaded BLAS overhead for this three-column sum.
        weighted_points = np.einsum("ij,jk->ik", weights, points, optimize=False)
        cartesian_gradient[start:stop] = -3 * (
            weights.sum(axis=1)[:, None] * points[start:stop] - weighted_points
        ) / pair_count
    theta_derivative = np.column_stack((cosine*cos_phi, cosine*sin_phi, -sine)) * axes
    phi_derivative = np.column_stack((-sine*sin_phi, sine*cos_phi, np.zeros(count))) * axes
    gradient = np.column_stack((np.sum(cartesian_gradient*theta_derivative, axis=1),
                                np.sum(cartesian_gradient*phi_derivative, axis=1)))
    return 0.5 * energy / pair_count, gradient.ravel()


def _relax_riesz3(points, axes, iterations):
    """Minimize E3 with an angular surface constraint and L-BFGS-B."""
    directions = points / axes
    angles = np.column_stack((np.arccos(np.clip(directions[:, 2], -1, 1)),
                              np.arctan2(directions[:, 1], directions[:, 0])))
    initial_energy = _riesz3_energy_gradient(angles.ravel(), axes)[0]
    result = minimize(_riesz3_energy_gradient, angles.ravel(), args=(axes,), jac=True,
                      method="L-BFGS-B", options={"maxiter": iterations, "maxls": 40,
                      "maxfun": 50*iterations + 1, "ftol": 2e-13, "gtol": 1e-8,
                      "maxcor": 15})
    if (not np.isfinite(result.fun) or not np.all(np.isfinite(result.x))
            or not np.all(np.isfinite(result.jac)) or result.fun > initial_energy * (1+1e-12)):
        raise ValueError("E3 optimization produced a nonfinite or energy-increasing result")
    theta, phi = result.x.reshape(-1, 2).T
    points = np.column_stack((np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi),
                              np.cos(theta))) * axes
    diagnostics = {
        "converged": bool(result.success), "message": str(result.message),
        "iterations": int(result.nit), "evaluations": int(result.nfev),
        "initial_normalized_mean_pair_energy": float(initial_energy),
        "final_normalized_mean_pair_energy": float(result.fun),
        "max_angular_gradient": float(np.max(np.abs(result.jac))),
    }
    return points, diagnostics


def generate_ellipsoid_points(
    num_points, axes, center=(0.0, 0.0, 0.0), *, iterations=None,
    samples_per_point=64, seed=0, method="riesz3", return_diagnostics=False,
):
    """Return an (N,3) array using E3 repulsion by default.

    E3=sum(1/r_ij**3) has repulsive force magnitude proportional to 1/r**4.
    Angular parameterization keeps every iterate on the ellipsoid. The old
    sampled-centroid method remains available with method='centroid'. Neither
    method certifies equal-area cells or a global energy minimum.
    With return_diagnostics=True, return (points, optimization_diagnostics).
    """
    axes = np.asarray(axes, dtype=float)
    center = np.asarray(center, dtype=float)
    if (not isinstance(num_points, (int, np.integer)) or num_points <= 0):
        raise ValueError("num_points must be a positive integer")
    if axes.shape != (3,) or not np.all(np.isfinite(axes)) or np.any(axes <= 0):
        raise ValueError("axes must contain three finite positive semi-axis lengths")
    if center.shape != (3,) or not np.all(np.isfinite(center)):
        raise ValueError("center must contain three finite coordinates")
    if method not in ("riesz3", "centroid"):
        raise ValueError("method must be 'riesz3' or 'centroid'")
    if iterations is None:
        iterations = 3000 if method == "riesz3" else 80
    if not isinstance(iterations, (int, np.integer)) or iterations < 0:
        raise ValueError("iterations must be a nonnegative integer")
    if (not isinstance(samples_per_point, (int, np.integer)) or samples_per_point < 8):
        raise ValueError("samples_per_point must be an integer >= 8")
    if not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")

    # Work in units of the largest semi-axis, independent of coordinate units.
    scale = axes.max()
    unit_axes = axes / scale
    points = _surface_samples(num_points, unit_axes, seed)
    diagnostics = {
        "method": method, "seed": int(seed), "requested_max_iterations": int(iterations),
        "iterations": 0, "converged": num_points == 1,
        "message": "Single marker; no pair interaction" if num_points == 1 else "Relaxation disabled",
        "global_minimum_certified": False,
    }
    if method == "riesz3":
        diagnostics.update(energy_exponent=3, force_magnitude_power=4,
                           energy_coordinate_scale=float(scale))
    if iterations and num_points > 1 and method == "riesz3":
        points, optimization = _relax_riesz3(points, unit_axes, iterations)
        diagnostics.update(optimization)
    elif iterations and num_points > 1:
        samples = _surface_samples(max(4096, samples_per_point * num_points),
                                   unit_axes, seed + 1)
        for iteration in range(iterations):
            _, labels = cKDTree(points).query(samples)
            counts = np.bincount(labels, minlength=num_points)
            centroids = points.copy()
            occupied = counts > 0
            for axis in range(3):
                sums = np.bincount(labels, weights=samples[:, axis], minlength=num_points)
                centroids[occupied, axis] = sums[occupied] / counts[occupied]
            normals = points / unit_axes**2
            normals /= np.linalg.norm(normals, axis=1)[:, None]
            movement = centroids - points
            movement -= np.einsum("ij,ij->i", movement, normals)[:, None] * normals
            # Limit the step in ellipsoid coordinates before radial retraction.
            lengths = np.linalg.norm(movement / unit_axes, axis=1)
            movement *= np.minimum(0.8, 0.15 / np.maximum(lengths, 1e-30))[:, None]
            updated = points + movement
            updated /= np.linalg.norm(updated / unit_axes, axis=1)[:, None]
            change = np.max(np.linalg.norm(updated - points, axis=1))
            points = updated
            if change < 1e-6 / np.sqrt(num_points):
                break
        diagnostics.update(iterations=iteration+1, converged=bool(change < 1e-6 / np.sqrt(num_points)),
                           message="Centroid displacement tolerance reached" if change < 1e-6 / np.sqrt(num_points)
                           else "Centroid iteration limit reached",
                           max_displacement_normalized=float(change))
    points = points * scale + center
    return (points, diagnostics) if return_diagnostics else points


def estimate_marker_areas(points, axes, center=(0.0, 0.0, 0.0), *,
                          samples_per_point=1024, replicates=4, seed=0):
    """Estimate physical surface Voronoi areas with independent area samples.

    A cell contains the surface locations nearest to its marker in Euclidean
    distance (not geodesic distance). New scrambled Sobol samples, distinct from
    relaxation samples, measure each cell's area fraction. Replicate variation
    estimates quadrature error; it is not an exact confidence bound.

    Returns (areas, area_standard_errors, report), ordered like input points.
    """
    points = np.asarray(points, dtype=float)
    axes = np.asarray(axes, dtype=float)
    center = np.asarray(center, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not len(points) or not np.all(np.isfinite(points)):
        raise ValueError("points must be a nonempty finite (N,3) array")
    if axes.shape != (3,) or not np.all(np.isfinite(axes)) or np.any(axes <= 0):
        raise ValueError("axes must contain three finite positive semi-axis lengths")
    if center.shape != (3,) or not np.all(np.isfinite(center)):
        raise ValueError("center must contain three finite coordinates")
    if not isinstance(samples_per_point, (int, np.integer)) or samples_per_point < 8:
        raise ValueError("area samples_per_point must be an integer >= 8")
    if not isinstance(replicates, (int, np.integer)) or replicates < 2:
        raise ValueError("area replicates must be an integer >= 2")
    if not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    if len(np.unique(points, axis=0)) != len(points):
        raise ValueError("area estimation requires distinct marker positions")

    scale = axes.max()
    unit_axes = axes / scale
    local_points = (points - center) / scale
    if np.max(np.abs(np.sum((local_points / unit_axes)**2, axis=1) - 1)) > 1e-8:
        raise ValueError("points must lie on the specified ellipsoid")
    # Carlson's symmetric elliptic integral gives the full smooth surface area:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprg.html
    total_area = float(4 * np.pi * np.prod(unit_axes)
                       * elliprg(*(1 / unit_axes**2)) * scale**2)
    tree = cKDTree(local_points)
    count = len(points)
    sample_count = max(65536, samples_per_point * count)
    per_replicate = (sample_count + replicates - 1) // replicates
    replicate_areas = []
    for replicate in range(replicates):
        hits = np.zeros(count, dtype=np.int64)
        for chunk, start in enumerate(range(0, per_replicate, 262144)):
            chunk_seed = int(np.random.SeedSequence([seed, 10000, replicate, chunk])
                             .generate_state(1)[0])
            samples = _surface_samples(min(262144, per_replicate - start), unit_axes,
                                       chunk_seed)
            _, labels = tree.query(samples)
            hits += np.bincount(labels, minlength=count)
        replicate_areas.append(total_area * (hits / per_replicate))
    replicate_areas = np.asarray(replicate_areas)
    areas = replicate_areas.mean(axis=0)
    standard_errors = replicate_areas.std(axis=0, ddof=1) / np.sqrt(replicates)
    uniform = total_area / count
    relative = areas / uniform - 1
    histogram_counts, histogram_edges = np.histogram(areas / uniform, bins="auto")
    report = {
        "definition": "Surface Voronoi cells defined by nearest Euclidean marker; area estimated by independent area-uniform Sobol samples.",
        "area_units": "coordinate_units_squared",
        "axes": axes.tolist(), "center": center.tolist(), "marker_count": count,
        "total_surface_area": total_area, "total_area_method": "Carlson elliprg",
        "uniform_area": uniform, "estimated_area_sum": float(areas.sum()),
        "area_min": float(areas.min()), "area_mean": float(areas.mean()),
        "area_max": float(areas.max()), "area_rms": float(np.sqrt(np.mean(areas**2))),
        "area_std": float(areas.std()), "area_cv": float(areas.std() / uniform),
        "relative_rms_deviation": float(np.sqrt(np.mean(relative**2))),
        "mean_absolute_relative_deviation": float(np.mean(np.abs(relative))),
        "max_absolute_relative_deviation": float(np.max(np.abs(relative))),
        "area_ratio_percentiles": {str(q): float(np.percentile(areas / uniform, q))
                                   for q in (0, 5, 25, 50, 75, 95, 100)},
        "fraction_within_relative_tolerance": {
            str(tolerance): float(np.mean(np.abs(relative) <= tolerance))
            for tolerance in (0.05, 0.1, 0.2)},
        "histogram": {"area_ratio_edges": histogram_edges.tolist(),
                      "marker_counts": histogram_counts.tolist()},
        "zero_sample_area_count": int(np.count_nonzero(areas == 0)),
        "quadrature": {
            "total_samples": int(per_replicate * replicates),
            "samples_per_replicate": int(per_replicate), "replicates": int(replicates),
            "seed": int(seed), "independent_of_relaxation_samples": True,
            "relative_standard_error_rms": float(np.sqrt(np.mean(standard_errors**2)) / uniform),
            "relative_standard_error_max": float(standard_errors.max() / uniform),
            "note": "Replicate-based estimates of quadrature error, not exact confidence bounds. Increase samples to check convergence.",
        },
    }
    return areas, standard_errors, report


def write_area_reports(output, points, areas, standard_errors, report, *, plot=False):
    """Write a marker-indexed CSV, summary JSON and optional distribution plot."""
    csv_path = output.with_suffix(".areas.csv")
    json_path = output.with_suffix(".areas.json")
    uniform = report["uniform_area"]
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("id", "x", "y", "z", "area", "uniform_area", "area_ratio",
                         "relative_deviation_percent", "area_standard_error"))
        for mid, (point, area, error) in enumerate(zip(points, areas, standard_errors)):
            writer.writerow((mid, *point, area, uniform, area / uniform,
                             100 * (area / uniform - 1), error))
    report = dict(report, coordinate_file=str(output), area_file=str(csv_path))
    json_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"Wrote marker areas to {csv_path}; summary to {json_path}")
    print(f"Area: total={report['total_surface_area']:.8g}, uniform S/N={uniform:.8g}, "
          f"min={areas.min():.8g}, max={areas.max():.8g}")
    print(f"Area uniformity: relative RMS (=CV)={report['relative_rms_deviation']:.2%}, "
          f"max deviation={report['max_absolute_relative_deviation']:.2%}; "
          f"quadrature standard-error RMS={report['quadrature']['relative_standard_error_rms']:.2%}")
    if report["zero_sample_area_count"]:
        print("WARNING: some cells received no area samples; increase --area-samples-per-point.")
    if plot:
        _plot_areas(output.with_suffix(".areas.png"), points, areas, report)


def _plot_areas(path, points, areas, report):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    ratios = areas / report["uniform_area"]
    deviations = 100 * (ratios - 1)
    fig = plt.figure(figsize=(10.5, 4.5), constrained_layout=True)
    histogram = fig.add_subplot(121)
    histogram.hist(ratios, bins=report["histogram"]["area_ratio_edges"], edgecolor="white")
    histogram.axvline(1, color="black", linestyle="--", label="Uniform area S/N")
    histogram.set(xlabel="Marker area / uniform area", ylabel="Number of markers",
                  title=f"Area CV = {report['area_cv']:.2%}")
    histogram.legend()
    spatial = fig.add_subplot(122, projection="3d")
    limit = max(1.0, float(np.max(np.abs(deviations))))
    scatter = spatial.scatter(*points.T, c=deviations, cmap="coolwarm", vmin=-limit,
                              vmax=limit, s=12, depthshade=False)
    spatial.set(xlabel="X", ylabel="Y", zlabel="Z", title="Local area deviation")
    spatial.set_box_aspect(report["axes"])
    for axis in (spatial.xaxis, spatial.yaxis, spatial.zaxis):
        axis.set_major_locator(MaxNLocator(nbins=3))
    fig.colorbar(scatter, ax=spatial, shrink=0.65, label="Deviation from S/N (%)")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"Wrote area distribution plot to {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-n", "--num-points", type=int, required=True,
                        help="exact number of Lagrangian markers")
    parser.add_argument("--axes", type=float, nargs=3, required=True,
                        metavar=("A", "B", "C"),
                        help="x/y/z semi-axis lengths in the supplied order, not diameters")
    parser.add_argument("--center", type=float, nargs=3, default=(0.0, 0.0, 0.0),
                        metavar=("X", "Y", "Z"), help="default: 0 0 0 (body coordinates)")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="destination coordinate text file")
    parser.add_argument("--method", choices=("riesz3", "centroid"), default="riesz3",
                        help="placement method (default: riesz3, energy sum(1/r**3))")
    parser.add_argument("--iterations", type=int, default=None,
                        help="iteration limit (default: E3 3000, centroid 80; 0 disables)")
    parser.add_argument("--samples-per-point", type=int, default=64,
                        help="centroid-method quadrature samples per marker (default: 64; minimum: 8)")
    parser.add_argument("--seed", type=int, default=0, help="reproducible sampling seed (default: 0)")
    parser.add_argument("--area-samples-per-point", type=int, default=1024,
                        help="independent area samples per marker, total across replicates (default: 1024)")
    parser.add_argument("--area-replicates", type=int, default=4,
                        help="independent area sampling replicates, minimum 2 (default: 4)")
    parser.add_argument("--plot-areas", action="store_true",
                        help="also save an area histogram and spatial deviation plot (requires Matplotlib)")
    args = parser.parse_args(argv)
    axes = np.asarray(args.axes, dtype=float)
    try:
        points, optimization = generate_ellipsoid_points(
            args.num_points, axes, args.center, iterations=args.iterations,
            samples_per_point=args.samples_per_point, seed=args.seed, method=args.method,
            return_diagnostics=True,
        )
        areas, standard_errors, report = estimate_marker_areas(
            points, axes, args.center, samples_per_point=args.area_samples_per_point,
            replicates=args.area_replicates, seed=args.seed,
        )
    except ValueError as exc:
        parser.error(str(exc))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # 17 significant digits retain double precision, including small geometries.
    np.savetxt(args.output, points, fmt="%.17g", header="X Y Z", comments="# ")
    print(f"Wrote {len(points)} points to {args.output}")
    print(f"Semi-axes (x, y, z): {axes.tolist()}")
    print(f"Placement: {args.method}, iterations={optimization['iterations']}, "
          f"converged={optimization['converged']}; {optimization['message']}")
    if not optimization["converged"]:
        print("WARNING: placement has not met its stopping criterion; inspect diagnostics or increase --iterations.")
    report["optimization"] = optimization
    if len(points) > 1:
        distances = cKDTree(points).query(points, k=2)[0][:, 1]
        print(f"Nearest-neighbor distance: min={distances.min():.8g}, "
              f"mean={distances.mean():.8g}, max={distances.max():.8g}, "
              f"CV={distances.std()/distances.mean():.2%}")
    write_area_reports(args.output, points, areas, standard_errors, report, plot=args.plot_areas)


if __name__ == "__main__":
    main()
