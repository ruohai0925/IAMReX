#!/usr/bin/env python3
"""Unified entry point for RKPM/ML weight generation and MPMD service."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from src import mapping
from src.grid_index import containing_cell_indices
from src.weight_solver import (
    DEFAULT_BATCH_SIZE,
    GridSpec,
    build_solver,
    parse_inputs,
)


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = SCRIPT_DIR / "unit_test" / "fixtures" / "real_case"
DEFAULT_INPUTS = DEFAULT_DATA_DIR / "inputs.3d.flow_past_ellipsoid"
DEFAULT_ID_FILE = DEFAULT_DATA_DIR / "rkpm_mappings.id"
DEFAULT_LAG_FILE = DEFAULT_DATA_DIR / "rkpm_mappings.lag"
DEFAULT_OUTPUT_PREFIX = DEFAULT_DATA_DIR / "rkpm_mappings"
DEFAULT_MODEL_ROOT = SCRIPT_DIR.parents[1] / "Transolver_v2"
DEFAULT_MODEL_DIR = DEFAULT_MODEL_ROOT / "out_inv"
DEFAULT_MODEL_CODE = DEFAULT_MODEL_ROOT / "transolver_slim.py"

CONFIG_PATH_KEYS = {
    "inputs": "inputs",
    "geometry": "geometry",
    "id": "id_file",
    "lag": "lag_file",
    "output_prefix": "output_prefix",
    "model_dir": "model_dir",
    "model_code": "model_code",
}
CONFIG_VALUE_KEYS = {
    "solver": ("solver", str),
    "transport": ("transport", str),
    "batch_size": ("batch_size", int),
    "device": ("device", str),
    "body_frame": ("body_frame", None),
    "angle": ("angle", float),
    "axes": ("axes", lambda value: list(map(float, value.split()))),
    "check_action": ("check_action", str),
    "check_interval": ("check_interval", int),
    "sum_tolerance": ("sum_tolerance", float),
    "moment_tolerance": ("moment_tolerance", float),
}


def _unquote(value):
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        return value[1:-1]
    return value


def _parse_bool(value, key):
    normalized = _unquote(value).lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"Configuration key {key} expects a boolean, got {value!r}")


def _load_run_config(path):
    """Load Python-tool arguments from an AMReX-style ``key = value`` file."""
    path = Path(path).resolve()
    raw = parse_inputs(path)
    base = path.parent
    defaults = {}
    known_keys = set(CONFIG_PATH_KEYS) | set(CONFIG_VALUE_KEYS)

    for raw_key, raw_value in raw.items():
        key = raw_key.removeprefix("rkpm.")
        if key not in known_keys:
            raise ValueError(f"Unknown RKPM configuration key: {raw_key}")
        value = _unquote(raw_value)
        if key in CONFIG_PATH_KEYS:
            destination = CONFIG_PATH_KEYS[key]
            if value.lower() in {"", "none", "null"}:
                defaults[destination] = None
            else:
                configured_path = Path(value).expanduser()
                if not configured_path.is_absolute():
                    configured_path = base / configured_path
                defaults[destination] = configured_path.resolve()
            continue

        destination, converter = CONFIG_VALUE_KEYS[key]
        defaults[destination] = (
            _parse_bool(value, raw_key) if converter is None else converter(value)
        )

    print(f"[config] Loaded {path}")
    return defaults


def _first_float(params, key):
    if key not in params:
        raise KeyError(f"Missing parameter in inputs file: {key}")
    return float(params[key].split()[0])


def _body_center(inputs_path):
    params = parse_inputs(inputs_path)
    return np.asarray([
        _first_float(params, "particle_inputs.x"),
        _first_float(params, "particle_inputs.y"),
        _first_float(params, "particle_inputs.z"),
    ])


def _print_grid(grid: GridSpec):
    print(f"[inputs] prob_lo = {grid.prob_lo.tolist()}")
    print(f"[inputs] prob_hi = {grid.prob_hi.tolist()}")
    print(
        f"[inputs] amr.n_cell = {grid.n_cell.tolist()}, "
        f"amr.max_level = {grid.max_level}"
    )
    print(f"[inputs] dx_finest = {grid.dx.tolist()}")


def _save_mapping(id_map, lag_map, output_prefix):
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    mapping.save_mappings_txt(id_map, lag_map, str(output_prefix))
    print(f"[file] Wrote {output_prefix}.id and {output_prefix}.lag")


def _weights_in_support_order(solved_lag, all_S_I, grid):
    """Restore .lag-sorted weights to the original ``all_S_I`` point order."""
    result = []
    for marker_id, support in enumerate(all_S_I):
        by_cell = {
            (row["i"], row["j"], row["k"]): row["w"]
            for row in solved_lag[marker_id]
        }
        marker_weights = []
        for point in support[:, :3]:
            cell = tuple(
                int(value)
                for value in containing_cell_indices(
                    point, grid.prob_lo, grid.dx
                )
            )
            marker_weights.append(by_cell[cell])
        result.append(np.asarray(marker_weights))
    return result


def generate_from_geometry(args, grid, solver):
    """Build stencils from a point cloud and solve and write their weights."""
    # These modules load matplotlib, so import them only for offline generation.
    from src import SI_generated, error, visual

    center = _body_center(args.inputs) if args.body_frame else None
    if center is not None:
        print(f"[inputs] center (particle_inputs) = {center.tolist()}")

    (
        eulerian_points,
        _,
        lagrangian_points,
        nearest_grid_points,
        delta_I,
        eta_I,
        theta_I,
        all_S_I,
        V_lag,
    ) = SI_generated.generate_grid(
        grid.prob_lo,
        grid.prob_hi,
        grid.dx,
        args.geometry,
        center=center,
        angle=args.angle,
        axes=args.axes,
    )

    if center is None:
        center = (
            lagrangian_points.min(axis=0) + lagrangian_points.max(axis=0)
        ) / 2.0
        print(f"[inputs] center (bbox) = {center.tolist()}")

    # visual.visualize_results(
    #     lagrangian_points,
    #     nearest_grid_points,
    #     delta_I,
    #     eta_I,
    #     theta_I,
    #     all_S_I,
    #     target_idx=0,
    # )

    id_map = mapping.build_lagrangian_id_to_coord_map(lagrangian_points)
    # Stencil/eps data is backend-independent. Use zeros until the solver runs.
    zero_weights = [np.zeros(len(support)) for support in all_S_I]
    lag_skeleton = mapping.build_lag_to_eul_map(
        lagrangian_points,
        all_S_I,
        zero_weights,
        grid.prob_lo,
        grid.dx,
        V_lag,
    )

    started = time.perf_counter()
    weights = solver.solve(id_map, lag_skeleton)
    print(f"[{solver.name}] Solve completed in {time.perf_counter()-started:.4f} s")
    solved_lag = mapping.replace_mapping_weights(lag_skeleton, weights)

    weights_by_marker = _weights_in_support_order(solved_lag, all_S_I, grid)
    error.compute_error_volume(
        lagrangian_points, all_S_I, weights_by_marker, V_lag
    )
    error.compute_conservation_check(
        center,
        eulerian_points,
        lagrangian_points,
        all_S_I,
        weights_by_marker,
        V_lag,
    )
    _save_mapping(id_map, solved_lag, args.output_prefix)


def recalculate_mapping(args, solver):
    """Read positions and stencils from .id/.lag and recompute only weights."""
    id_map = mapping.load_id_map(args.id_file)
    lag_map = mapping.load_lag_map(args.lag_file)
    marker_ids = mapping.validate_mapping_ids(id_map, lag_map)
    print(
        f"[file] Loaded {len(marker_ids)} markers from "
        f"{args.id_file}, {args.lag_file}"
    )

    started = time.perf_counter()
    weights = solver.solve(id_map, lag_map)
    print(f"[{solver.name}] Solve completed in {time.perf_counter()-started:.4f} s")
    solved_lag = mapping.replace_mapping_weights(lag_map, weights)
    _save_mapping(id_map, solved_lag, args.output_prefix)


def run(args):
    grid = GridSpec.from_inputs(args.inputs)
    _print_grid(grid)
    # initialize solver
    solver = build_solver(
        args.solver,
        grid,
        model_dir=args.model_dir,
        model_code=args.model_code,
        batch_size=args.batch_size,
        device=args.device,
    )
    device_text = (
        f", device={solver.device_name}" if solver.name == "ml" else ""
    )
    print(
        f"[solver] backend={solver.name}, batch_size={solver.batch_size}"
        f"{device_text}"
    )

    if args.transport == "mpmd":
        if args.geometry is not None:
            raise ValueError("MPMD receives positions from C++; do not use --geometry")
        lag_map = mapping.load_lag_map(args.lag_file)
        from src.mpmd_transport import serve_mpmd

        serve_mpmd(
            solver,
            lag_map,
            grid,
            check_action=args.check_action,
            check_interval=args.check_interval,
            sum_tolerance=args.sum_tolerance,
            moment_tolerance=args.moment_tolerance,
        )
        return

    if args.geometry is not None:
        generate_from_geometry(args, grid, solver)
    else:
        recalculate_mapping(args, solver)


def parse_args(argv=None):
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config", type=Path)
    config_args, _ = config_parser.parse_known_args(argv)
    config_defaults = (
        _load_run_config(config_args.config) if config_args.config else {}
    )

    parser = argparse.ArgumentParser(
        description="Generate RKPM weight mappings or run an IAMReX MPMD weight service"
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="AMReX-style key=value file containing Python-tool arguments",
    )
    parser.add_argument(
        "--solver",
        choices=("rkpm", "ml"),
        default=config_defaults.get("solver", "rkpm"),
        help="Weight solver backend (default: rkpm)",
    )
    parser.add_argument(
        "--transport",
        choices=("file", "mpmd"),
        default=config_defaults.get("transport", "file"),
        help="Output transport: .id/.lag files or direct MPMD transfer (default: file)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=config_defaults.get("batch_size", DEFAULT_BATCH_SIZE),
        help=(
            "Markers per vectorized RKPM/ML batch; reduce this value to limit "
            f"memory use (default: {DEFAULT_BATCH_SIZE})"
        ),
    )
    parser.add_argument(
        "--device",
        choices=("cpu", "cuda", "auto"),
        default=config_defaults.get("device", "cpu"),
        help="PyTorch inference device for the ML backend (default: cpu)",
    )
    parser.add_argument(
        "--inputs",
        type=Path,
        default=config_defaults.get("inputs", DEFAULT_INPUTS),
        help=f"AMReX inputs file (default: {DEFAULT_INPUTS})",
    )
    parser.add_argument(
        "--geometry",
        type=Path,
        default=config_defaults.get("geometry"),
        help="Three-column point cloud for file mode; omit to reuse --id/--lag",
    )
    parser.add_argument(
        "--id",
        dest="id_file",
        type=Path,
        default=config_defaults.get("id_file", DEFAULT_ID_FILE),
        help="Existing marker-coordinate file used when recomputing file output",
    )
    parser.add_argument(
        "--lag",
        dest="lag_file",
        type=Path,
        default=config_defaults.get("lag_file", DEFAULT_LAG_FILE),
        help="Existing stencil file used by file-recompute and MPMD modes",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=config_defaults.get("output_prefix", DEFAULT_OUTPUT_PREFIX),
        help="Output prefix in file mode (default: rkpm_mappings)",
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=config_defaults.get("model_dir", DEFAULT_MODEL_DIR),
        help="Directory containing the ML checkpoint and normalization statistics",
    )
    parser.add_argument(
        "--model-code",
        type=Path,
        default=config_defaults.get("model_code", DEFAULT_MODEL_CODE),
        help="Python file defining TransolverSlim",
    )
    frame_group = parser.add_mutually_exclusive_group()
    frame_group.add_argument(
        "--body-frame",
        action="store_true",
        dest="body_frame",
        help="Treat --geometry as body-frame data and translate it using particle_inputs",
    )
    frame_group.add_argument(
        "--world-frame",
        action="store_false",
        dest="body_frame",
        help="Treat --geometry coordinates as world-frame data",
    )
    # Keep an unspecified frame distinct from an explicit world-frame choice.
    # Geometry generation needs this distinction because body-frame point clouds
    # must be translated by the particle center before their grid bounds are used.
    parser.set_defaults(body_frame=config_defaults.get("body_frame"))
    parser.add_argument(
        "--angle",
        type=float,
        default=config_defaults.get("angle", 0.0),
        help="Rotation of the --geometry point cloud around the z-axis in degrees",
    )
    parser.add_argument(
        "--axes",
        nargs=3,
        type=float,
        metavar=("A", "B", "C"),
        default=config_defaults.get("axes"),
        help="Semi-axes of the ellipsoid carrying the markers, in grid length units",
    )
    parser.add_argument(
        "--check-action",
        choices=("off", "warn", "abort"),
        default=config_defaults.get("check_action", "off"),
        help="MPMD weight-moment check behavior (default: off)",
    )
    parser.add_argument(
        "--check-interval",
        type=int,
        default=config_defaults.get("check_interval", 1),
        help="Check the first MPMD exchange and every N exchanges thereafter",
    )
    parser.add_argument(
        "--sum-tolerance",
        type=float,
        default=config_defaults.get("sum_tolerance", 1.0e-6),
        help="Maximum MPMD zeroth-moment residual",
    )
    parser.add_argument(
        "--moment-tolerance",
        type=float,
        default=config_defaults.get("moment_tolerance", 1.0e-6),
        help="Maximum MPMD first-moment residual in cell units",
    )
    args = parser.parse_args(argv)
    if args.solver not in {"rkpm", "ml"}:
        parser.error("--solver must be either rkpm or ml")
    if args.transport not in {"file", "mpmd"}:
        parser.error("--transport must be either file or mpmd")
    if args.axes is not None and (args.transport != "file" or args.geometry is None):
        parser.error("--axes requires --geometry and --transport file")
    if (
        args.transport == "file"
        and args.geometry is not None
        and args.body_frame is None
    ):
        parser.error(
            "file-mode geometry generation requires an explicit coordinate "
            "frame: use --body-frame for coordinates relative to the particle "
            "center, or --world-frame for coordinates already in the AMReX domain"
        )
    if args.batch_size <= 0:
        parser.error("--batch-size must be a positive integer")
    if args.device not in {"cpu", "cuda", "auto"}:
        parser.error("--device must be cpu, cuda, or auto")
    if args.check_action not in {"off", "warn", "abort"}:
        parser.error("--check-action must be off, warn, or abort")
    if args.check_interval <= 0:
        parser.error("--check-interval must be a positive integer")
    if args.sum_tolerance <= 0.0 or args.moment_tolerance <= 0.0:
        parser.error("moment-check tolerances must be positive")
    return args


if __name__ == "__main__":
    program_start = time.perf_counter()
    run(parse_args())
    print(f"Total runtime: {time.perf_counter()-program_start:.4f} s")
