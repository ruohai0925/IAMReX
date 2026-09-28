"""Runtime backends for RKPM weights.

Solvers only convert ``marker positions + fixed Euler stencil`` into weights.
The caller selects file output or MPMD transport, allowing traditional RKPM
and ML to share validation and coordinate conversion.
"""

from __future__ import annotations

import importlib.util
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Sequence

import numpy as np

from . import window
from .grid_index import containing_cell_indices


DEFAULT_BATCH_SIZE = 4096


@dataclass(frozen=True)
class PackedStencilGroup:
    """Rectangular NumPy representation of markers with one stencil size."""

    marker_indices: np.ndarray
    marker_ids: tuple[int, ...]
    cell_indices: np.ndarray
    eps: np.ndarray


@dataclass(frozen=True)
class PackedStencils:
    """Static stencil metadata cached by a runtime weight solver."""

    marker_ids: tuple[int, ...]
    groups: tuple[PackedStencilGroup, ...]


@dataclass(frozen=True)
class GridSpec:
    prob_lo: np.ndarray
    prob_hi: np.ndarray
    n_cell: np.ndarray
    max_level: int
    dx: np.ndarray

    @classmethod
    def from_inputs(cls, path: Path | str) -> "GridSpec":
        params = parse_inputs(path)
        prob_lo = _get_vector(params, "geometry.prob_lo", float)
        prob_hi = _get_vector(params, "geometry.prob_hi", float)
        n_cell = _get_vector(params, "amr.n_cell", int)
        max_level = int(_get_tokens(params, "amr.max_level")[0])

        if not (len(prob_lo) == len(prob_hi) == len(n_cell) == 3):
            raise ValueError("RKPM_weight currently supports only 3D grid parameters")
        if np.any(n_cell <= 0) or max_level < 0:
            raise ValueError("amr.n_cell must be positive and amr.max_level nonnegative")

        dx = (prob_hi - prob_lo) / (n_cell * (2 ** max_level))
        return cls(prob_lo, prob_hi, n_cell, max_level, dx)


def parse_inputs(path: Path | str) -> Dict[str, str]:
    """Parse an AMReX ``key = value`` inputs file."""
    params = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            key, _, value = line.partition("=")
            params[key.strip()] = value.strip()
    return params


def _get_tokens(params: Mapping[str, str], key: str) -> Sequence[str]:
    if key not in params:
        raise KeyError(f"Missing parameter in inputs file: {key}")
    return params[key].split()


def _get_vector(params: Mapping[str, str], key: str, dtype) -> np.ndarray:
    return np.asarray([dtype(value) for value in _get_tokens(params, key)])


def validate_stencils(lag_map, *, expected_size: int | None = None) -> list[int]:
    """Validate marker ID continuity, stencil sizes and grid-index fields."""
    marker_ids = sorted(lag_map)
    if marker_ids != list(range(len(marker_ids))):
        raise ValueError("Marker IDs must be contiguous from zero to match C++ ordering")

    for marker_id in marker_ids:
        rows = lag_map[marker_id]
        if not rows:
            raise ValueError(f"Marker {marker_id} has an empty stencil")
        if expected_size is not None and len(rows) != expected_size:
            raise ValueError(
                f"MPMD requires exactly {expected_size} stencil points per marker; "
                f"marker {marker_id} has {len(rows)}"
            )
        for row in rows:
            for key in ("i", "j", "k"):
                if key not in row:
                    raise ValueError(f"Marker {marker_id} stencil is missing {key}")
    return marker_ids


def pack_stencils(lag_map, *, expected_size: int | None = None) -> PackedStencils:
    """Validate and pack dictionary-based stencils into rectangular groups."""
    marker_ids = tuple(validate_stencils(lag_map, expected_size=expected_size))
    grouped_marker_indices = {}
    for marker_index, marker_id in enumerate(marker_ids):
        grouped_marker_indices.setdefault(len(lag_map[marker_id]), []).append(
            marker_index
        )

    groups = []
    for marker_indices_list in grouped_marker_indices.values():
        group_ids = tuple(marker_ids[index] for index in marker_indices_list)
        stencil_size = len(lag_map[group_ids[0]])
        cell_indices = np.empty(
            (len(group_ids), stencil_size, 3), dtype=np.int64
        )
        eps = np.empty((len(group_ids), stencil_size), dtype=float)
        # Fill preallocated arrays one marker at a time. This avoids constructing
        # a second, potentially very large, nested Python list during MPMD setup.
        for group_index, marker_id in enumerate(group_ids):
            rows = lag_map[marker_id]
            cell_indices[group_index] = [
                (row["i"], row["j"], row["k"]) for row in rows
            ]
            eps[group_index] = [row.get("eps", np.nan) for row in rows]
        groups.append(
            PackedStencilGroup(
                marker_indices=np.asarray(marker_indices_list, dtype=int),
                marker_ids=group_ids,
                cell_indices=cell_indices,
                eps=eps,
            )
        )
    return PackedStencils(marker_ids=marker_ids, groups=tuple(groups))


def _validate_batch_size(batch_size: int) -> int:
    batch_size = int(batch_size)
    if batch_size <= 0:
        raise ValueError(f"batch_size must be positive, got {batch_size}")
    return batch_size


def _batch_ranges(item_count: int, batch_size: int):
    for start in range(0, item_count, batch_size):
        yield start, min(start + batch_size, item_count)


def positions_as_array(positions, marker_ids: Sequence[int]) -> np.ndarray:
    """Convert an ID mapping or array to an ID-ordered ``[N,3]`` array."""
    if isinstance(positions, Mapping):
        missing = set(marker_ids) - set(positions)
        extra = set(positions) - set(marker_ids)
        if missing or extra:
            raise ValueError(
                f"Position and stencil IDs differ; missing {sorted(missing)[:5]}, "
                f"extra {sorted(extra)[:5]}"
            )
        array = np.asarray(
            [positions[marker_id] for marker_id in marker_ids], dtype=float
        )
    else:
        array = np.asarray(positions, dtype=float)

    if array.shape != (len(marker_ids), 3):
        raise ValueError(
            f"Marker positions must have shape ({len(marker_ids)}, 3), got {array.shape}"
        )
    if not np.all(np.isfinite(array)):
        raise ValueError("Marker positions contain NaN or Inf")
    return array


def stencil_centers(lag_map, grid: GridSpec, marker_ids: Sequence[int]):
    """Construct world-coordinate cell centers from global cell indices."""
    try:
        indices = np.asarray(
            [
                [[row["i"], row["j"], row["k"]] for row in lag_map[marker_id]]
                for marker_id in marker_ids
            ],
            dtype=float,
        )
    except ValueError as exc:
        raise ValueError(
            "Stencil centers require equal stencil sizes for all markers"
        ) from exc
    return grid.prob_lo + (indices + 0.5) * grid.dx


def fixed_stencil_center_indices(
    lag_map,
    marker_ids: Sequence[int],
    *,
    cell_indices=None,
) -> np.ndarray:
    """Validate complete 3x3x3 stencils and return their center-cell indices."""
    if cell_indices is None:
        for marker_id in marker_ids:
            if len(lag_map[marker_id]) != 27:
                raise ValueError(
                    f"Marker {marker_id} is not a complete 3x3x3 stencil"
                )
        indices = pack_stencils(lag_map).groups[0].cell_indices
    else:
        indices = np.asarray(cell_indices, dtype=np.int64)
        if indices.shape != (len(marker_ids), 27, 3):
            raise ValueError(
                "Fixed-stencil indices must have shape "
                f"({len(marker_ids)}, 27, 3), got {indices.shape}"
            )
    centers = np.median(indices, axis=1).astype(np.int64)
    relative = indices - centers[:, None, :]
    inside = np.all((relative >= -1) & (relative <= 1), axis=(1, 2))
    codes = (
        (relative[..., 0] + 1) * 9
        + (relative[..., 1] + 1) * 3
        + relative[..., 2]
        + 1
    )
    complete = inside & np.all(
        np.sort(codes, axis=1) == np.arange(27), axis=1
    )
    if not np.all(complete):
        marker_index = int(np.flatnonzero(~complete)[0])
        raise ValueError(
            f"Marker {marker_ids[marker_index]} is not a complete 3x3x3 stencil"
        )
    return centers


def validate_fixed_stencils(
    positions,
    lag_map,
    grid: GridSpec,
    marker_ids,
    *,
    expected_centers=None,
):
    """Ensure each fixed 3x3x3 MPMD stencil still encloses its marker.

    The current C++ protocol returns weights but not ``i/j/k``. Once a marker
    crosses a cell boundary, the old stencil is invalid. Fail explicitly instead
    of silently applying weights to the wrong cells.
    """
    containing = containing_cell_indices(positions, grid.prob_lo, grid.dx)
    if expected_centers is None:
        expected_centers = fixed_stencil_center_indices(lag_map, marker_ids)
    else:
        expected_centers = np.asarray(expected_centers, dtype=np.int64)
        if expected_centers.shape != (len(marker_ids), 3):
            raise ValueError(
                "Expected fixed-stencil centers must have shape "
                f"({len(marker_ids)}, 3), got {expected_centers.shape}"
            )

    mismatched = np.any(containing != expected_centers, axis=1)
    if np.any(mismatched):
        row_index = int(np.flatnonzero(mismatched)[0])
        marker_id = marker_ids[row_index]
        raise ValueError(
            f"Marker {marker_id} left the center cell of its fixed stencil: "
            f"current cell {containing[row_index].tolist()}, stencil center "
            f"{expected_centers[row_index].tolist()}. The current MPMD protocol "
            "transfers weights only and cannot update i/j/k."
        )


class RKPMSolver:
    name = "rkpm"

    def __init__(self, grid: GridSpec, batch_size: int = DEFAULT_BATCH_SIZE):
        self.grid = grid
        self.batch_size = _validate_batch_size(batch_size)
        self._cached_lag_map = None
        self._packed = None
        self._lagrangian_volumes = None

    def prepare(self, lag_map) -> PackedStencils:
        """Pack and validate static stencil data once for repeated MPMD solves."""
        if self._cached_lag_map is lag_map:
            return self._packed

        packed = pack_stencils(lag_map)
        cell_volume = float(np.prod(self.grid.dx))
        lagrangian_volumes = []
        for group in packed.groups:
            eps = group.eps
            valid = np.all(np.isfinite(eps) & (eps > 0.0), axis=1)
            if not np.all(valid):
                bad = int(np.flatnonzero(~valid)[0])
                raise ValueError(
                    f"Marker {group.marker_ids[bad]} has invalid eps values; "
                    "cannot recover V_lag"
                )
            consistent = np.all(
                np.isclose(eps, eps[:, :1], rtol=1.0e-12, atol=0.0), axis=1
            )
            if not np.all(consistent):
                bad = int(np.flatnonzero(~consistent)[0])
                raise ValueError(
                    f"Marker {group.marker_ids[bad]} has inconsistent eps values "
                    "in its stencil"
                )
            lagrangian_volumes.append(eps[:, 0] * cell_volume)

        self._cached_lag_map = lag_map
        self._packed = packed
        self._lagrangian_volumes = tuple(lagrangian_volumes)
        return packed

    def _solve_rows(self, positions, lag_map, stencil_indices=None):
        packed = self.prepare(lag_map)
        marker_ids = packed.marker_ids
        marker_positions = positions_as_array(positions, marker_ids)
        if stencil_indices is not None:
            if (len(packed.groups) != 1):
                raise ValueError("unsupported multiple packed groups")
            if (stencil_indices.shape != (len(marker_ids), 27, 3)):
                raise ValueError("invalid stencil indices dimension")
            if (packed.groups[0].cell_indices.shape[1] != 27):
                raise ValueError("lag_map should have 27 stencils")

        # The .lag file stores eps = V_lag / Delta_V. Recover the physical cell
        # volume from the finest-grid spacing and then recover V_lag per marker.
        # Do not use row["Vcell"] here: in the current C++ mapping contract that
        # field is an interpolation/spreading multiplier fixed at 1.0, not the
        # physical Eulerian cell volume used by the RKPM moment equations.
        cell_volume = float(np.prod(self.grid.dx))
        scales = self.grid.dx * 1.001
        rectangular = len(packed.groups) == 1
        solved_by_marker = (
            np.empty(
                (len(marker_ids), packed.groups[0].cell_indices.shape[1]),
                dtype=float,
            )
            if rectangular
            else [None] * len(marker_ids)
        )

        for group, lagrangian_volumes in zip(
            packed.groups, self._lagrangian_volumes
        ):
            for start, stop in _batch_ranges(len(group.marker_ids), self.batch_size):
                global_indices = group.marker_indices[start:stop]
                if stencil_indices is None:
                    indices = group.cell_indices
                    centers = self.grid.prob_lo + (
                        indices[start:stop] + 0.5
                    ) * self.grid.dx
                else:
                    indices = stencil_indices
                    centers = self.grid.prob_lo + (
                                indices[global_indices] + 0.5
                            ) * self.grid.dx
                support_domains = np.empty(
                    (stop - start, centers.shape[1], 4), dtype=float
                )
                support_domains[..., :3] = centers
                support_domains[..., 3] = cell_volume
                solved = window.compute_modified_window_functions_batch(
                    support_domains,
                    marker_positions[global_indices],
                    scales[0],
                    scales[1],
                    scales[2],
                    lagrangian_volumes[start:stop],
                )
                if rectangular:
                    solved_by_marker[global_indices] = solved
                else:
                    for local_index, marker_index in enumerate(global_indices):
                        solved_by_marker[int(marker_index)] = solved[local_index]

        return marker_ids, solved_by_marker

    def solve_array(self, positions, lag_map, stencil_indices=None):
        """Return ID-ordered weights directly as a rectangular NumPy array."""
        marker_ids, solved_by_marker = self._solve_rows(positions, lag_map, stencil_indices)
        if not isinstance(solved_by_marker, np.ndarray):
            raise ValueError(
                "Array output requires equal stencil sizes for all markers"
            )
        return marker_ids, solved_by_marker

    def solve(self, positions, lag_map) -> Dict[int, np.ndarray]:
        marker_ids, solved_by_marker = self._solve_rows(positions, lag_map)
        return {
            marker_id: np.asarray(solved_by_marker[index], dtype=float)
            for index, marker_id in enumerate(marker_ids)
        }


class MLWeightSolver:
    name = "ml"

    def __init__(
        self,
        grid: GridSpec,
        model_dir: Path | str,
        model_code: Path | str,
        batch_size: int = DEFAULT_BATCH_SIZE,
        device: str = "cpu",
    ):
        self.grid = grid
        self.model_dir = Path(model_dir).resolve()
        self.model_code = Path(model_code).resolve()
        self.batch_size = _validate_batch_size(batch_size)
        self.device_name = device
        self._cached_lag_map = None
        self._packed = None
        self._load_model()

    def _load_model(self):
        try:
            import torch
        except ImportError as exc:
            raise RuntimeError("The ML solver requires PyTorch") from exc

        if not self.model_code.is_file():
            raise FileNotFoundError(
                f"Transolver model definition not found: {self.model_code}"
            )
        spec = importlib.util.spec_from_file_location(
            "rkpm_weight_transolver_slim", self.model_code
        )
        if spec is None or spec.loader is None:
            raise ImportError(f"Unable to load model definition: {self.model_code}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        if self.device_name == "auto":
            self.device_name = "cuda" if torch.cuda.is_available() else "cpu"
        if self.device_name == "cuda" and not torch.cuda.is_available():
            raise RuntimeError(
                "The ML solver requested CUDA, but PyTorch cannot access a CUDA device"
            )
        if self.device_name not in {"cpu", "cuda"}:
            raise ValueError(
                f"ML device must be cpu, cuda, or auto; got {self.device_name!r}"
            )
        self.device = torch.device(self.device_name)

        checkpoint_path = self.model_dir / "model_best.pt"
        stats_path = self.model_dir / "norm_stats.npz"
        checkpoint = torch.load(checkpoint_path, map_location="cpu")
        self.model = module.TransolverSlim(**checkpoint["cfg"])
        self.model.load_state_dict(checkpoint["state_dict"])
        self.model.to(self.device)
        self.model.eval()

        with np.load(stats_path) as stats:
            self.xm = stats["xm"].copy()
            self.xs = stats["xs"].copy()
            self.ym = stats["ym"].copy()
            self.ys = stats["ys"].copy()
        if np.any(self.xs == 0):
            raise ValueError(f"Model normalization parameter xs contains zero: {stats_path}")
        self.torch = torch

    def prepare(self, lag_map) -> PackedStencils:
        """Pack the fixed 27-point stencil once for repeated MPMD inference."""
        if self._cached_lag_map is not lag_map:
            self._packed = pack_stencils(lag_map, expected_size=27)
            self._cached_lag_map = lag_map
        return self._packed

    def solve_array(self, positions, lag_map, stencil_indices=None):
        """Return ID-ordered Transolver weights as one NumPy array."""
        packed = self.prepare(lag_map)
        if stencil_indices is not None:
            nb_marker = len(positions)
            if (stencil_indices.shape != (nb_marker,27,3)):
                raise ValueError("invalid stencil indices dimension")
        marker_ids = packed.marker_ids
        marker_positions = positions_as_array(positions, marker_ids)
        group = packed.groups[0]
        weights = np.empty((len(marker_ids), 27), dtype=float)

        with self.torch.inference_mode():
            for start, stop in _batch_ranges(len(marker_ids), self.batch_size):
                if stencil_indices is None:
                    indices = group.cell_indices
                else:
                    indices = stencil_indices
                centers = self.grid.prob_lo + (
                    indices[start:stop] + 0.5
                ) * self.grid.dx
                relative = np.asarray(
                    (centers - marker_positions[start:stop, None, :])
                    / self.grid.dx,
                    dtype=np.float32,
                )
                normalized = (relative - self.xm) / self.xs
                batch = self.torch.as_tensor(
                    normalized,
                    dtype=self.torch.float32,
                    device=self.device,
                )
                prediction = self.model(fx=batch, embedding=batch)
                predicted = (
                    prediction.detach().cpu().numpy()[..., 0] * self.ys + self.ym
                )
                expected_shape = (stop - start, 27)
                if predicted.shape != expected_shape:
                    raise ValueError(
                        f"ML output must have shape {expected_shape}, "
                        f"got {predicted.shape}"
                    )
                sums = predicted.sum(axis=1, keepdims=True)
                if not np.all(np.isfinite(predicted)) or np.any(
                    np.abs(sums) < 1.0e-12
                ):
                    raise ValueError(
                        "ML output contains NaN/Inf or a near-zero marker weight sum"
                    )
                weights[start:stop] = predicted / sums

        return marker_ids, weights

    def solve(self, positions, lag_map) -> Dict[int, np.ndarray]:
        marker_ids, weights = self.solve_array(positions, lag_map)
        return {
            marker_id: np.asarray(weights[index], dtype=float)
            for index, marker_id in enumerate(marker_ids)
        }


def build_solver(
    name: str,
    grid: GridSpec,
    *,
    model_dir=None,
    model_code=None,
    batch_size: int = DEFAULT_BATCH_SIZE,
    device: str = "cpu",
):
    if name == "rkpm":
        return RKPMSolver(grid, batch_size=batch_size)
    if name == "ml":
        if model_dir is None or model_code is None:
            raise ValueError("The ML solver requires --model-dir and --model-code")
        return MLWeightSolver(
            grid,
            model_dir,
            model_code,
            batch_size=batch_size,
            device=device,
        )
    raise ValueError(f"Unknown solver: {name}")
