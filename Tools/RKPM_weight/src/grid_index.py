# SPDX-FileCopyrightText: 2026 IAMReX contributors
# SPDX-License-Identifier: BSD-3-Clause

"""Shared conversion from physical coordinates to Eulerian cell indices."""

from __future__ import annotations

import numpy as np


CELL_INDEX_SNAP_ULPS = 16


def containing_cell_indices(positions, prob_lo, dx) -> np.ndarray:
    """Return cell indices using a snap-to-grid-line, then floor convention.

    Cells use half-open intervals, so an exact internal grid-line coordinate is
    assigned to the cell on its high-index side. Floating-point arithmetic can
    place that coordinate a few representable values below the integer cell
    coordinate; values within a small roundoff-scaled tolerance are therefore
    snapped to the nearest integer before ``floor`` is applied.
    """
    points = np.asarray(positions, dtype=np.float64)
    lower = np.asarray(prob_lo, dtype=np.float64)
    spacing = np.asarray(dx, dtype=np.float64)

    if points.ndim == 0 or points.shape[-1] != 3:
        raise ValueError(
            f"Positions must have a final dimension of size 3, got {points.shape}"
        )
    if lower.shape != (3,) or spacing.shape != (3,):
        raise ValueError(
            "prob_lo and dx must each contain exactly three coordinates"
        )
    if not (
        np.all(np.isfinite(points))
        and np.all(np.isfinite(lower))
        and np.all(np.isfinite(spacing))
    ):
        raise ValueError("Positions, prob_lo and dx must be finite")
    if np.any(spacing <= 0.0):
        raise ValueError(f"Grid spacing must be positive, got {spacing}")

    normalized = (points - lower) / spacing
    nearest_integer = np.rint(normalized)

    # Include the magnitudes of both operands before subtraction. This keeps
    # the tolerance meaningful when the domain origin is large relative to dx.
    scale = np.maximum(1.0, np.abs(normalized))
    scale = np.maximum(scale, np.abs(points / spacing))
    scale = np.maximum(scale, np.abs(lower / spacing))
    tolerance = CELL_INDEX_SNAP_ULPS * np.finfo(np.float64).eps * scale

    snapped = np.where(
        np.abs(normalized - nearest_integer) <= tolerance,
        nearest_integer,
        normalized,
    )
    return np.floor(snapped).astype(np.int64)

# inputs :
# positions (N, 3) : all markers' position in world frame
# prob_lo, dx
# outputs : 27 stencil ijk of each marker (N, 27, 3)
def build_stencil_indices(positions, prob_lo, dx):
    indices_ijk = containing_cell_indices(positions, prob_lo, dx)  # cell where markers locate
    if (indices_ijk.ndim != 2):
        raise ValueError("invalid indices_ijk dimension")
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
    return support_indices