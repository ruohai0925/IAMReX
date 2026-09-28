# SPDX-FileCopyrightText: 2026 IAMReX contributors
# SPDX-License-Identifier: BSD-3-Clause

"""Deterministic MPMD server process used by the protocol test."""

from __future__ import annotations

import numpy as np

from common import build_small_fixture

from src.mpmd_transport import serve_mpmd


class DeterministicSolver:
    """Return position- and exchange-dependent weights with an obvious order."""

    name = "deterministic-test"

    def __init__(self):
        self.exchange = 0

    def solve_array(self, positions, lag_map, stencil_indices=None):
        exchange = self.exchange
        self.exchange += 1
        marker_ids = tuple(sorted(lag_map))
        weights = np.asarray(
            [
                np.arange(27, dtype=float)
                + 100.0 * marker_id
                + 1000.0 * exchange
                + positions[marker_id, 0]
                for marker_id in marker_ids
            ]
        )
        return marker_ids, weights


if __name__ == "__main__":
    fixture_grid, _, _, fixture_lag = build_small_fixture()
    serve_mpmd(
        DeterministicSolver(),
        fixture_lag,
        fixture_grid,
        check_action="off",
    )
