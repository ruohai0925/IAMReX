# SPDX-FileCopyrightText: 2026 IAMReX contributors
# SPDX-License-Identifier: BSD-3-Clause

"""MPI client that emulates the C++ side of the RKPM MPMD protocol."""

from __future__ import annotations

import numpy as np
from mpi4py import MPI

from common import build_small_fixture

from src.grid_index import build_stencil_indices


TAG_MARKER_COUNT = 200
TAG_POSITIONS = 201
TAG_WEIGHTS = 202
STENCIL_INDICES = 203


def expected_weights(positions, exchange):
    return np.asarray([
        np.arange(27, dtype=float)
        + 100.0 * marker_id
        + 1000.0 * exchange
        + positions[marker_id, 0]
        for marker_id in range(len(positions))
    ], dtype=np.float32)


if __name__ == "__main__":
    world = MPI.COMM_WORLD
    appnum = world.Get_attr(MPI.APPNUM)
    appnum_send = np.asarray([appnum], dtype=np.int32)
    appnums = np.empty(world.Get_size(), dtype=np.int32)
    world.Allgather([appnum_send, MPI.INT], [appnums, MPI.INT])
    local = world.Split(color=appnum)

    grid, initial_positions, _, _ = build_small_fixture()
    server_rank = world.Get_size() - 1
    for exchange in range(2):
        positions = initial_positions.copy()
        # Put marker 0 exactly on the lower x-face of its expected cell. The
        # double-precision protocol plus snap-to-grid-line indexing must keep it
        # on the high-index side of that face.
        positions[0, 0] = grid.prob_lo[0] + 4 * grid.dx[0]
        positions[:, 0] += exchange * 1.01 * grid.dx[0]
        positions = positions.astype(np.float64)
        marker_count = np.asarray([len(positions)], dtype=np.int32)

        world.Send(
            [marker_count, MPI.INT], dest=server_rank, tag=TAG_MARKER_COUNT
        )
        world.Send(
            [positions.reshape(-1), MPI.DOUBLE],
            dest=server_rank,
            tag=TAG_POSITIONS,
        )
        received = np.empty(len(positions) * 27, dtype=np.float32)
        stencil_indices = np.empty(len(positions) * 27 * 3, dtype=np.int32)
        world.Recv(
            [received, MPI.FLOAT], source=server_rank, tag=TAG_WEIGHTS
        )
        world.Recv(
            [stencil_indices, MPI.INT], source=server_rank, tag=STENCIL_INDICES
        )
        expected_indices = build_stencil_indices(positions, grid.prob_lo, grid.dx)
        stencil_indices = stencil_indices.reshape(len(positions), 27, 3)
        if not np.array_equal(stencil_indices, expected_indices):
            raise ValueError("stencil_indices is not what expected")

        expected = expected_weights(positions, exchange).reshape(-1)
        if not np.array_equal(received, expected):
            max_error = float(np.max(np.abs(received - expected)))
            raise AssertionError(
                f"MPMD exchange {exchange} changed weight values or order; "
                f"max_error={max_error}"
            )

    stop = np.asarray([0], dtype=np.int32)
    world.Send([stop, MPI.INT], dest=server_rank, tag=TAG_MARKER_COUNT)
    local.Free()
    print("MPMD_CLIENT_PASS", flush=True)
