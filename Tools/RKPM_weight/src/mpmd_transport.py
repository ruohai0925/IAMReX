# SPDX-FileCopyrightText: 2026 IAMReX contributors
# SPDX-License-Identifier: BSD-3-Clause

"""MPMD transport loop between IAMReX C++ and an RKPM weight backend."""

from __future__ import annotations

import time

import numpy as np

from .weight_solver import (
    fixed_stencil_center_indices,
    positions_as_array,
    stencil_centers,
    validate_stencils,
)

from .grid_index import (
    build_stencil_indices
)
from .validation import compute_weight_moment_errors, weights_as_array


TAG_MARKER_COUNT = 200
TAG_POSITIONS = 201
TAG_WEIGHTS = 202
STENCIL_INDICES = 203
RKPM_STENCIL_SIZE = 27


def serve_mpmd(
    solver,
    lag_map,
    grid,
    *,
    check_action="off",
    check_interval=1,
    sum_tolerance=1.0e-6,
    moment_tolerance=1.0e-6,
) -> None:
    """Receive positions, solve and return weights until C++ sends ``Nm=0``."""
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("MPMD transport requires mpi4py") from exc

    world = MPI.COMM_WORLD
    rank = world.Get_rank()
    size = world.Get_size()
    appnum = world.Get_attr(MPI.APPNUM)
    print(
        f"[rkpm-server] world rank {rank}/{size}, appnum={appnum}, "
        f"solver={solver.name}",
        flush=True,
    )
    if appnum != 1:
        raise RuntimeError("The RKPM server must be the second mpiexec app context")

    # Match the collective call order and data type in C++ MPMD::Initialize.
    appnum_send = np.asarray([appnum], dtype=np.int32)
    appnums = np.empty(size, dtype=np.int32)
    world.Allgather([appnum_send, MPI.INT], [appnums, MPI.INT])
    local = world.Split(color=appnum)

    try:
        if rank != size - 1:
            raise RuntimeError(
                "The C++ protocol uses the last world rank as the RKPM server; "
                "launch exactly one MPI rank for the Python server"
            )
        if local.Get_size() != 1:
            raise RuntimeError("The Python RKPM server supports exactly one MPI rank")

        prepared = solver.prepare(lag_map) if hasattr(solver, "prepare") else None
        if prepared is None:
            marker_ids = validate_stencils(
                lag_map, expected_size=RKPM_STENCIL_SIZE
            )
            cell_indices = None
            centers = stencil_centers(lag_map, grid, marker_ids)
        else:
            marker_ids = list(prepared.marker_ids)
            if (
                len(prepared.groups) != 1
                or prepared.groups[0].cell_indices.shape[1] != RKPM_STENCIL_SIZE
            ):
                raise ValueError("MPMD requires one rectangular 27-point stencil")
            cell_indices = prepared.groups[0].cell_indices
            centers = grid.prob_lo + (cell_indices + 0.5) * grid.dx
        expected_markers = len(marker_ids)
        fixed_stencil_center_indices(
            lag_map, marker_ids, cell_indices=cell_indices
        )
        cfd_root = 0
        step = 0
        print(
            f"[rkpm-server] Loaded {expected_markers} markers, "
            f"dx={grid.dx.tolist()}, check={check_action}, "
            f"check_interval={check_interval}",
            flush=True,
        )

        while True:
            start = time.perf_counter()
            count = np.empty(1, dtype=np.int32)
            world.Recv([count, MPI.INT], source=cfd_root, tag=TAG_MARKER_COUNT)
            marker_count = int(count[0])
            if marker_count == 0:
                print("[rkpm-server] Received shutdown signal", flush=True)
                break
            if marker_count != expected_markers:
                raise ValueError(
                    f"C++ sent {marker_count} markers, but the lag file contains "
                    f"{expected_markers}"
                )
            # Receive the flat C++ position array.
            position_buffer = np.empty(marker_count * 3, dtype=np.float64)
            world.Recv(
                [position_buffer, MPI.DOUBLE], source=cfd_root, tag=TAG_POSITIONS
            )
            # Reshape the double-precision wire buffer to (Nm, 3).
            positions = positions_as_array(
                position_buffer.reshape(marker_count, 3), marker_ids
            )
            received = time.perf_counter()
            stencil_indices = build_stencil_indices(positions, grid.prob_lo, grid.dx)
            centers = grid.prob_lo + (stencil_indices + 0.5) * grid.dx

            stencil_ready_at = time.perf_counter()
            if hasattr(solver, "solve_array"):
                solved_ids, weight_array = solver.solve_array(positions, lag_map, stencil_indices)
                if tuple(solved_ids) != tuple(marker_ids):
                    raise ValueError("Solver output marker IDs changed during MPMD")
            else:
                weights = solver.solve(positions, lag_map)
                weight_array = weights_as_array(weights, marker_ids)
            solved = time.perf_counter()
            weight_array = np.asarray(weight_array, dtype=np.float32)
            if weight_array.shape != (marker_count, RKPM_STENCIL_SIZE):
                raise ValueError(
                    "MPMD weights must have shape "
                    f"({marker_count}, {RKPM_STENCIL_SIZE}), "
                    f"got {weight_array.shape}"
                )
            packed = time.perf_counter()
            step += 1
            check_due = (
                check_action != "off" and (step - 1) % check_interval == 0
            )
            check_text = ""
            if check_due:
                checked_at = time.perf_counter()
                errors = compute_weight_moment_errors(
                    positions, centers, weight_array, grid.dx
                )
                check_seconds = time.perf_counter() - checked_at
                check_text = (
                    f", check={check_seconds:.4f}s, "
                    f"sum_error={errors.sum_error:.3e}, "
                    f"moment_error={errors.moment_error:.3e} cells"
                )
                if not errors.passes(sum_tolerance, moment_tolerance):
                    message = (
                        f"Weight moment check failed at step {step}: "
                        f"sum_error={errors.sum_error:.3e} "
                        f"(tolerance {sum_tolerance:.3e}), "
                        f"moment_error={errors.moment_error:.3e} cells "
                        f"(tolerance {moment_tolerance:.3e})"
                    )
                    if check_action == "abort":
                        raise ValueError(message)
                    print(f"[rkpm-server] WARNING: {message}", flush=True)

            flat_weights = weight_array.reshape(-1)
            ready_to_send = time.perf_counter()

            world.Send(
                [flat_weights, MPI.FLOAT], dest=cfd_root, tag=TAG_WEIGHTS
            )
            world.Send(
                [stencil_indices.astype(np.int32), MPI.INT], dest=cfd_root, tag=STENCIL_INDICES
            )
            sent = time.perf_counter()
            print(
                f"[rkpm-server] {marker_count} markers: "
                f"receive={received-start:.4f}s, "
                f"stencil={stencil_ready_at-received:.4f}s, "
                f"{solver.name}={solved-stencil_ready_at:.4f}s, "
                f"pack={packed-solved:.4f}s, send={sent-ready_to_send:.4f}s, "
                f"sum={float(flat_weights.sum()):.6f}"
                f"{check_text}",
                flush=True,
            )
    except Exception as exc:
        # Abort the job because otherwise C++ remains blocked in MPI_Recv.
        print(f"[rkpm-server] Fatal error: {exc}", flush=True)
        world.Abort(1)
        raise
    finally:
        local.Free()
