# SPDX-FileCopyrightText: 2026 IAMReX contributors
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for vectorized 3x3x3 support-domain construction."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from src import SI_generated
from src.grid_index import containing_cell_indices


class SupportGenerationTests(unittest.TestCase):
    def test_generated_supports_are_complete_and_use_expected_cells(self):
        prob_lo = np.asarray([0.2, -0.4, 1.0])
        dx = np.asarray([0.1, 0.2, 0.25])
        prob_hi = prob_lo + np.asarray([16, 14, 12]) * dx
        points = np.asarray(
            [
                [0.75, 0.72, 2.38],
                [0.82, 0.61, 2.47],
                [0.69, 0.79, 2.31],
                [0.77, 0.68, 2.55],
            ]
        )

        with tempfile.TemporaryDirectory() as directory:
            geometry = Path(directory) / "markers.txt"
            np.savetxt(geometry, points)
            with patch.object(SI_generated.visual, "PointCloud"):
                generated = SI_generated.generate_grid(
                    prob_lo, prob_hi, dx, geometry_file=geometry
                )

        eulerian_points = generated[0]
        lagrangian_points = generated[2]
        nearest = generated[3]
        supports = generated[7]
        containing = containing_cell_indices(lagrangian_points, prob_lo, dx)

        np.testing.assert_allclose(
            nearest, prob_lo + (containing + 0.5) * dx, rtol=0.0, atol=0.0
        )
        expected_offsets = {
            (i, j, k)
            for i in (-1, 0, 1)
            for j in (-1, 0, 1)
            for k in (-1, 0, 1)
        }
        all_cells = []
        for marker_index, support in enumerate(supports):
            self.assertEqual(support.shape, (27, 4))
            support_cells = containing_cell_indices(support[:, :3], prob_lo, dx)
            offsets = {
                tuple(cell - containing[marker_index]) for cell in support_cells
            }
            self.assertEqual(offsets, expected_offsets)
            np.testing.assert_allclose(support[:, 3], np.prod(dx))
            all_cells.append(support_cells)

        affected_cells = containing_cell_indices(eulerian_points, prob_lo, dx)
        expected_affected = np.unique(np.concatenate(all_cells), axis=0)
        np.testing.assert_array_equal(affected_cells, expected_affected)

    def test_marker_volume_is_invariant_under_rotation(self):
        cloud_axes = (0.38, 0.24, 0.15)
        positive_points = np.diag(cloud_axes)
        points = np.concatenate((positive_points, -positive_points))
        prob_lo, prob_hi = [-2.0] * 3, [2.0] * 3
        dx = [0.1, 0.125, 0.2]
        explicit_axes = (0.4, 0.25, 0.2)

        with tempfile.TemporaryDirectory() as directory:
            geometry = Path(directory) / "markers.txt"
            np.savetxt(geometry, points)
            for axes, expected_axes in ((None, cloud_axes), (explicit_axes, explicit_axes)):
                surface_area = SI_generated.ellipsoid_area_approx(*expected_axes)
                expected_volume = surface_area * min(dx) / len(points)
                for angle in (0.0, 45.0, 90.0):
                    with self.subTest(axes=axes, angle=angle):
                        generated = SI_generated.generate_grid(
                            prob_lo,
                            prob_hi,
                            dx,
                            geometry_file=geometry,
                            center=[0.2, 0.1, -0.05],
                            angle=angle,
                            axes=axes,
                        )
                        np.testing.assert_allclose(generated[-1], expected_volume, rtol=1e-14)


if __name__ == "__main__":
    unittest.main()
