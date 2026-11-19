# SPDX-FileCopyrightText: 2026 IAMReX contributors
# SPDX-License-Identifier: BSD-3-Clause

"""Check surface geometry, physical-area distribution and coordinate output."""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist

from generate_ellipsoid_points import (
    _riesz3_energy_gradient, estimate_marker_areas, generate_ellipsoid_points, main,
)


class EllipsoidPointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.points = generate_ellipsoid_points(512, (1.0, 1.0, 3.0))

    def test_count_surface_and_spacing(self):
        self.assertEqual(self.points.shape, (512, 3))
        np.testing.assert_allclose(np.sum((self.points / (1, 1, 3))**2, axis=1),
                                   1.0, rtol=0, atol=2e-14)
        nearest = cKDTree(self.points).query(self.points, k=2)[0][:, 1]
        self.assertGreater(nearest.min(), 0.65 * nearest.mean())
        self.assertLess(nearest.std() / nearest.mean(), 0.12)

    def test_polar_population_matches_surface_area_not_solid_angle(self):
        # For x=(sqrt(1-z²)cos(phi), sqrt(1-z²)sin(phi), 3z),
        # dS/dz/dphi=sqrt(9-8z²). A scaled uniform sphere would put 20%
        # of its points in |z|>0.8; the ellipsoid area fraction is 12.65%.
        density = lambda z: np.sqrt(9.0 - 8.0 * z * z)
        expected = quad(density, 0.8, 1.0)[0] / quad(density, 0.0, 1.0)[0]
        observed = np.mean(np.abs(self.points[:, 2] / 3.0) > 0.8)
        self.assertLess(abs(observed - expected), 0.025)

    def test_translation_and_reproducibility(self):
        arguments = dict(num_points=31, axes=(0.04, 0.06, 0.02), iterations=4, seed=7)
        body = generate_ellipsoid_points(**arguments)
        np.testing.assert_array_equal(body, generate_ellipsoid_points(**arguments))
        center = np.asarray((2.0, 2.25, 2.25))
        world = generate_ellipsoid_points(**arguments, center=center)
        np.testing.assert_allclose(world - center, body, rtol=0, atol=3e-16)

    def test_cli_coordinate_file_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "markers.txt"
            with contextlib.redirect_stdout(io.StringIO()):
                main(["-n", "17", "--axes", "0.04", "0.06", "0.04",
                      "--iterations", "0", "-o", str(output)])
            self.assertEqual(output.read_text().splitlines()[0], "# X Y Z")
            loaded = np.loadtxt(output)
            self.assertEqual(loaded.shape, (17, 3))
            np.testing.assert_allclose(np.sum((loaded / (.04, .06, .04))**2, axis=1),
                                       1.0, rtol=0, atol=2e-14)
            area_rows = np.genfromtxt(output.with_suffix(".areas.csv"), delimiter=",", names=True)
            np.testing.assert_array_equal(area_rows["id"], np.arange(17))
            np.testing.assert_array_equal(area_rows["x"], loaded[:, 0])
            report = json.loads(output.with_suffix(".areas.json").read_text())
            self.assertEqual(report["axes"], [.04, .06, .04])
            self.assertEqual(report["marker_count"], 17)
            self.assertEqual(report["optimization"]["method"], "riesz3")
            self.assertEqual(report["optimization"]["energy_exponent"], 3)
            self.assertAlmostEqual(float(area_rows["area"].sum()), report["total_surface_area"])

    def test_cli_preserves_axis_order_with_translation(self):
        center = np.asarray((2, -1, 0.5))
        with tempfile.TemporaryDirectory() as directory:
            for axes in ((.06, .02, .04), (.02, .04, .06)):
                with self.subTest(axes=axes):
                    output = Path(directory) / "markers.txt"
                    arguments = ["-n", "24", "--axes", *map(str, axes),
                                 "--center", *map(str, center), "--iterations", "2",
                                 "-o", str(output)]
                    with contextlib.redirect_stdout(io.StringIO()):
                        main(arguments)
                    points = np.loadtxt(output)
                    np.testing.assert_allclose(np.sum(((points-center)/axes)**2, axis=1),
                                               1, rtol=0, atol=3e-14)
                    report = json.loads(output.with_suffix(".areas.json").read_text())
                    self.assertEqual(report["axes"], list(axes))
                    rows = np.genfromtxt(output.with_suffix(".areas.csv"), delimiter=",", names=True)
                    np.testing.assert_array_equal(rows["x"], points[:, 0])
                    self.assertAlmostEqual(float(rows["area"].sum()), report["total_surface_area"])

    def test_cli_produces_mirror_pairs_around_translated_center(self):
        center = np.array([2, 2.25, 2.25])
        axes = np.array([0.04, 0.06, 0.04])
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "markers.txt"
            with contextlib.redirect_stdout(io.StringIO()):
                main(["-n", "32", "--axes", *map(str, axes), "--center", *map(str, center),
                      "--mirror-axis", "y", "--iterations", "10", "-o", str(output)])
            points = np.loadtxt(output)
        positive, negative = np.split(points - center, 2)
        self.assertEqual(points.shape, (32, 3))
        self.assertEqual(len(np.unique(points, axis=0)), 32)
        self.assertTrue(np.all(positive[:, 1] > 0))
        np.testing.assert_allclose(negative, positive * (1, -1, 1), atol=5e-16, rtol=0)
        np.testing.assert_allclose(np.sum(((points - center) / axes)**2, axis=1),
                                   1, atol=3e-14, rtol=0)

    def test_mirror_axis_works_with_centroid_and_requires_even_count(self):
        for normal_axis, mirror_axis in enumerate("xyz"):
            with self.subTest(mirror_axis=mirror_axis):
                points = generate_ellipsoid_points(24, (1, 1.5, 0.8), method="centroid",
                                                   mirror_axis=mirror_axis, iterations=3)
                positive, negative = np.split(points, 2)
                sign = np.ones(3)
                sign[normal_axis] = -1
                np.testing.assert_array_equal(negative, positive * sign)
                self.assertTrue(np.all(positive[:, normal_axis] > 0))
        with self.assertRaisesRegex(ValueError, "even number of points"):
            generate_ellipsoid_points(25, (1, 1.5, 0.8), mirror_axis="y")

    def test_e3_energy_and_gradient_against_pair_distances(self):
        axes = np.asarray((1.0, 1.5, 0.8))
        angles = np.random.default_rng(15).uniform(0.2, 1.4, 16)

        def reference(vector, mirror):
            theta, phi = vector.reshape(-1, 2).T
            points = np.column_stack((np.sin(theta)*np.cos(phi),
                                      np.sin(theta)*np.sin(phi), np.cos(theta))) * axes
            if mirror:
                points = np.concatenate((points, points * (1, 1, -1)))
            return np.mean(pdist(points)**-3)

        for mirror in (False, True):
            with self.subTest(mirror=mirror):
                energy, gradient = _riesz3_energy_gradient(angles, axes, block_size=3, mirror=mirror)
                self.assertAlmostEqual(energy, reference(angles, mirror))
                finite = np.empty_like(angles)
                for index in range(len(angles)):
                    plus, minus = angles.copy(), angles.copy()
                    plus[index] += 1e-6
                    minus[index] -= 1e-6
                    finite[index] = (reference(plus, mirror)-reference(minus, mirror))/2e-6
                self.assertLess(np.linalg.norm(gradient-finite)/np.linalg.norm(finite), 2e-7)

    def test_e3_reports_energy_descent_and_iteration_limit(self):
        _, diagnostics = generate_ellipsoid_points(24, (1, 1.5, 1), iterations=1,
                                                   return_diagnostics=True)
        self.assertEqual(diagnostics["method"], "riesz3")
        self.assertEqual(diagnostics["force_magnitude_power"], 4)
        self.assertEqual(diagnostics["iterations"], 1)
        self.assertFalse(diagnostics["converged"])
        self.assertLess(diagnostics["final_normalized_mean_pair_energy"],
                        diagnostics["initial_normalized_mean_pair_energy"])

    def test_centroid_remains_available_as_explicit_comparison(self):
        points, diagnostics = generate_ellipsoid_points(32, (1, 1.5, 1), method="centroid",
                                                       iterations=3, return_diagnostics=True)
        self.assertEqual(diagnostics["method"], "centroid")
        self.assertEqual(points.shape, (32, 3))
        np.testing.assert_allclose(np.sum((points/(1, 1.5, 1))**2, axis=1), 1, atol=2e-14)

    def test_area_estimator_on_symmetric_sphere_cells(self):
        center = np.asarray((2.0, -1.0, 0.5))
        points = 2 * np.concatenate((np.eye(3), -np.eye(3))) + center
        areas, errors, report = estimate_marker_areas(points, (2, 2, 2), center)
        expected = 4 * np.pi * 2**2
        self.assertAlmostEqual(report["total_surface_area"], expected)
        self.assertAlmostEqual(float(areas.sum()), expected)
        np.testing.assert_allclose(areas, expected / 6, rtol=0.01)
        self.assertTrue(np.all(errors >= 0))
        self.assertAlmostEqual(report["relative_rms_deviation"], report["area_cv"])
        scaled, _, _ = estimate_marker_areas(2 * (points-center), (4, 4, 4))
        np.testing.assert_allclose(scaled, 4 * areas, rtol=1e-14)

    def test_relaxation_improves_independently_measured_area_uniformity(self):
        _, _, relaxed = estimate_marker_areas(self.points, (1, 1, 3), samples_per_point=256)
        initial = generate_ellipsoid_points(512, (1, 1, 3), iterations=0)
        _, _, unrelaxed = estimate_marker_areas(initial, (1, 1, 3), samples_per_point=256)
        self.assertLess(relaxed["area_cv"], 0.7 * unrelaxed["area_cv"])
        self.assertLess(relaxed["area_cv"], 0.04)
        expected = 4 * np.pi * quad(lambda z: np.sqrt(9-8*z*z), 0, 1)[0]
        self.assertAlmostEqual(relaxed["total_surface_area"], expected, places=11)

    def test_area_estimation_rejects_off_surface_and_duplicate_points(self):
        for points in (np.zeros((1, 3)), np.asarray(((1, 0, 0), (1, 0, 0)))):
            with self.subTest(points=points), self.assertRaises(ValueError):
                estimate_marker_areas(points, (1, 1, 1))

    def test_single_marker_and_invalid_geometry(self):
        point = generate_ellipsoid_points(1, (1, 2, 3))
        self.assertEqual(point.shape, (1, 3))
        self.assertAlmostEqual(float(np.sum((point / (1, 2, 3))**2)), 1.0)
        for arguments in ({"num_points": 0, "axes": (1, 2, 3)},
                          {"num_points": 3, "axes": (1, 0, 3)},
                          {"num_points": 3, "axes": (1, np.nan, 3)}):
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                generate_ellipsoid_points(**arguments)


if __name__ == "__main__":
    unittest.main()
