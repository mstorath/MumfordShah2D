# Tests for mumfordshah2d.utils — small array helpers ported from
# MATLAB Auxiliary/. Authored 2026 by Claude Sonnet coding agent (Anthropic).

import numpy as np
import pytest

from mumfordshah2d import (
    expand_weights,
    hard_threshold,
    psnr,
    rotate90,
    soft_threshold,
)


class TestSoftThreshold:
    def test_below_threshold_is_zero(self):
        x = np.array([0.1, -0.2, 0.3, -0.4])
        out = soft_threshold(x, 0.5)
        np.testing.assert_array_equal(out, np.zeros(4))

    def test_above_threshold_shrinks_by_tau(self):
        x = np.array([1.0, -2.0, 3.0])
        out = soft_threshold(x, 0.5)
        np.testing.assert_allclose(out, [0.5, -1.5, 2.5])

    def test_per_element_tau(self):
        x = np.array([1.0, 1.0, 1.0])
        tau = np.array([0.1, 0.5, 0.9])
        out = soft_threshold(x, tau)
        np.testing.assert_allclose(out, [0.9, 0.5, 0.1])

    def test_zero_input(self):
        np.testing.assert_array_equal(soft_threshold(np.zeros(5), 0.3), np.zeros(5))


class TestHardThreshold:
    def test_passes_above_threshold(self):
        x = np.array([0.1, 0.4, 1.0, -2.0, 0.5])
        out = hard_threshold(x, 0.5)
        # |x| >= 0.5 → keep; otherwise zero. 0.5 itself counts as above.
        np.testing.assert_array_equal(out, [0.0, 0.0, 1.0, -2.0, 0.5])

    def test_zero_below(self):
        x = np.array([0.1, -0.2, 0.3])
        np.testing.assert_array_equal(hard_threshold(x, 1.0), np.zeros(3))


class TestExpandWeights:
    def test_2d_weights_to_3d_image(self):
        w = np.ones((4, 5))
        f = np.zeros((4, 5, 3))
        out = expand_weights(w, f)
        assert out.shape == (4, 5, 3)
        np.testing.assert_array_equal(out, np.ones((4, 5, 3)))

    def test_already_matching_passes_through(self):
        w = np.full((4, 5, 3), 0.7)
        f = np.zeros((4, 5, 3))
        out = expand_weights(w, f)
        np.testing.assert_array_equal(out, w)

    def test_grayscale_image_returns_weights_unchanged(self):
        # Both inputs have ndim < 3 (so MATLAB's size(., 3) == 1 for both),
        # the channel counts already match, and expand_weights is a no-op.
        w = np.full((4, 5), 0.7)
        f = np.zeros((4, 5))
        out = expand_weights(w, f)
        assert out.shape == (4, 5)
        np.testing.assert_array_equal(out, w)

    def test_1d_inputs_match_matlab(self):
        # MATLAB happily takes vectors here; size(v, 3) is 1 either way and
        # the function is a no-op.
        w = np.array([1.0, 1.0, 1.0])
        f = np.array([0.0, 0.0, 0.0])
        out = expand_weights(w, f)
        np.testing.assert_array_equal(out, w)


class TestRotate90:
    def test_round_trip_4_rotations(self):
        f = np.arange(2 * 3 * 4, dtype=np.float64).reshape(2, 3, 4)
        rotated = f
        for _ in range(4):
            rotated = rotate90(rotated, 1)
        np.testing.assert_array_equal(rotated, f)

    def test_minus_one_inverts_one(self):
        f = np.arange(2 * 3 * 4, dtype=np.float64).reshape(2, 3, 4)
        np.testing.assert_array_equal(rotate90(rotate90(f, 1), -1), f)

    def test_requires_3d(self):
        with pytest.raises(ValueError):
            rotate90(np.zeros((4, 4)), 1)


class TestPsnr:
    def test_perfect_reconstruction_is_inf(self):
        gt = np.ones(10)
        assert np.isinf(psnr(gt, gt))

    def test_known_value(self):
        # gt = ones of length 10, signal differs by 0.1 at every entry.
        # mse = 0.01, max(gt) = 1, psnr = 10 log10(1 / 0.01) = 20 dB.
        gt = np.ones(10)
        sig = gt + 0.1
        np.testing.assert_allclose(psnr(gt, sig), 20.0, atol=1e-9)
