# Tests for the four prox handles in mumfordshah2d.prox. Each handle is a
# direct port of a MATLAB closed-form lambda in
# `MumfordShah2D/Auxiliary/makeProx*.m`. We test:
#   * the closed-form on hand-derived inputs,
#   * boundary behaviour as tau→0 and tau→∞,
#   * weight-handling (uniform, per-pixel, multichannel),
#   * the documented "fixed point at f" property.
#
# Authored 2026 by Claude Sonnet coding agent (Anthropic).

import numpy as np
import pytest

from mumfordshah2d import (
    make_prox_inpaint,
    make_prox_l0w,
    make_prox_l1w,
    make_prox_l2w,
)


# ---------------------------------------------------------------------------
# make_prox_l2w  —  MATLAB: (weights .* f + tau * v) ./ (weights + tau)
# ---------------------------------------------------------------------------


class TestProxL2w:
    def test_closed_form_on_hand_input(self):
        f = np.array([1.0, 2.0, 3.0, 4.0])
        w = np.array([1.0, 1.0, 1.0, 1.0])
        z = np.array([5.0, 5.0, 5.0, 5.0])
        prox = make_prox_l2w(f, w)
        # Expected: (1*f + 2*z) / (1+2) = (f + 2z)/3
        out = prox(z, tau=2.0)
        expected = (1 * f + 2 * z) / 3.0
        np.testing.assert_allclose(out, expected, atol=1e-12)

    def test_tau_zero_recovers_f(self):
        # tau → 0: the data fidelity dominates; result → f.
        f = np.array([1.0, 2.0, 3.0])
        w = np.ones(3)
        z = np.array([99.0, 99.0, 99.0])
        prox = make_prox_l2w(f, w)
        out = prox(z, tau=1e-12)
        np.testing.assert_allclose(out, f, atol=1e-9)

    def test_large_tau_recovers_z(self):
        # tau → ∞: the proximal pull dominates; result → z.
        f = np.array([1.0, 2.0, 3.0])
        w = np.ones(3)
        z = np.array([99.0, 99.0, 99.0])
        prox = make_prox_l2w(f, w)
        out = prox(z, tau=1e12)
        np.testing.assert_allclose(out, z, atol=1e-6)

    def test_zero_weights_drops_data_term(self):
        # weights = 0 means f is irrelevant; output should equal z.
        f = np.array([1.0, 2.0, 3.0])
        w = np.zeros(3)
        z = np.array([7.0, 8.0, 9.0])
        prox = make_prox_l2w(f, w)
        out = prox(z, tau=2.0)
        np.testing.assert_allclose(out, z, atol=1e-12)

    def test_per_pixel_weights(self):
        # Heterogeneous weights: each pixel has its own pull strength.
        f = np.array([1.0, 1.0, 1.0])
        w = np.array([0.0, 1.0, 1e6])
        z = np.array([0.0, 0.0, 0.0])
        prox = make_prox_l2w(f, w)
        out = prox(z, tau=1.0)
        # Pixel 0: w=0 → output = z = 0
        # Pixel 1: w=1 → (f + z) / 2 = 0.5
        # Pixel 2: w=1e6 → output ≈ f = 1
        np.testing.assert_allclose(out[0], 0.0, atol=1e-12)
        np.testing.assert_allclose(out[1], 0.5, atol=1e-12)
        np.testing.assert_allclose(out[2], 1.0, atol=1e-5)

    def test_multichannel_with_2d_weights(self):
        # 2-D weights broadcast across the channel axis.
        f = np.zeros((4, 5, 3))
        f[..., 0] = 1.0
        f[..., 1] = 2.0
        f[..., 2] = 3.0
        w = np.ones((4, 5))  # 2-D → broadcast to (4, 5, 3)
        z = np.full((4, 5, 3), 0.0)
        prox = make_prox_l2w(f, w)
        out = prox(z, tau=1.0)
        # (1*f + 1*z) / 2 = f/2
        np.testing.assert_allclose(out, f / 2.0, atol=1e-12)

    def test_init_kwarg_is_ignored(self):
        # The closed-form prox doesn't need a warm-start; passing init must
        # not change the result.
        f = np.array([1.0, 2.0, 3.0])
        w = np.ones(3)
        z = np.array([5.0, 6.0, 7.0])
        prox = make_prox_l2w(f, w)
        a = prox(z, tau=1.0)
        b = prox(z, tau=1.0, init=np.array([99.0, 99.0, 99.0]))
        np.testing.assert_array_equal(a, b)


# ---------------------------------------------------------------------------
# make_prox_l1w  —  MATLAB: softThreshold(z - f, w/tau) + f
# ---------------------------------------------------------------------------


class TestProxL1w:
    def test_z_equals_f_returns_f(self):
        # When z = f, z - f = 0; soft_threshold of 0 is 0; output = f.
        f = np.array([0.0, 1.0, -2.0, 3.0])
        prox = make_prox_l1w(f, np.ones(4))
        out = prox(f, tau=1.0)
        np.testing.assert_array_equal(out, f)

    def test_within_tube_clamps_to_f(self):
        # |z - f| <= w/tau → soft-threshold output is 0; final = f.
        f = np.array([0.0, 0.0, 0.0])
        w = np.array([1.0, 1.0, 1.0])  # w/tau = 0.5 with tau=2
        z = np.array([0.4, -0.3, 0.49])
        prox = make_prox_l1w(f, w)
        out = prox(z, tau=2.0)  # threshold = 0.5
        np.testing.assert_allclose(out, np.zeros(3), atol=1e-12)

    def test_outside_tube_shrinks_by_w_over_tau(self):
        f = np.array([0.0, 0.0, 0.0])
        w = np.array([1.0, 1.0, 1.0])
        z = np.array([1.0, -1.0, 2.5])
        prox = make_prox_l1w(f, w)
        out = prox(z, tau=2.0)  # threshold = 0.5
        # soft(1.0, 0.5) = 0.5; soft(-1.0, 0.5) = -0.5; soft(2.5, 0.5) = 2.0
        np.testing.assert_allclose(out, [0.5, -0.5, 2.0], atol=1e-12)

    def test_per_pixel_weights(self):
        # Different weights → different thresholds per pixel.
        f = np.zeros(3)
        w = np.array([0.0, 1.0, 100.0])
        z = np.array([0.7, 0.7, 0.7])
        prox = make_prox_l1w(f, w)
        out = prox(z, tau=1.0)  # thresholds = [0, 1, 100]
        # Pixel 0: thr=0, soft(0.7, 0) = 0.7
        # Pixel 1: thr=1, soft(0.7, 1) = 0
        # Pixel 2: thr=100, soft(0.7, 100) = 0
        np.testing.assert_allclose(out, [0.7, 0.0, 0.0], atol=1e-12)


# ---------------------------------------------------------------------------
# make_prox_l0w  —  MATLAB: hardThreshold(z - f, sqrt(2*w/tau)) + f
# ---------------------------------------------------------------------------


class TestProxL0w:
    def test_z_equals_f_returns_f(self):
        f = np.array([1.0, 2.0, 3.0])
        prox = make_prox_l0w(f, np.ones(3))
        out = prox(f, tau=1.0)
        np.testing.assert_array_equal(out, f)

    def test_within_threshold_clamps_to_f(self):
        # threshold = sqrt(2 * w / tau) = sqrt(2 * 1 / 2) = 1.0
        f = np.array([0.0, 0.0, 0.0])
        w = np.array([1.0, 1.0, 1.0])
        z = np.array([0.5, -0.99, 0.99])
        prox = make_prox_l0w(f, w)
        out = prox(z, tau=2.0)
        np.testing.assert_allclose(out, np.zeros(3), atol=1e-12)

    def test_above_threshold_passes_through_unchanged(self):
        f = np.array([0.0, 0.0])
        w = np.array([1.0, 1.0])
        z = np.array([1.5, -2.0])  # both > threshold = 1.0
        prox = make_prox_l0w(f, w)
        out = prox(z, tau=2.0)
        # hard_threshold passes through z-f when |z-f| >= threshold; output
        # = (z - f) + f = z.
        np.testing.assert_allclose(out, z, atol=1e-12)

    def test_threshold_scales_with_weight_sqrt(self):
        # threshold = sqrt(2 w / tau). Doubling w → threshold * sqrt(2).
        f = np.zeros(1)
        z = np.array([1.4])
        # tau=2 → threshold = sqrt(w). With w=1, threshold = 1. With w=2,
        # threshold = sqrt(2) ≈ 1.414.
        out_w1 = make_prox_l0w(f, np.array([1.0]))(z, tau=2.0)
        out_w2 = make_prox_l0w(f, np.array([2.0]))(z, tau=2.0)
        # |1.4| >= 1: passes; |1.4| < sqrt(2): clipped to 0.
        np.testing.assert_allclose(out_w1, [1.4], atol=1e-12)
        np.testing.assert_allclose(out_w2, [0.0], atol=1e-12)


# ---------------------------------------------------------------------------
# make_prox_inpaint  —  MATLAB: f .* mask + z .* ~mask
# ---------------------------------------------------------------------------


class TestProxInpaint:
    def test_full_mask_returns_f(self):
        f = np.array([1.0, 2.0, 3.0])
        mask = np.array([1.0, 1.0, 1.0])
        z = np.array([99.0, 99.0, 99.0])
        prox = make_prox_inpaint(f, mask)
        np.testing.assert_array_equal(prox(z, tau=1.0), f)

    def test_empty_mask_returns_z(self):
        f = np.array([1.0, 2.0, 3.0])
        mask = np.array([0.0, 0.0, 0.0])
        z = np.array([99.0, 99.0, 99.0])
        prox = make_prox_inpaint(f, mask)
        np.testing.assert_array_equal(prox(z, tau=1.0), z)

    def test_mixed_mask_per_pixel(self):
        f = np.array([1.0, 2.0, 3.0, 4.0])
        mask = np.array([1.0, 0.0, 1.0, 0.0])
        z = np.array([10.0, 20.0, 30.0, 40.0])
        prox = make_prox_inpaint(f, mask)
        out = prox(z, tau=1.0)
        np.testing.assert_array_equal(out, [1.0, 20.0, 3.0, 40.0])

    def test_tau_is_ignored(self):
        # Inpainting prox is a hard masked overwrite; tau has no effect.
        f = np.array([1.0, 2.0, 3.0])
        mask = np.array([1.0, 0.0, 1.0])
        z = np.array([10.0, 20.0, 30.0])
        prox = make_prox_inpaint(f, mask)
        a = prox(z, tau=1e-9)
        b = prox(z, tau=1e9)
        np.testing.assert_array_equal(a, b)
        np.testing.assert_array_equal(a, [1.0, 20.0, 3.0])

    def test_2d_mask_broadcasts_to_3d_image(self):
        f = np.full((4, 5, 3), 1.0)
        mask = np.zeros((4, 5))
        mask[0, 0] = 1.0  # one known pixel only
        z = np.full((4, 5, 3), 7.0)
        prox = make_prox_inpaint(f, mask)
        out = prox(z, tau=1.0)
        # Pixel (0,0): all channels come from f.
        np.testing.assert_array_equal(out[0, 0, :], [1.0, 1.0, 1.0])
        # All other pixels: come from z.
        assert np.all(out[1:, :, :] == 7.0)
        assert np.all(out[0, 1:, :] == 7.0)


# ---------------------------------------------------------------------------
# Cross-prox: shared invariants
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory",
    [make_prox_l2w, make_prox_l1w, make_prox_l0w],
)
def test_prox_at_f_returns_f(factory):
    """Sanity: when z = f and the prox is well-behaved, output = f."""
    f = np.array([1.0, -2.0, 3.5, 0.0])
    w = np.ones(4)
    prox = factory(f, w)
    out = prox(f, tau=1.0)
    np.testing.assert_allclose(out, f, atol=1e-12)
