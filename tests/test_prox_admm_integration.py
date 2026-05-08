# Integration tests for the user-supplied prox path through min_l2_mum_2d.
#
# These exercise the PyO3 callback bridge (PyProx in src/lib.rs): a Python
# callable wrapped as a Rust `Prox` impl. Bit-identical equivalence with
# the Rust-native L2DataProx path is the strongest regression check —
# anything that breaks the bridge will surface here as nonzero diffs.

from __future__ import annotations

import numpy as np
import pytest

from mumfordshah2d import min_l2_mum_2d
from mumfordshah2d.prox import make_prox_inpaint, make_prox_l2w


@pytest.fixture
def image_2d():
    rng = np.random.default_rng(0)
    return rng.standard_normal((20, 20))


@pytest.fixture
def image_3d():
    rng = np.random.default_rng(1)
    return rng.standard_normal((15, 15, 3))


def test_user_l2_prox_with_uniform_weights_matches_default(image_2d):
    """prox=make_prox_l2w(f, ones) MUST be bit-identical to no prox at all,
    because both reduce to the same closed form (f + λ·z) / (1 + λ).
    Any drift here is a bug in the PyProx bridge or the shape-shim layer."""
    out_default = min_l2_mum_2d(image_2d, gamma=2.0, alpha=1.0, max_iter=200)
    prox = make_prox_l2w(image_2d, np.ones_like(image_2d))
    out_user = min_l2_mum_2d(image_2d, gamma=2.0, alpha=1.0, max_iter=200, prox=prox)
    np.testing.assert_array_equal(out_default, out_user)


def test_user_l2_prox_with_weights_matches_weights_path(image_2d):
    """prox=make_prox_l2w(f, w) MUST match the equivalent weights= path
    (which constructs the same closed form Rust-side). The two paths
    compute the same prox; only the dispatch is different."""
    rng = np.random.default_rng(42)
    w = rng.uniform(0.5, 2.0, size=image_2d.shape)
    out_weights = min_l2_mum_2d(image_2d, gamma=2.0, alpha=1.0, max_iter=200, weights=w)
    prox = make_prox_l2w(image_2d, w)
    out_user = min_l2_mum_2d(image_2d, gamma=2.0, alpha=1.0, max_iter=200, prox=prox)
    np.testing.assert_allclose(out_weights, out_user, atol=1e-12, rtol=1e-12)


def test_user_l2_prox_multichannel(image_3d):
    """Shape shim must round-trip (rows, cols, channels) through the Rust
    core's (channels, rows, cols) layout cleanly."""
    prox = make_prox_l2w(image_3d, np.ones(image_3d.shape[:2]))
    out = min_l2_mum_2d(image_3d, gamma=2.0, alpha=1.0, max_iter=200, prox=prox)
    assert out.shape == image_3d.shape
    assert np.all(np.isfinite(out))


def test_inpaint_prox_fills_masked_region(image_2d):
    """make_prox_inpaint clamps known pixels to f and lets z drive missing
    ones. The output inside the mask must therefore differ from the input
    after the ADMM has had a chance to spread information across the hole."""
    mask = np.ones_like(image_2d)
    mask[5:15, 5:15] = 0  # 10×10 hole — much larger than any natural feature
    prox = make_prox_inpaint(image_2d, mask)
    out = min_l2_mum_2d(image_2d, gamma=1.0, alpha=1.0, max_iter=2000, tol=1e-5, prox=prox)
    # Inside the hole, the inpainted values should NOT match the (random) input.
    inside_diff = np.abs(out[5:15, 5:15] - image_2d[5:15, 5:15]).max()
    assert inside_diff > 0.1, f"inpainting did not spread into hole (max diff {inside_diff})"


def test_prox_and_weights_conflict_raises(image_2d):
    prox = make_prox_l2w(image_2d, np.ones_like(image_2d))
    with pytest.raises(ValueError, match="either"):
        min_l2_mum_2d(
            image_2d,
            gamma=2.0,
            alpha=1.0,
            prox=prox,
            weights=np.ones_like(image_2d),
        )


def test_non_callable_prox_raises(image_2d):
    with pytest.raises(TypeError, match="callable"):
        min_l2_mum_2d(image_2d, gamma=2.0, alpha=1.0, prox="not a callable")
