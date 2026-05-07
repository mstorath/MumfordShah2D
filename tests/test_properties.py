# Property-based tests for Phase 1 primitives, using Hypothesis. These
# complement the closed-form tests in `test_utils.py` / `test_gauss_l2.py`
# / `test_prox.py` by exercising the algorithms over distributions of
# random inputs and asserting invariants that must hold mathematically.
#
# An invariant violated here often points to a real port-vs-original
# mismatch: invariants are insensitive to scaling, translation, etc., so
# any drift between a literal-translation bug and the original semantics
# tends to show up.
#
# Authored 2026 by Claude Sonnet coding agent (Anthropic).

import numpy as np
import pytest
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from hypothesis.extra import numpy as hnp

from mumfordshah2d import (
    expand_weights,
    gauss_l2_mum_cost,
    gauss_l2_mum_solve,
    hard_threshold,
    make_prox_inpaint,
    make_prox_l1w,
    make_prox_l2w,
    psnr,
    rotate90,
    soft_threshold,
)


HSETTINGS = settings(
    max_examples=120,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture],
)


finite_floats = st.floats(
    allow_nan=False,
    allow_infinity=False,
    min_value=-1e3,
    max_value=1e3,
    width=64,
)
positive_floats = st.floats(
    allow_nan=False,
    allow_infinity=False,
    min_value=1e-6,
    max_value=1e3,
    width=64,
)


def _array_1d(min_n: int = 1, max_n: int = 64):
    return hnp.arrays(
        dtype=np.float64,
        shape=hnp.array_shapes(min_dims=1, max_dims=1, min_side=min_n, max_side=max_n),
        elements=finite_floats,
    )


# ---------------------------------------------------------------------------
# soft_threshold properties
# ---------------------------------------------------------------------------


@HSETTINGS
@given(
    x=_array_1d(),
    tau1=positive_floats,
    tau2=positive_floats,
)
def test_soft_threshold_additive_composition(x, tau1, tau2):
    """Soft-thresholding twice (with arbitrary positive thresholds) equals
    soft-thresholding once with the *sum* of the thresholds:
    ``soft(soft(x, tau1), tau2) == soft(x, tau1 + tau2)``.

    This is the proximal-operator composition law for the L1 norm and is
    far more sensitive to off-by-one bugs than mere idempotency (which
    soft does NOT satisfy at the same threshold).
    """
    twice = soft_threshold(soft_threshold(x, tau1), tau2)
    once = soft_threshold(x, tau1 + tau2)
    np.testing.assert_allclose(twice, once, atol=1e-10, rtol=1e-10)


@HSETTINGS
@given(x=_array_1d(), tau=positive_floats)
def test_soft_threshold_one_step_pulls_toward_zero(x, tau):
    """A single soft-threshold step never pushes any entry past zero
    (the pull is one-sided): sign(soft(x)) is sign(x) or 0."""
    out = soft_threshold(x, tau)
    pulled_too_far = (np.sign(out) != 0) & (np.sign(out) != np.sign(x))
    assert not pulled_too_far.any()


@HSETTINGS
@given(x=_array_1d(), tau=positive_floats)
def test_soft_threshold_shrinks(x, tau):
    """soft_threshold never increases magnitude: |soft(x, tau)| <= |x|."""
    out = soft_threshold(x, tau)
    assert np.all(np.abs(out) <= np.abs(x) + 1e-12)


@HSETTINGS
@given(x=_array_1d(), tau=positive_floats)
def test_soft_threshold_preserves_sign_or_zeros(x, tau):
    """Soft-thresholded entries either have the same sign as x or are
    exactly zero (within FP)."""
    out = soft_threshold(x, tau)
    same_sign_or_zero = (np.sign(out) == np.sign(x)) | (out == 0.0)
    assert np.all(same_sign_or_zero)


@HSETTINGS
@given(x=_array_1d(), tau=positive_floats, scale=st.floats(0.5, 4.0, allow_nan=False))
def test_soft_threshold_homogeneity(x, tau, scale):
    """soft(s x, s tau) = s soft(x, tau) for s > 0 (positive homogeneity)."""
    np.testing.assert_allclose(
        soft_threshold(scale * x, scale * tau),
        scale * soft_threshold(x, tau),
        atol=1e-10,
        rtol=1e-10,
    )


# ---------------------------------------------------------------------------
# hard_threshold properties
# ---------------------------------------------------------------------------


@HSETTINGS
@given(x=_array_1d(), tau=positive_floats)
def test_hard_threshold_is_passthrough_or_zero(x, tau):
    out = hard_threshold(x, tau)
    is_pass = out == x
    is_zero = out == 0.0
    assert np.all(is_pass | is_zero)


@HSETTINGS
@given(x=_array_1d(), tau=positive_floats)
def test_hard_threshold_idempotent(x, tau):
    once = hard_threshold(x, tau)
    twice = hard_threshold(once, tau)
    np.testing.assert_array_equal(twice, once)


# ---------------------------------------------------------------------------
# expand_weights properties
# ---------------------------------------------------------------------------


@HSETTINGS
@given(
    rows=st.integers(min_value=1, max_value=8),
    cols=st.integers(min_value=1, max_value=8),
    chans=st.integers(min_value=2, max_value=4),
)
def test_expand_weights_replicates_to_channel_count(rows, cols, chans):
    """When the channel counts differ, weights are replicated along the
    new trailing axis."""
    f = np.zeros((rows, cols, chans))
    w = np.ones((rows, cols))
    out = expand_weights(w, f)
    assert out.shape == (rows, cols, chans)


@HSETTINGS
@given(
    rows=st.integers(min_value=1, max_value=8),
    cols=st.integers(min_value=1, max_value=8),
)
def test_expand_weights_no_op_when_channel_counts_match(rows, cols):
    """When the channel counts already agree, weights pass through
    unchanged — even if both arrays are 2-D (1 channel by MATLAB
    semantics)."""
    f = np.zeros((rows, cols))
    w = np.ones((rows, cols))
    out = expand_weights(w, f)
    assert out.shape == w.shape
    np.testing.assert_array_equal(out, w)


# ---------------------------------------------------------------------------
# rotate90 properties
# ---------------------------------------------------------------------------


@HSETTINGS
@given(
    arr=hnp.arrays(
        dtype=np.float64,
        shape=hnp.array_shapes(min_dims=3, max_dims=3, min_side=2, max_side=8),
        elements=finite_floats,
    ),
)
def test_rotate90_period_four(arr):
    """rotate90(arr, 4) == arr."""
    rotated = arr.copy()
    for _ in range(4):
        rotated = rotate90(rotated, 1)
    np.testing.assert_array_equal(rotated, arr)


@HSETTINGS
@given(
    arr=hnp.arrays(
        dtype=np.float64,
        shape=hnp.array_shapes(min_dims=3, max_dims=3, min_side=2, max_side=8),
        elements=finite_floats,
    ),
    k=st.integers(min_value=-4, max_value=4),
)
def test_rotate90_inverse(arr, k):
    np.testing.assert_array_equal(
        rotate90(rotate90(arr, k), -k),
        arr,
    )


# ---------------------------------------------------------------------------
# psnr properties
# ---------------------------------------------------------------------------


@HSETTINGS
@given(arr=hnp.arrays(dtype=np.float64, shape=10, elements=positive_floats))
def test_psnr_self_is_inf(arr):
    assert np.isinf(psnr(arr, arr))


# ---------------------------------------------------------------------------
# gauss_l2_mum_solve properties — these are the algorithm's defining laws
# ---------------------------------------------------------------------------


@HSETTINGS
@given(
    y=_array_1d(min_n=1, max_n=32),
    alpha=st.floats(min_value=0.0, max_value=1e3, allow_nan=False),
    c=finite_floats,
)
def test_gauss_l2_translation_equivariance(y, alpha, c):
    """Solver is translation-equivariant: solve(y + c) = solve(y) + c.

    Adds a constant to both data term and the smoothness target — the
    objective shifts by a constant, so the optimum shifts by the same
    constant. Highly sensitive to off-by-one bugs in either index.
    """
    mu = gauss_l2_mum_solve(y, alpha)
    mu_shifted = gauss_l2_mum_solve(y + c, alpha)
    np.testing.assert_allclose(mu_shifted, mu + c, atol=1e-9, rtol=1e-9)


@HSETTINGS
@given(
    y=_array_1d(min_n=1, max_n=32),
    alpha=st.floats(min_value=0.0, max_value=1e3, allow_nan=False),
    s=st.floats(min_value=-50.0, max_value=50.0, allow_nan=False),
)
def test_gauss_l2_scale_equivariance(y, alpha, s):
    """Solver is scale-equivariant: solve(s y) = s solve(y) for any s.

    Both data and smoothness terms are quadratic in mu, so multiplying y
    by s rescales the optimum by s.
    """
    assume(abs(s) > 1e-6)
    mu = gauss_l2_mum_solve(y, alpha)
    mu_scaled = gauss_l2_mum_solve(s * y, alpha)
    np.testing.assert_allclose(mu_scaled, s * mu, atol=1e-9, rtol=1e-9)


@HSETTINGS
@given(y=_array_1d(min_n=2, max_n=16))
def test_gauss_l2_alpha_zero_identity(y):
    """At alpha=0 there is no smoothness coupling and mu = y exactly."""
    np.testing.assert_allclose(gauss_l2_mum_solve(y, 0.0), y, atol=1e-12)


@HSETTINGS
@given(y=_array_1d(min_n=2, max_n=20))
def test_gauss_l2_solution_is_average_at_large_alpha(y):
    """As alpha → ∞ the optimum is the global mean."""
    mu = gauss_l2_mum_solve(y, 1e8)
    np.testing.assert_allclose(mu.mean(), y.mean(), atol=1e-6)
    assert (mu.max() - mu.min()) < 1e-3


@HSETTINGS
@given(
    y=_array_1d(min_n=2, max_n=16),
    alpha=st.floats(min_value=1e-3, max_value=10.0, allow_nan=False),
)
def test_gauss_l2_cost_matches_manual(y, alpha):
    """gauss_l2_mum_cost == sum (y-mu)^2 + alpha * sum diff(mu)^2 at mu*."""
    mu = gauss_l2_mum_solve(y, alpha)
    manual = float(np.sum((y - mu) ** 2)) + alpha * float(np.sum(np.diff(mu) ** 2))
    np.testing.assert_allclose(gauss_l2_mum_cost(y, alpha), manual, atol=1e-10, rtol=1e-10)


@HSETTINGS
@given(
    y=_array_1d(min_n=2, max_n=12),
    alpha=st.floats(min_value=1e-3, max_value=10.0, allow_nan=False),
)
def test_gauss_l2_minimum_via_perturbation(y, alpha):
    """The cost at the solver's mu must be lower than (or equal to) the
    cost at mu plus a small random perturbation. This is an empirical
    optimality check — does not require an analytic gradient."""
    mu = gauss_l2_mum_solve(y, alpha)
    cost = gauss_l2_mum_cost(y, alpha)

    # Manually evaluate the same cost at mu + eps perturbation.
    rng = np.random.default_rng(int(np.sum(np.abs(y)) * 1e6) % 2**31)
    eps = rng.standard_normal(mu.shape) * 1e-3
    mu_p = mu + eps
    cost_p = float(np.sum((y - mu_p) ** 2)) + alpha * float(np.sum(np.diff(mu_p) ** 2))
    # cost_p must be >= cost (within FP tolerance).
    assert cost_p >= cost - 1e-9


# ---------------------------------------------------------------------------
# Prox properties: agreement with their analytic definitions
# ---------------------------------------------------------------------------


@st.composite
def _same_shape_pair(draw, min_n=1, max_n=8):
    """Generate a pair of 1-D arrays guaranteed to share their shape — so
    Hypothesis doesn't waste cycles filtering out mismatched draws."""
    n = draw(st.integers(min_value=min_n, max_value=max_n))
    a = draw(hnp.arrays(dtype=np.float64, shape=n, elements=finite_floats))
    b = draw(hnp.arrays(dtype=np.float64, shape=n, elements=finite_floats))
    return a, b


@HSETTINGS
@given(pair=_same_shape_pair(), tau=positive_floats, w=positive_floats)
def test_prox_l2w_closed_form(pair, tau, w):
    """make_prox_l2w(f, w)(z, tau) == (w f + tau z) / (w + tau)."""
    f, z = pair
    weights = np.full_like(f, w)
    out = make_prox_l2w(f, weights)(z, tau)
    expected = (weights * f + tau * z) / (weights + tau)
    np.testing.assert_allclose(out, expected, atol=1e-10, rtol=1e-10)


@HSETTINGS
@given(pair=_same_shape_pair(), tau=positive_floats, w=positive_floats)
def test_prox_l1w_closed_form(pair, tau, w):
    """make_prox_l1w(f, w)(z, tau) == soft(z - f, w/tau) + f."""
    f, z = pair
    weights = np.full_like(f, w)
    out = make_prox_l1w(f, weights)(z, tau)
    expected = soft_threshold(z - f, weights / tau) + f
    np.testing.assert_allclose(out, expected, atol=1e-10, rtol=1e-10)


@HSETTINGS
@given(pair=_same_shape_pair(), mask_seed=st.integers(min_value=0, max_value=2**31 - 1))
def test_prox_inpaint_is_pixelwise_select(pair, mask_seed):
    f, z = pair
    mask = (
        np.random.default_rng(mask_seed)
        .integers(0, 2, size=f.shape)
        .astype(np.float64)
    )
    out = make_prox_inpaint(f, mask)(z, tau=1.0)
    expected = np.where(mask.astype(bool), f, z)
    np.testing.assert_array_equal(out, expected)
