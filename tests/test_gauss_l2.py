# Tests for the L2-Mumford-Shah within-segment tridiagonal solver
# (`gauss_l2_mum_solve` / `gauss_l2_mum_cost`). Cross-checks the cached
# LU factorisation against `numpy.linalg.solve` on the dense system.
#
# Authored 2026 by Claude Sonnet coding agent (Anthropic).

import numpy as np
import pytest

from mumfordshah2d import gauss_l2_mum_cost, gauss_l2_mum_solve


def dense_l2_ms_matrix(n: int, alpha: float) -> np.ndarray:
    """Build the dense L2-MS tridiagonal matrix for a segment of length n.

    The matrix that ``gauss_l2_mum_solve`` factors implicitly:

        [1+α   -α                       ]
        [-α    1+2α  -α                 ]
        [      -α    1+2α  -α           ]
        [            ...                ]
        [                  -α    1+α    ]
    """
    if n == 0:
        return np.zeros((0, 0))
    if n == 1:
        return np.array([[1.0]])
    diag_interior = 1.0 + 2.0 * alpha
    diag_boundary = 1.0 + alpha
    off = -alpha
    m = np.zeros((n, n))
    for i in range(n):
        m[i, i] = diag_boundary if (i == 0 or i == n - 1) else diag_interior
        if i > 0:
            m[i, i - 1] = off
        if i < n - 1:
            m[i, i + 1] = off
    return m


@pytest.mark.parametrize("alpha", [0.1, 0.5, 1.0, 2.5, 10.0])
@pytest.mark.parametrize("n", [1, 2, 3, 5, 8, 16])
def test_matches_dense_solve(alpha: float, n: int):
    rng = np.random.default_rng(0xC0FFEE + n + int(alpha * 1000))
    y = rng.standard_normal(n)
    mu = gauss_l2_mum_solve(y, alpha)
    if n == 0:
        assert mu.shape == (0,)
        return
    if n == 1:
        np.testing.assert_array_equal(mu, y)
        return
    want = np.linalg.solve(dense_l2_ms_matrix(n, alpha), y)
    np.testing.assert_allclose(mu, want, atol=1e-10, rtol=1e-10)


def test_alpha_zero_is_identity():
    # With alpha=0 the system is the identity and mu == y.
    for n in [1, 4, 7, 32]:
        y = np.linspace(-2.0, 3.0, n)
        np.testing.assert_allclose(gauss_l2_mum_solve(y, 0.0), y, atol=1e-12)


def test_constant_input_is_fixed_point():
    # If y_i = c for all i, the optimum is mu = y for every alpha.
    for alpha in [0.0, 0.5, 5.0, 100.0]:
        y = np.full(20, 2.5)
        mu = gauss_l2_mum_solve(y, alpha)
        np.testing.assert_allclose(mu, y, atol=1e-10)


def test_large_alpha_drives_toward_constant():
    # As alpha → ∞, the solution tends to the global mean.
    y = np.linspace(-1.0, 1.0, 40)
    mu = gauss_l2_mum_solve(y, 1e8)
    mean = float(y.mean())
    spread = float(mu.max() - mu.min())
    np.testing.assert_allclose(mu.mean(), mean, atol=1e-6)
    assert spread < 1e-3, f"Expected near-constant solution, got spread {spread}"


def test_cost_zero_for_constant_input():
    # Constant input → mu = y → data error 0, smoothness 0 → total 0.
    y = np.full(10, 1.5)
    for alpha in [0.0, 0.5, 5.0]:
        cost = gauss_l2_mum_cost(y, alpha)
        assert abs(cost) < 1e-12, f"alpha={alpha}: cost={cost}"


def test_cost_alpha_zero_equals_zero():
    # With alpha=0 the optimum is mu = y, so the L2 data error is 0 too.
    y = np.array([1.0, -2.0, 3.0, 0.5, -0.1])
    cost = gauss_l2_mum_cost(y, 0.0)
    assert abs(cost) < 1e-12


def test_cost_matches_manual_evaluation():
    # Verify that gauss_l2_mum_cost returns the L2-MS energy at the
    # solver's mu, by computing the same energy by hand.
    y = np.array([0.0, 1.0, 1.0, 0.0, 2.0])
    alpha = 0.7
    mu = gauss_l2_mum_solve(y, alpha)
    manual = float(np.sum((y - mu) ** 2)) + alpha * float(
        np.sum(np.diff(mu) ** 2)
    )
    cost = gauss_l2_mum_cost(y, alpha)
    np.testing.assert_allclose(cost, manual, atol=1e-12)


def test_cost_is_positive_for_noisy_input():
    rng = np.random.default_rng(7)
    y = rng.standard_normal(15)
    cost = gauss_l2_mum_cost(y, 1.0)
    assert cost > 0.0
