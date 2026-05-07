# Pure-Python fallback for the Rust `_core` extension. Used automatically
# at import time if the compiled extension is missing.
#
# Ported from Java/src/mumfordShah/{GaussElim, GaussL2Mum}.java
# (Hohm, Storath & Weinmann) by Claude Sonnet coding agent, Anthropic, 2026.
"""Pure-Python implementations of the Phase 1 Rust primitives.

Only used as a fallback when ``mumfordshah2d._core`` is not importable.
The real Rust implementation in ``src/gauss_elim.rs`` is the
authoritative version; this file mirrors it line-for-line so that tests
pass identically with either backend.
"""

from __future__ import annotations

import numpy as np


def _build_lu(n: int, alpha: float) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Pre-compute the LU factors used by ``gauss_l2_mum_solve``.

    Mirrors ``GaussL2Mum.initialize`` in the Java source. ``diagonals[k]``
    has length ``k+1`` and ``factors[k]`` has length ``k`` (empty for
    ``k == 0``).
    """
    a = alpha + 1.0
    b = -alpha
    c = 2.0 * alpha + 1.0

    diagonals: list[np.ndarray] = [np.array([1.0])]
    factors: list[np.ndarray] = [np.zeros(0)]

    if n == 1:
        return diagonals, factors

    diagonals.append(np.array([a, a - b * b / a]))
    factors.append(np.array([-b / a]))

    for idx in range(2, n):
        d = np.zeros(idx + 1)
        cf = np.zeros(idx)
        d[0] = a
        for cur in range(1, idx):
            multi = -b / d[cur - 1]
            cf[cur - 1] = multi
            d[cur] = c + b * multi
        multi = -b / d[idx - 1]
        cf[idx - 1] = multi
        d[idx] = a + b * multi
        diagonals.append(d)
        factors.append(cf)

    return diagonals, factors


def _compute_mu_with_lu(
    y: np.ndarray, alpha: float, diagonals: list[np.ndarray], factors: list[np.ndarray]
) -> np.ndarray:
    """Apply forward + back substitution given pre-computed LU."""
    length = y.shape[0]
    if length == 0:
        return np.zeros(0)
    if length == 1:
        return y.astype(np.float64).copy()
    diag = diagonals[length - 1]
    factor = factors[length - 1]

    bvec = np.empty(length)
    bvec[0] = y[0]
    for i in range(1, length):
        bvec[i] = y[i] + factor[i - 1] * bvec[i - 1]

    result = np.empty(length)
    result[length - 1] = bvec[length - 1] / diag[length - 1]
    # Java uses `this.b == -alpha`; result[i] = (b[i] - this.b * result[i+1]) / diagonal[i]
    # → result[i] = (b[i] + alpha * result[i+1]) / diagonal[i]
    for i in range(length - 2, -1, -1):
        result[i] = (bvec[i] + alpha * result[i + 1]) / diag[i]
    return result


def gauss_l2_mum_solve(y: np.ndarray, alpha: float) -> np.ndarray:
    """Solve the L2-MS within-segment problem (pure-Python fallback)."""
    y = np.ascontiguousarray(np.asarray(y, dtype=np.float64))
    n = max(int(y.shape[0]), 1)
    diagonals, factors = _build_lu(n, float(alpha))
    return _compute_mu_with_lu(y, float(alpha), diagonals, factors)


def gauss_l2_mum_cost(y: np.ndarray, alpha: float) -> float:
    """L2-MS within-segment cost at the optimum (pure-Python fallback)."""
    y = np.ascontiguousarray(np.asarray(y, dtype=np.float64))
    if y.shape[0] == 0:
        return 0.0
    mu = gauss_l2_mum_solve(y, alpha)
    a = float(alpha)
    cost = 0.0
    for i in range(mu.shape[0] - 1):
        dy = y[i] - mu[i]
        dm = mu[i + 1] - mu[i]
        cost += dy * dy + a * dm * dm
    last = mu.shape[0] - 1
    dy = y[last] - mu[last]
    cost += dy * dy
    return float(cost)
