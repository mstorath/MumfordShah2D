# Rust core vs pure-Python fallback parity. The fallback in `_core_py.py`
# mirrors the Rust algorithm in `src/gauss_elim.rs` line-for-line; tests
# must produce bit-identical (or near-bit-identical, modulo FP rounding
# order) results regardless of which backend is loaded.
#
# This file imports both backends *explicitly* and compares — independent
# of which one `mumfordshah2d` happens to have selected at import time.
#
# Authored 2026 by Claude Sonnet coding agent (Anthropic).

import numpy as np
import pytest

import mumfordshah2d
from mumfordshah2d import _core_py


def _try_import_rust():
    try:
        from mumfordshah2d import _core as rust
        return rust
    except ImportError:
        return None


_rust = _try_import_rust()


# Skip this whole module if the Rust extension hasn't been built — the
# tests have nothing to compare to. (`pytest -v` will still show the
# skip reason.)
pytestmark = pytest.mark.skipif(
    _rust is None,
    reason="Rust extension not built; nothing to compare against the pure-Python fallback",
)


def test_active_backend_is_consistent():
    """`mumfordshah2d.gauss_l2_mum_solve` should be the Rust binding when
    the extension is available."""
    assert mumfordshah2d._USING_RUST_CORE is True


@pytest.mark.parametrize("alpha", [0.0, 0.05, 1.0, 10.0, 1000.0])
@pytest.mark.parametrize("n", [1, 2, 3, 5, 9, 16, 33, 100])
def test_solve_rust_matches_pure_python(alpha: float, n: int):
    """For a wide range of (n, alpha), Rust and pure-Python solve give
    bit-identical (atol 1e-13) outputs on random inputs.

    This is the single best protection against future divergence between
    the two implementations: any tweak to the algorithm that lands in only
    one of them will surface here immediately.
    """
    rng = np.random.default_rng(seed=alpha.__hash__() ^ n)
    y = rng.standard_normal(n)
    rust_out = _rust.gauss_l2_mum_solve(y, alpha)
    py_out = _core_py.gauss_l2_mum_solve(y, alpha)
    np.testing.assert_allclose(rust_out, py_out, atol=1e-13, rtol=1e-13)


@pytest.mark.parametrize("alpha", [0.0, 0.5, 5.0, 100.0])
@pytest.mark.parametrize("n", [1, 4, 16, 64])
def test_cost_rust_matches_pure_python(alpha: float, n: int):
    rng = np.random.default_rng(seed=42 + n + int(alpha * 13))
    y = rng.standard_normal(n)
    rust_cost = _rust.gauss_l2_mum_cost(y, alpha)
    py_cost = _core_py.gauss_l2_mum_cost(y, alpha)
    np.testing.assert_allclose(rust_cost, py_cost, atol=1e-12, rtol=1e-12)


def test_extreme_alpha_does_not_diverge():
    """At alpha = 1e10 the LU has nearly singular blocks; both backends
    should still agree to high precision (or both degrade together)."""
    y = np.array([0.0, 1.0, -1.0, 2.0, -2.0, 0.5])
    rust_out = _rust.gauss_l2_mum_solve(y, 1e10)
    py_out = _core_py.gauss_l2_mum_solve(y, 1e10)
    np.testing.assert_allclose(rust_out, py_out, atol=1e-9, rtol=1e-9)


def test_special_inputs_give_identical_results():
    """All-zero, all-equal, and single-element inputs are tested exactly."""
    cases = [
        np.zeros(7),
        np.full(10, 3.5),
        np.array([42.0]),
        np.array([0.0, 0.0]),
        np.array([1.0, 1.0, 1.0, 1.0]),
    ]
    for y in cases:
        for alpha in [0.0, 1.0, 100.0]:
            rust_out = _rust.gauss_l2_mum_solve(y, alpha)
            py_out = _core_py.gauss_l2_mum_solve(y, alpha)
            np.testing.assert_array_equal(
                rust_out,
                py_out,
                err_msg=f"y={y}, alpha={alpha}",
            )
