# Changelog

## 0.1.0 — 2026-05-06 (Phase 1 — scaffolding + small primitives)

First alpha of the Python/Rust port. **Not feature-complete** — the algorithmic
2D Mumford-Shah solver is still in progress. See `PORTED_BY.md` for the phasing
plan.

### What's included

- **Build system.** `Cargo.toml` (PyO3 0.22 + ndarray + rayon) and
  `pyproject.toml` (maturin backend, Python ≥3.9, numpy ≥1.20, scipy ≥1.7).
- **Rust core skeleton.** `src/lib.rs` with PyO3 bindings, `src/image.rs`
  (`Point` + `Image`), `src/gauss_elim.rs` (`GaussElim` trait + `GaussL2Mum`
  impl).
- **Python API skeleton.** `mumfordshah2d/__init__.py` with versioned exports,
  `utils.py` (5 utility ports), `prox.py` (4 prox handles), `_core_py.py`
  (pure-Python fallbacks for primitives).
- **Tests.** `tests/conftest.py` with fixtures, `tests/test_utils.py`,
  `tests/test_gauss_l2.py` — covering thresholds, weight expansion, image
  array round-trips, and the tridiagonal solve against `numpy.linalg.solve`.

### What's NOT included yet

- `TautString` (Phase 2)
- 1D Mumford-Shah DPs (Phase 3)
- 2D direction processor + outer ADMM (Phase 4)
- Inverse problems (deconv / inpaint / Radon) and the four MATLAB demos as
  Python equivalents (Phase 5)

### Verification

- `pytest tests/` green (≥10 tests).
- `pip install -e .` succeeds via maturin.
- No edits to existing `.m` or `.java` files.
