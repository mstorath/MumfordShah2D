# Changelog

## [Unreleased]

## 0.2.0 — 2026-05-07 (Phase 4 — 4-connected 2-D ADMM)

The 2-D L2 Mumford-Shah solver lands. The package transitions from "alpha,
1-D primitives only" to "beta, the full 2-D denoising pipeline that the
MATLAB driver exposes (subset of cases)."

### Added

- **Rust core.**
  - `src/mumfordshah_1d.rs` — 1-D L2 Mumford-Shah solver (DP partition
    search + within-segment via `GaussL2Mum`). Direct port of
    `Java/src/mumfordShah/{AbstractPotts1D,L2MumfordShah1DFaster}`.
  - `src/direction_processor.rs` — drives the 1-D solver over every stripe
    of a given direction. Compute parallelised across stripes via rayon.
  - `src/admm.rs` — 4-connected (anisotropic) 2-D ADMM driver
    (`admm_4connected_l2_ms`) with the `Prox` trait and `L2DataProx`
    default (matches MATLAB's `makeProxL2w` at uniform weights). Default
    μ schedule `μ(k) = k².⁰¹ · 10⁻⁶` matches `mumfordShah2D.m`.
  - `src/image.rs` — added `rotate90_cw/ccw`, `n_stripes`, and named
    direction constants (`directions::{ROW_STEP, NW_SE, KNIGHT_2_1,
    KNIGHT_1_2}`).

- **Python facade.**
  - `mumfordshah2d.min_l2_mum_2d(f, gamma, alpha, …)` accepts 2-D
    (grayscale) or 3-D (multi-channel) input, validates parameters,
    routes to the Rust core, and returns the result with matching
    dimensionality.

- **Tests.**
  - 37 Rust unit tests (`mumfordshah_1d`, `direction_processor`, `admm`,
    `image` extensions).
  - 4 live MATLAB-parity tests at `atol=1e-9` for the 4-connected ADMM
    on 4×4 inputs, max_iter ∈ {20, 100} (`tests/test_matlab_parity_2d.py`).

### Changed

- `Cargo.toml`: `crate-type` extended to `["cdylib", "rlib"]` and
  `extension-module` is now an optional default-on feature, so `cargo test
  --no-default-features` can link unit tests against libpython directly.
  The shipped wheel is unchanged.

### Not yet supported (planned for 0.3.0+)

- 8-connected (near-isotropic) and knight-move neighbourhoods (S=4 / S=8).
- The ρ-coupled MATLAB-default behaviour (rhoFlag=true). The current Rust
  port matches MATLAB only when `nuSeq = @(k) 0` is passed; the default
  MATLAB call uses inter-direction ρ dual variables that the Rust port
  doesn't yet implement.
- Custom μ schedules and arbitrary prox handles via the Python API.
  (The Rust core already accepts `mu_seq: impl Fn(usize) -> f64` and
  any `impl Prox`; only the PyO3 wrapper is constrained to the defaults.)

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
