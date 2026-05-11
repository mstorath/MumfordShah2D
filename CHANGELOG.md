# Changelog

## [Unreleased]

## 0.5.2 — 2026-05-11 (metadata maintenance release)

No algorithmic changes.

### Changed

- README rewritten and harmonised with the lab-repo family (badge block, Quickstart-first ordering, "See also" section linking the five sibling repos, License footer). The previous README pre-dated the Python port and described only the MATLAB workflow.
- **Internal version-string skew fixed.** Before this release, four different version strings were live simultaneously: `pyproject.toml` had `0.5.1`, `mumfordshah2d/__init__.py::__version__` had `"0.4.0"`, the `__init__.py` docstring claimed `0.2.0 — Phase 4 (beta)`, and `README_PYTHON.md` claimed `Phase 1 alpha (0.1.0)`. All four are now aligned at `0.5.2 — Phase 4 (beta)`.
- Andreas Weinmann's email updated to `andreas.weinmann@thws.de` in `CITATION.cff`.
- `CITATION.cff` populated with `version` and `date-released` fields so citation tooling reads from the file directly rather than falling back to GitHub Releases.

## 0.5.1 — 2026-05-08 (release fix — aarch64 cross-compile)

### Fixed

- **`pyo3` now built with the `abi3-py39` feature.** Without it, the
  release workflow's `linux (aarch64)` job failed with
  `Couldn't find any python interpreters` because maturin's QEMU
  cross-compile container has no aarch64 Python interpreter to detect.
  `abi3-py39` targets the stable Python ABI from 3.9 onward, removing
  the version-detection requirement and producing a single wheel that
  works on every Python ≥ 3.9 (matches the CSSD release pattern).
- The `0.5.0` tag exists on origin but produced no PyPI artifact
  (publish job skipped after the wheel job failed). `0.5.1` is the
  first version to actually reach PyPI.

## 0.5.0 — 2026-05-08 (Phase 7 — user-supplied prox + first inpainting demo)

Closes the v0.4.0 "Not yet supported" item. Users can now pass any of the
four MATLAB-canonical prox handles (or any custom callable with the same
contract) directly to `min_l2_mum_2d`. The Rust core's `Prox` trait is
now reachable from Python via a PyO3 callback bridge.

### Added

- **`prox=` kwarg on `min_l2_mum_2d`.** Accepts any callable with
  signature `prox(z, lambda) -> ndarray`, matching the factories in
  `mumfordshah2d.prox` (`make_prox_l2w`, `make_prox_l1w`, `make_prox_l0w`,
  `make_prox_inpaint`). The Python facade transparently bridges the
  user-facing array layout (`(rows, cols, channels)` or 2-D grayscale)
  to the Rust core's `(channels, rows, cols)` layout, so user-supplied
  prox factories can close over `f` in whatever shape they were
  constructed against.
- **`PyProx` Rust adapter** (in `src/lib.rs`) wrapping a `Py<PyAny>` as a
  `Prox` impl. One GIL acquisition per outer ADMM iteration. Acceptable
  for the four MATLAB-canonical proxes (closed-form NumPy ops); for
  inner loops the Rust-native `L2DataProx` (via `weights=`) remains
  available with no GIL crossings.
- **`demos_python/demo_inpainting.py`** — Python port of
  `Demos/demoInpainting.m`. Loads `fruitsColor.jpg`, generates a 60%
  random binary mask, restores via `make_prox_l2w(corrupted, mask)` +
  knight-move ADMM. Includes a `--smoke` mode for CI.

### Tested

- New `tests/test_prox_admm_integration.py` with 6 tests including a
  bit-identical regression check: `prox=make_prox_l2w(f, ones)` MUST
  produce the same array as the default `weights=None` path. Any drift
  in the PyProx bridge or the shape-shim layer surfaces immediately.
- Conflicts (`prox=` and `weights=` both passed) raise `ValueError` at
  the Python layer; non-callable `prox` raises `TypeError`.
- 155 pre-existing tests still pass (full suite minus matlab-parity
  which is host-shim-dependent).

### Security

- Built against `pyo3 >= 0.24.1` (closes GHSA-pph8-gcv7-4qj5,
  `PyString::from_object` buffer-overread). Also pins `numpy = 0.24`.

### Not yet supported (planned for 0.6.0+)

- **`makeProxL2Linop` / operator-stack prox** for deconvolution. Needs
  an iterative linear solver (CG/LSQR) per outer iteration and a
  warm-start contract on the prox `init` argument. Open questions on
  solver choice (SciPy callback vs. Rust-native `nalgebra-sparse`) and
  operator API (callable vs. matrix-free FFT for convolution).
- **Port of `Demos/demoDeconv*.m`** — blocked on the above.

## 0.4.0 — 2026-05-07 (Phase 6 — non-uniform weights + custom schedules)

Closes the last two MATLAB-feature-parity gaps in the Python facade. The
underlying Rust core already exposed both via traits/closures; this
release wires them through to Python kwargs.

### Added

- **`weights=` kwarg on `min_l2_mum_2d`.** Per-pixel non-negative
  weights, shape `(rows, cols)`, broadcast across channels. Mirrors
  MATLAB's `makeProxL2w(f, weights)`. Pixels with `weights[i,j] == 0`
  are unobserved (inpainting); the L2 prox at those pixels returns the
  consensus `z` directly. Validated for shape and non-negativity.
- **`mu_schedule=`, `nu_schedule=` kwargs**: pre-computed NumPy arrays
  of length ≤ `max_iter` indexed at iteration `k`. Shorter arrays clamp
  to the last entry (MATLAB semantics). Default schedules
  (`k².⁰¹·10⁻⁶` and `μ(k)·S/C(S,2)`) still apply when omitted. Useful
  for diagnostics, replicating published parameter sweeps, or warm
  starts. Validated for non-empty / finite / non-negative.
- **`L2DataProx::with_weights(f, weights)`** Rust constructor.

### Verified

19/19 live MATLAB-parity tests pass:
  - 16 prior cases (4-conn / 8-conn / knight-move × ρ × max_iter)
  - 2 weighted-prox cases (max_iter ∈ {20, 100})
  - 1 custom-μ-schedule case (linear ramp)

### Not yet supported (planned for 0.5.0+)

- **Custom `Prox` via Python callable.** The Rust `admm_l2_ms` accepts
  any `impl Prox`; the PyO3 wrapper is constrained to `L2DataProx`.
  Wiring an arbitrary Python callable through PyO3 requires per-iter
  GIL acquisition and careful error propagation — engineering effort
  on the order of the weights+schedule work but with smaller user demand
  (most users denoise; few do custom inverse problems). Defer until a
  concrete external use case lands.
- **L1 data fidelity** (`'method', 'L1'` in MATLAB), L0 hard-thresholding,
  and the inpaint-mask prox. The Auxiliary/makeProx*.m family except L2.

## 0.3.0 — 2026-05-07 (Phase 5 — full neighbourhood coverage)

Closes the gap between the Rust port and `mumfordShah2D.m`'s feature set:
the third (knight-move) neighbourhood lands and is parity-tested.

### Added

- **Knight-move neighbourhood** (`isotropic=2`, S=8). Uses the four
  primary stencils `[1,0]`, `[1,1]`, `[2,1]`, `[1,2]` with the
  MATLAB-default near-isotropic weights (ω ≈ {0.236, 0.115, 0.083, 0.083}).
  The closest-to-isotropic discretisation in the family.
- `admm.admm_knight_l2_ms` Rust entry point (thin wrapper over
  `admm_l2_ms` with `nhood_knight_move()`).
- `min_l2_mum_2d(..., isotropic=2)` Python kwarg.

### Verified

16/16 live MATLAB-parity cases pass (4 isotropic modes ✕ ρ ∈ {True, False}
× max_iter ∈ {20, 100}, plus the constant-image control). Tight
`atol=1e-9` everywhere except the known transient ρ=True/max_iter=20
case (`atol=1e-5` — converges to 1e-9 by max_iter=100).

### Not yet supported (planned for 0.4.0+)

- Custom μ schedules and arbitrary prox handles via the Python API.
  (The Rust `admm_l2_ms` already accepts `mu_seq: impl Fn(usize) -> f64`,
  `nu_seq: impl Fn(usize) -> f64`, and any `impl Prox`; only the PyO3
  wrapper hardcodes `default_mu_seq` and `L2DataProx`.)
- Non-uniform data weights (`weights ≠ 1` in the L2 prox; the Rust trait
  is general, the wrapper isn't).

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

### Added (post-initial draft of this CHANGELOG entry)

- **8-connected (near-isotropic) ADMM driver** (`isotropic=1`, S=4: rows
  + cols + NW-SE diagonals + NE-SW diagonals). Uses the MATLAB-default
  ω = {√2−1, 1−√2/2} weights for near-isotropic smoothing.
- **ρ-coupled ADMM variant** matching MATLAB's `rhoFlag=true` default.
  `min_l2_mum_2d(..., rho_coupling=True)` engages C(S, 2) inter-direction
  dual variables; `rho_coupling=False` is equivalent to the MATLAB
  `'nuSeq', @(k) 0` simplification.
- **`isotropic` and `rho_coupling` kwargs** on the Python facade.
- **Smoke notebook** at `demos_python/release_smoke_test.ipynb` exercises
  4-connected vs 8-connected, varying γ, and multi-channel RGB input.

### Verified end-to-end

12 live MATLAB-parity cases pass:
  - `test_step_4x4` × ρ ∈ {True, False} × max_iter ∈ {20, 100}: 4 cases
  - `test_step_4x4_isotropic_1` × ρ ∈ {True, False} × max_iter ∈ {20, 100}: 4 cases
  - `test_constant_image_4x4` × ρ ∈ {True, False} × max_iter ∈ {20, 100}: 4 cases
Tight `atol=1e-9` for ρ off and at convergence (max_iter=100); intermediate
`atol=1e-5` for ρ on at low iter counts (transient ADMM iterate divergence
from floating-point ordering, tightens to 1e-9 by max_iter=100).

### Not yet supported (planned for 0.3.0+)

- Knight-move neighbourhood (`isotropic=2`, S=8). Direction stencils [2,1]
  and [1,2] are already implemented in the direction processor; only the
  ADMM wiring is missing.
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
