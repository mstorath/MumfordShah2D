# Port Attribution

This Python/Rust package (`mumfordshah2d`) is an **automated port** of the original
MATLAB/Java *MumfordShah2D* library:

| | |
|---|---|
| **Original authors** | Kilian Hohm, Martin Storath, Andreas Weinmann |
| **Original license** | MIT |
| **Original language** | MATLAB (high-level), Java (performance core) |

## Port details

| | |
|---|---|
| **Port author** | Claude Sonnet coding agent (Anthropic, 2026) |
| **Port language** | Python (API layer) + Rust/PyO3 (performance core) |
| **Port license** | MIT (unchanged) |

## Phasing

The port is being delivered in five phases. Each phase ships independently on
its own `claude/mumfordshah-port-phase-N-2026-MM` branch.

| Phase | Scope | Status |
|---|---|---|
| 1 | Scaffolding + small primitives (`Point`, `GaussElim`/`GaussL2Mum`, `Image`, MATLAB utility ports, prox handles) | **shipped (this commit set)** |
| 2 | `TautString` (1D weighted-TV via lower-envelope) | pending |
| 3 | 1D Mumford-Shah DP solvers (L1 and L2 variants) | pending |
| 4 | 2D direction processor + outer ADMM loop | pending |
| 5 | Inverse problems (deconv / inpaint / Radon), demos, full test suite | pending |

## What was ported in Phase 1

- **`src/image.rs::Point`** — Direct port of `Java/src/mumfordShah/Point.java`
  (Kilian Hohm). Plain struct holding `(position, value, gradient)` for
  `TautString` book-keeping in later phases.
- **`src/image.rs::Image`** — Direct port of `Java/src/mumfordShah/Image.java`
  (Hohm/Storath/Weinmann). 3D `(channels, rows, cols)` image with row, column,
  and direction-slice accessors. Used by the 2D ADMM in Phase 4.
- **`src/gauss_elim.rs::GaussElim` trait + `GaussL2Mum` impl** — Direct port of
  `Java/src/mumfordShah/GaussElim.java` and `GaussL2Mum.java`
  (Hohm/Storath/Weinmann). Pre-computed tridiagonal LU factorisations for the
  L2-Mumford-Shah within-segment smoothing problem
  `min_μ Σ(y_i - μ_i)² + α Σ(μ_{i+1} - μ_i)²`.
- **`mumfordshah2d/utils.py`** — Port of `Auxiliary/{softThreshold,
  hardThreshold, expandWeights, rotate90, plpsnr}.m`.
- **`mumfordshah2d/prox.py`** — Port of `Auxiliary/{makeProxL2w, makeProxL1w,
  makeProxL0w, makeProxInpaint}.m`. The remaining proxes (`makeProxL2Linop`)
  ship in Phase 5 with the operator stack.

## What changed

- API names follow Python conventions: `make_prox_l2w` instead of `makeProxL2w`,
  `soft_threshold` instead of `plSoftThreshold`.
- Data is passed as NumPy arrays (`float64`) instead of MATLAB matrices.
- The Python API uses `(rows, cols)` for grayscale and `(rows, cols, channels)`
  for colour/multichannel images (numpy/imageio convention). The internal Rust
  representation is `(channels, rows, cols)` (Java/MATLAB convention); a single
  conversion happens at the PyO3 boundary in `src/lib.rs`.
- Performance-critical loops are in Rust (replacing the Java JAR).
- The Python API uses keyword-only arguments with the same defaults as MATLAB.

## What this port does NOT change

- The original MATLAB and Java sources are untouched. Both code paths
  remain runnable with MATLAB + the existing `pottslab-standalone.jar`.
- The original MATLAB demos in `Demos/*.m` work exactly as before.
- The repository LICENSE (MIT, copyright Hohm/Storath/Weinmann) is unchanged.

## Academic reference (from original authors)

K. Hohm, M. Storath, A. Weinmann.
*An algorithmic framework for Mumford–Shah regularization of inverse problems
in imaging.* Inverse Problems 31 (11), 115011, 2015.
[DOI 10.1088/0266-5611/31/11/115011](https://doi.org/10.1088/0266-5611/31/11/115011).
[Preprint PDF](http://bigwww.epfl.ch/publications/hohm1501.pdf).
