# mumfordshah2d (Python / Rust)

Python / Rust port of the MATLAB / Java [MumfordShah2D](README.md) library.
Edge-preserving image restoration via the piecewise-smooth Mumford-Shah model
of Hohm, Storath, and Weinmann (Inverse Problems 31(11), 2015).

> **Status: 0.5.2 — Phase 4 (beta).** The 4-connected (anisotropic) 2-D
> Mumford-Shah ADMM driver is ported and verified at 1e-9 element-wise
> against MATLAB. 8-connected and ρ-coupled variants remain on the
> roadmap. See [`PORTED_BY.md`](PORTED_BY.md) for the phasing plan.

## Installation

From PyPI:

```bash
pip install mumfordshah2d
```

From source (requires Rust + maturin):

```bash
pip install maturin
maturin develop --release
```

Or as a regular editable install (also builds the Rust core):

```bash
pip install -e .
```

## Quick check

```python
import mumfordshah2d
print(mumfordshah2d.__version__)             # 0.1.0
print(mumfordshah2d.__original_authors__)    # Kilian Hohm, Martin Storath, Andreas Weinmann

import numpy as np
from mumfordshah2d import gauss_l2_mum_solve

# L2-Mumford-Shah within-segment smoothing on a noisy step:
y = np.array([0.0, 0.1, -0.05, 1.05, 0.95, 1.1])
mu = gauss_l2_mum_solve(y, alpha=2.0)
print(mu)
```

## What's currently exported

```python
# Utilities
soft_threshold, hard_threshold, expand_weights, rotate90, psnr

# Prox handles for the ADMM data-fidelity term
make_prox_l2w, make_prox_l1w, make_prox_l0w, make_prox_inpaint

# Phase 1 Rust primitive (exposed for testing)
gauss_l2_mum_solve

# Metadata
__version__, __original_authors__, __ported_by__
```

The full public API (`min_l2_l2_mumford_shah_2d` etc.) ships in Phase 4. Until
then, **use the original MATLAB code** in `mumfordShah2D.m` for end-to-end
restoration — that path is unchanged by this port.

## Array conventions

- **Python API**: grayscale arrays are `(rows, cols)`; colour/multichannel
  arrays are `(rows, cols, channels)` (numpy / imageio convention).
- **Internal Rust core**: `(channels, rows, cols)` (Java / MATLAB convention).
- A single conversion happens at the PyO3 boundary in `src/lib.rs`.

## Verification approach (Phases 2–5)

The Java `.class` files in `Java/bin/mumfordShah/` are kept on disk as a
reference oracle. From Phase 2 onwards, every Rust algorithm is fuzzed
against the Java implementation via a small `TestHarness.java` subprocess
shim, ensuring bit-equivalent (≤ 1e-12) results on hundreds of random
inputs.

## License

MIT, copyright Kilian Hohm, Martin Storath, Andreas Weinmann (original) and
the Claude Sonnet coding agent contribution (Anthropic, 2026, port). See
[`LICENSE`](LICENSE) and [`PORTED_BY.md`](PORTED_BY.md).
