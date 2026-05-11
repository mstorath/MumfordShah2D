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
print(mumfordshah2d.__version__)             # 0.5.2
print(mumfordshah2d.__original_authors__)    # Kilian Hohm, Martin Storath, Andreas Weinmann

import numpy as np
from mumfordshah2d import min_l2_mum_2d

# 2-D L2 Mumford-Shah ADMM on a noisy step image:
H, W = 32, 32
f = np.zeros((H, W))
f[:, W // 2:] = 1.0
f += np.random.default_rng(0).normal(0, 0.1, size=f.shape)

u = min_l2_mum_2d(f, gamma=0.5, alpha=1.0)
print(u.shape, u.dtype)
```

## Currently exported

```python
# 2-D ADMM driver (port of mumfordShah2D.m)
min_l2_mum_2d(f, gamma, alpha, *, tol=1e-3, max_iter=50000)

# Utilities (port of MATLAB Auxiliary/*.m)
soft_threshold, hard_threshold, expand_weights, rotate90, psnr

# Prox handles for the ADMM data-fidelity term (port of Auxiliary/makeProx*.m)
make_prox_l2w, make_prox_l1w, make_prox_l0w, make_prox_inpaint

# Rust primitives (exposed for testing)
gauss_l2_mum_solve, gauss_l2_mum_cost

# Metadata
__version__, __original_authors__, __ported_by__
```

The 8-connected (near-isotropic) and ρ-coupled variants from the original
MATLAB driver have not been ported yet — for those, fall back to the
MATLAB code in `mumfordShah2D.m`.

## Array conventions

- **Python API**: grayscale arrays are `(rows, cols)`; colour/multichannel
  arrays are `(rows, cols, channels)` (numpy / imageio convention).
- **Internal Rust core**: `(channels, rows, cols)` (Java / MATLAB convention).
- A single conversion happens at the PyO3 boundary in `src/lib.rs`.

## Verification approach

The Java `.class` files in `Java/bin/mumfordShah/` are kept on disk as a
reference oracle. Every Rust algorithm is fuzzed against the Java
implementation via a small `TestHarness.java` subprocess shim, ensuring
bit-equivalent (≤ 1e-12) results on hundreds of random inputs.

## License

MIT, copyright Kilian Hohm, Martin Storath, Andreas Weinmann (original) and
the Claude Sonnet coding agent contribution (Anthropic, 2026, port). See
[`LICENSE`](LICENSE) and [`PORTED_BY.md`](PORTED_BY.md).
