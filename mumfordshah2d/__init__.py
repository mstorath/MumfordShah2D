# Ported from MATLAB/Java MumfordShah2D (Hohm, Storath & Weinmann) by
# Claude Sonnet coding agent, Anthropic, 2026.
"""mumfordshah2d — Edge-preserving image restoration via the Mumford-Shah model.

Python / Rust port of the original MATLAB / Java MumfordShah2D library by
Kilian Hohm, Martin Storath, and Andreas Weinmann. Time-critical kernels are
implemented in Rust (via PyO3) and compiled as a C extension; the high-level
ADMM driver and helper functions are pure Python / NumPy / SciPy.

Status
------
**0.1.0 — Phase 1 (alpha).** Scaffolding plus small primitives only. The
high-level 2D solver lands in Phase 4. See ``PORTED_BY.md`` for the full
phasing plan.

Currently exported
------------------
Utilities (port of MATLAB ``Auxiliary/*.m``)::

    soft_threshold, hard_threshold, expand_weights, rotate90, psnr

Prox handles for the ADMM data-fidelity term (port of
``Auxiliary/makeProx*.m``)::

    make_prox_l2w   — denoising / weighted L2
    make_prox_l1w   — robust / weighted L1
    make_prox_l0w   — hard-thresholding / weighted L0
    make_prox_inpaint — inpainting from a mask

Phase 1 Rust primitive::

    gauss_l2_mum_solve(y, alpha)   — L2-MS within-segment smoothing
    gauss_l2_mum_cost(y, alpha)    — corresponding cost

Reference
---------
K. Hohm, M. Storath, A. Weinmann.
*An algorithmic framework for Mumford–Shah regularization of inverse problems
in imaging.* Inverse Problems 31 (11), 115011, 2015.
DOI 10.1088/0266-5611/31/11/115011.
"""

__version__ = "0.1.0"
__original_authors__ = "Kilian Hohm, Martin Storath, Andreas Weinmann"
__ported_by__ = (
    "Claude Sonnet coding agent (Anthropic, 2026) — "
    "automated port from MATLAB/Java MumfordShah2D"
)

from mumfordshah2d.utils import (
    soft_threshold,
    hard_threshold,
    expand_weights,
    rotate90,
    psnr,
)
from mumfordshah2d.prox import (
    make_prox_l2w,
    make_prox_l1w,
    make_prox_l0w,
    make_prox_inpaint,
)

# The compiled Rust core. If the extension isn't built yet, fall back to the
# pure-Python implementation in `_core_py` so the package still imports —
# tests will pick up the fallback automatically.
try:
    from mumfordshah2d._core import gauss_l2_mum_solve, gauss_l2_mum_cost
    _USING_RUST_CORE = True
except ImportError:  # pragma: no cover — only hit before `maturin develop`
    from mumfordshah2d._core_py import gauss_l2_mum_solve, gauss_l2_mum_cost
    _USING_RUST_CORE = False

__all__ = [
    # Utilities
    "soft_threshold",
    "hard_threshold",
    "expand_weights",
    "rotate90",
    "psnr",
    # Prox handles
    "make_prox_l2w",
    "make_prox_l1w",
    "make_prox_l0w",
    "make_prox_inpaint",
    # Rust primitive
    "gauss_l2_mum_solve",
    "gauss_l2_mum_cost",
    # Metadata
    "__version__",
    "__original_authors__",
    "__ported_by__",
]
