# Ported from MATLAB Auxiliary/{makeProxL2w, makeProxL1w, makeProxL0w,
# makeProxInpaint}.m (MumfordShah2D — Hohm/Storath/Weinmann) by Claude
# Sonnet coding agent, Anthropic, 2026.
"""Prox operators for the data-fidelity term of the outer ADMM loop.

Each ``make_prox_*`` factory returns a callable ``prox(z, tau, init=None)``
that solves the corresponding proximal subproblem in closed form. The
signature mirrors the MATLAB ``handle = @(v, tau, opts) ...`` convention but
exposes ``init`` as a keyword instead of an opaque ``opts`` struct (the
denoising / L1 / L0 proxes ignore ``init`` entirely; ``make_prox_l2_linop``,
which takes a real solver, will land in Phase 5).
"""

from __future__ import annotations

from typing import Callable, Optional

import numpy as np

from mumfordshah2d.utils import expand_weights, hard_threshold, soft_threshold


ProxHandle = Callable[..., np.ndarray]


def make_prox_l2w(f: np.ndarray, weights: np.ndarray) -> ProxHandle:
    """Weighted-L2 (denoising) prox.

    Solves ``min_u  (1/2) Σ w_i (u_i - f_i)² + (tau/2) Σ (u_i - z_i)²``::

        u = (w * f + tau * z) / (w + tau)

    Direct port of ``Auxiliary/makeProxL2w.m``.
    """
    f = np.asarray(f, dtype=np.float64)
    w = expand_weights(np.asarray(weights, dtype=np.float64), f)

    def prox(z: np.ndarray, tau: float, init: Optional[np.ndarray] = None) -> np.ndarray:
        del init  # unused — the closed form needs no warm start
        z = np.asarray(z, dtype=np.float64)
        return (w * f + tau * z) / (w + tau)

    return prox


def make_prox_l1w(f: np.ndarray, weights: np.ndarray) -> ProxHandle:
    """Weighted-L1 prox.

    Solves ``min_u Σ w_i |u_i - f_i| + (tau/2) Σ (u_i - z_i)²``::

        u = soft_threshold(z - f, w / tau) + f

    Direct port of ``Auxiliary/makeProxL1w.m``.
    """
    f = np.asarray(f, dtype=np.float64)
    w = expand_weights(np.asarray(weights, dtype=np.float64), f)

    def prox(z: np.ndarray, tau: float, init: Optional[np.ndarray] = None) -> np.ndarray:
        del init
        z = np.asarray(z, dtype=np.float64)
        return soft_threshold(z - f, w / tau) + f

    return prox


def make_prox_l0w(f: np.ndarray, weights: np.ndarray) -> ProxHandle:
    """Weighted-L0 (hard-thresholding) prox.

    Solves ``min_u Σ w_i 1{u_i ≠ f_i} + (tau/2) Σ (u_i - z_i)²``::

        u = hard_threshold(z - f, sqrt(2 w / tau)) + f

    Direct port of ``Auxiliary/makeProxL0w.m``.
    """
    f = np.asarray(f, dtype=np.float64)
    w = expand_weights(np.asarray(weights, dtype=np.float64), f)

    def prox(z: np.ndarray, tau: float, init: Optional[np.ndarray] = None) -> np.ndarray:
        del init
        z = np.asarray(z, dtype=np.float64)
        return hard_threshold(z - f, np.sqrt(2.0 * w / tau)) + f

    return prox


def make_prox_inpaint(f: np.ndarray, mask: np.ndarray) -> ProxHandle:
    """Inpainting prox: clamp known pixels to ``f``, leave unknowns to ``z``.

    Direct port of ``Auxiliary/makeProxInpaint.m``::

        u = f .* mask + z .* ~mask

    The MATLAB code passes ``mask`` through ``plExpandWeights`` first; we
    do the same so a 2-D mask broadcasts to a multi-channel image.
    """
    f = np.asarray(f, dtype=np.float64)
    m = expand_weights(np.asarray(mask, dtype=np.float64), f)
    # MATLAB's `~mask` is a logical-not; we want the boolean complement of
    # whatever the mask "is non-zero" predicate yields.
    m_bool = m.astype(bool)

    def prox(z: np.ndarray, tau: float, init: Optional[np.ndarray] = None) -> np.ndarray:
        del init, tau  # both unused — the prox is just a masked overwrite
        z = np.asarray(z, dtype=np.float64)
        return np.where(m_bool, f, z)

    return prox
