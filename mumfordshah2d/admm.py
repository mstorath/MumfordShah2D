# Ported from MATLAB/Java MumfordShah2D (Hohm, Storath & Weinmann) by
# Claude Sonnet coding agent, Anthropic, 2026.
"""High-level ADMM driver for the 2-D L2 Mumford-Shah model.

Wraps the Rust core ``mumfordshah2d._core.min_l2_mum_2d`` with input
validation, dimension auto-detection (2-D grayscale vs 3-D multi-channel),
and a friendlier kwarg-only signature.

Currently implements only the 4-connected (anisotropic) variant with the
default L2 data-fidelity prox (``makeProxL2w`` analogue at uniform weights).
The 8-connected (near-isotropic) variant and the rho-coupled MATLAB-default
behaviour are deferred to a later release.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray


def min_l2_mum_2d(
    f: ArrayLike,
    gamma: float,
    alpha: float,
    *,
    tol: float = 1e-3,
    max_iter: int = 50000,
    verbose: bool = False,
    rho_coupling: bool = True,
    isotropic: int = 0,
    weights: Optional[ArrayLike] = None,
    mu_schedule: Optional[ArrayLike] = None,
    nu_schedule: Optional[ArrayLike] = None,
) -> NDArray[np.float64]:
    """Edge-preserving image restoration via the L2 Mumford-Shah model.

    Solves the consensus ADMM splitting of the discrete L2 Mumford-Shah
    objective with a Potts (jump) penalty and quadratic within-segment
    smoothness:

        min_u  γ |J|  +  α Σ_{i ∉ J} (u_{i+1} − u_i)²
                       +  Σ_i ||u_i − f_i||²

    on a 4-connected (rows + columns) grid. Within-segment smoothness is
    solved exactly per stripe via the cached tridiagonal LU; the partition
    via a 1-D dynamic program.

    Parameters
    ----------
    f
        Input image, shape ``(rows, cols)`` for grayscale or
        ``(rows, cols, channels)`` for multi-channel. Float64 internally.
    gamma
        Potts penalty weight (jump count cost). Larger ⇒ fewer jumps.
    alpha
        Within-segment smoothness weight. Larger ⇒ smoother segments.
        ``alpha = 0`` reduces to the L2 Potts (piecewise-constant) model.
    tol
        Stopping tolerance on the relative discrepancy
        ``‖u₀ − u₁‖ / (‖u₀‖ + ‖u₁‖)`` between the two ADMM directions.
    max_iter
        Iteration cap. The MATLAB default is 50000.
    verbose
        Print per-iteration ``(μ, error, γ', α')`` to stderr.

    Returns
    -------
    out
        Smoothed image, same shape and dtype convention as ``f``.

    Notes
    -----
    The Rust core uses a fixed μ schedule ``μ(k) = k².⁰¹ · 10⁻⁶`` (the
    MATLAB default) and a uniform L2 data-fidelity prox
    ``prox(z, λ) = (f + λ z) / (1 + λ)``. Custom μ schedules and weighted
    prox handles are deferred to a future release.
    """
    arr = np.asarray(f, dtype=np.float64)
    if arr.ndim == 2:
        # (rows, cols) → (1, rows, cols)
        arr3 = arr[np.newaxis, :, :]
        squeeze = True
    elif arr.ndim == 3:
        # (rows, cols, channels) → (channels, rows, cols)
        arr3 = np.transpose(arr, (2, 0, 1))
        squeeze = False
    else:
        raise ValueError(
            f"input must be 2-D (grayscale) or 3-D (multi-channel); got {arr.ndim}-D"
        )
    arr3 = np.ascontiguousarray(arr3)

    if not np.isfinite(arr3).all():
        raise ValueError("input contains NaN or Inf")
    if gamma < 0:
        raise ValueError(f"gamma must be ≥ 0, got {gamma}")
    if alpha < 0:
        raise ValueError(f"alpha must be ≥ 0, got {alpha}")
    if tol <= 0:
        raise ValueError(f"tol must be > 0, got {tol}")
    if max_iter < 1:
        raise ValueError(f"max_iter must be ≥ 1, got {max_iter}")

    # Late-bind the Rust core so the package still imports if the extension
    # isn't yet built (a fallback path in __init__.py routes 1-D primitives
    # through pure Python; the 2-D ADMM has no pure-Python fallback yet).
    try:
        from mumfordshah2d._core import min_l2_mum_2d as _rust_min_l2_mum_2d
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "mumfordshah2d._core is not built; run `maturin develop --release` "
            "from the repo root before calling min_l2_mum_2d."
        ) from e

    if isotropic not in (0, 1, 2):
        raise ValueError(
            f"isotropic must be 0 (4-connected), 1 (8-connected), or 2 (knight-move), got {isotropic}"
        )

    weights_arr = None
    if weights is not None:
        weights_arr = np.ascontiguousarray(np.asarray(weights, dtype=np.float64))
        if weights_arr.shape != arr3.shape[1:]:
            raise ValueError(
                f"weights must have shape (rows, cols) = {arr3.shape[1:]}, "
                f"got {weights_arr.shape}"
            )
        if (weights_arr < 0).any():
            raise ValueError("weights must be non-negative")

    def _validate_schedule(name: str, val: Optional[ArrayLike]) -> Optional[NDArray[np.float64]]:
        if val is None:
            return None
        a = np.ascontiguousarray(np.asarray(val, dtype=np.float64).ravel())
        if a.size == 0:
            raise ValueError(f"{name} must be non-empty")
        if not np.isfinite(a).all():
            raise ValueError(f"{name} contains NaN or Inf")
        if (a < 0).any():
            raise ValueError(f"{name} must be non-negative")
        return a

    mu_schedule_arr = _validate_schedule("mu_schedule", mu_schedule)
    nu_schedule_arr = _validate_schedule("nu_schedule", nu_schedule)

    out3 = _rust_min_l2_mum_2d(
        arr3,
        gamma=float(gamma),
        alpha=float(alpha),
        tol=float(tol),
        max_iter=int(max_iter),
        verbose=bool(verbose),
        rho_coupling=bool(rho_coupling),
        isotropic=int(isotropic),
        weights=weights_arr,
        mu_schedule=mu_schedule_arr,
        nu_schedule=nu_schedule_arr,
    )
    if squeeze:
        return out3[0]
    return np.transpose(out3, (1, 2, 0))


__all__ = ["min_l2_mum_2d"]
