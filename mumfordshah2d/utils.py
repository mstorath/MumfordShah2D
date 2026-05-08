# Ported from MATLAB Auxiliary/{softThreshold, hardThreshold, expandWeights,
# rotate90, plpsnr}.m (MumfordShah2D — Hohm/Storath/Weinmann) by Claude
# Sonnet coding agent, Anthropic, 2026.
"""Small array utilities used across the package.

Each function corresponds one-for-one to a MATLAB helper in
``MumfordShah2D/Auxiliary/`` and follows the same semantics. Where the
MATLAB version is dimension-flexible we accept the same shapes and return
arrays with the same broadcasting behaviour.
"""

from __future__ import annotations

import numpy as np


def soft_threshold(x: np.ndarray, tau: float | np.ndarray) -> np.ndarray:
    """Element-wise soft threshold ``sign(x) * max(|x| - tau, 0)``.

    Port of ``Auxiliary/softThreshold.m`` (the MATLAB code names the
    function ``plSoftThreshold`` in the body — same routine).
    """
    x = np.asarray(x, dtype=np.float64)
    tau_arr = np.asarray(tau, dtype=np.float64)
    return np.maximum(np.abs(x) - tau_arr, 0.0) * np.sign(x)


def hard_threshold(x: np.ndarray, tau: float | np.ndarray) -> np.ndarray:
    """Element-wise hard threshold ``x * (|x| >= tau)``.

    Port of ``Auxiliary/hardThreshold.m``.
    """
    x = np.asarray(x, dtype=np.float64)
    tau_arr = np.asarray(tau, dtype=np.float64)
    return x * (np.abs(x) >= tau_arr)


def expand_weights(weights: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Replicate a weight map along the channel axis of a multichannel image.

    Mirrors ``Auxiliary/expandWeights.m``::

        if size(f, 3) ~= size(weights, 3)
            weights = repmat(weights, [1 1 size(f, 3)]);
        end

    In MATLAB ``size(any_array, 3)`` returns 1 for a 1-D vector or 2-D
    matrix, so this is effectively a no-op unless ``f`` is genuinely 3-D
    (multichannel) and ``weights`` is not yet aligned to its channel
    count. We mirror the same semantics for arrays of any dimensionality.

    Parameters
    ----------
    weights : ndarray
        Weight map. Any ndim ≥ 1 is accepted.
    f : ndarray
        Reference array whose channel count we want to match. ``ndim < 3``
        is treated as "1 channel"; ``ndim == 3`` uses ``f.shape[2]``.

    Returns
    -------
    ndarray
        Either ``weights`` unchanged (when the channel counts already match)
        or ``weights`` replicated along a freshly added trailing axis to
        match ``f.shape[2]``.
    """
    weights = np.asarray(weights, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    # MATLAB: size(arr, 3) returns 1 if arr has fewer than 3 dims.
    f_channels = f.shape[2] if f.ndim >= 3 else 1
    w_channels = weights.shape[2] if weights.ndim >= 3 else 1
    if w_channels == f_channels:
        return weights
    if w_channels == 1:
        # Promote to (..., f_channels): add a trailing axis if needed,
        # then replicate.
        if weights.ndim < 3:
            weights = weights[..., np.newaxis]
        return np.repeat(weights, f_channels, axis=2)
    raise ValueError(
        f"Cannot expand weights with {w_channels} channels to match image "
        f"with {f_channels} channels."
    )


def rotate90(f: np.ndarray, k: int) -> np.ndarray:
    """Rotate a ``(channels, rows, cols)`` array by ``k * 90`` degrees per
    channel.

    Port of ``Auxiliary/rotate90.m``. The MATLAB version requires a 3-D
    input (``ndims(f) == 3``) and applies ``rot90`` to each ``squeeze(f(i,
    :, :))``. We do the same with ``np.rot90`` along axes ``(1, 2)``.
    """
    f = np.asarray(f, dtype=np.float64)
    if f.ndim != 3:
        raise ValueError(f"rotate90 expects a 3-D array; got ndim={f.ndim}")
    return np.rot90(f, k=k, axes=(1, 2)).copy()


def psnr(ground_truth: np.ndarray, signal: np.ndarray) -> float:
    """Peak signal-to-noise ratio in dB.

    Port of ``Auxiliary/plpsnr.m``::

        mse = mean(|signal - ground_truth|^2)
        psnr = 10 * log10(max(ground_truth)^2 / mse)

    Note that the MATLAB version uses the per-pixel max of ``ground_truth``
    (not the conventional dynamic range), which we preserve.
    """
    gt = np.asarray(ground_truth, dtype=np.float64)
    sig = np.asarray(signal, dtype=np.float64)
    mse = float(np.mean(np.abs(sig.ravel() - gt.ravel()) ** 2))
    if mse == 0.0:
        return float("inf")
    return 10.0 * float(np.log10((gt.max() ** 2) / mse))
