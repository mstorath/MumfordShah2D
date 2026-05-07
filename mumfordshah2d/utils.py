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
    """Broadcast a 2-D weight map across the channel axis of a 3-D image.

    Mirrors ``Auxiliary/expandWeights.m`` (alias ``plExpandWeights``). The
    MATLAB version replicates the weights along ``size(f, 3)`` if the
    channel counts disagree. We replicate to ``f.shape[2]``.

    Parameters
    ----------
    weights : ndarray
        Either ``(rows, cols)`` or ``(rows, cols, channels)``.
    f : ndarray
        Reference image whose channel count we match. Either ``(rows, cols)``
        (treated as 1 channel) or ``(rows, cols, channels)``.

    Returns
    -------
    ndarray
        ``(rows, cols, channels)`` weight map.
    """
    weights = np.asarray(weights, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    f_channels = 1 if f.ndim == 2 else f.shape[2]
    w_channels = 1 if weights.ndim == 2 else weights.shape[2]
    if w_channels == f_channels:
        if weights.ndim == 2 and f_channels == 1:
            return weights[..., np.newaxis]
        return weights
    if weights.ndim == 2:
        return np.repeat(weights[..., np.newaxis], f_channels, axis=2)
    if w_channels == 1:
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
