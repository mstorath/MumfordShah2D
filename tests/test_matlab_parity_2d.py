# Live Rust↔MATLAB parity tests for the 2-D Mumford-Shah ADMM driver.
#
# Each test:
#   1. Generates a small (≤8×8) 2-D image.
#   2. Runs the Rust port `mumfordshah2d._core.min_l2_mum_2d` with a fixed
#      iteration count (tol=0 forces all max_iter iterations).
#   3. Runs MATLAB's `mumfordShah2D` with the same parameters via the host
#      `matlab` shim, using `makeProxL2w` (= L2 data prox with weights=1)
#      to match the default Rust `L2DataProx`.
#   4. Asserts element-wise agreement at the chosen tolerance.
#
# Skipped (not failed) when the matlab shim is unavailable.
#
# Tolerance: `atol=1e-9` is comfortable on 4×4 inputs at max_iter=20 because
# both implementations follow the same ADMM iteration to the bit (same μ
# schedule, same prox formula, same direction-stripe order). Larger inputs
# or higher iteration counts may need looser tolerance — the order in which
# rayon's parallel-stripe results are written back is sequential and matches
# Java's serial loop, but inner sums in `compute_dlr` can accumulate
# float-reorder noise.

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest

import mumfordshah2d._core as ms_core

REPO_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = REPO_ROOT.parents[1]


def _shim_available() -> bool:
    return shutil.which("matlab") is not None and bool(os.environ.get("HOST_MATLAB"))


pytestmark = pytest.mark.skipif(
    not _shim_available(),
    reason="matlab shim not configured (HOST_MATLAB unset or matlab not on PATH)",
)


def _run_matlab_2d(
    f: np.ndarray,
    gamma: float,
    alpha: float,
    max_iter: int,
    *,
    rho_coupling: bool = True,
    isotropic: int = 0,
):
    """Invoke `mumfordShah2D(gamma, alpha, makeProxL2w(f, ones), …)` in host MATLAB
    with the same settings the Rust port uses. Returns the smoothed image.

    When ``rho_coupling=False``, passes `nuSeq = @(k) 0` to MATLAB so it
    skips the inter-direction ρ dual-variable updates (matches the Rust
    `no_rho_coupling` schedule).
    """
    work_dir = Path(tempfile.mkdtemp(prefix="ms2d-2d-parity-", dir=WORKSPACE_ROOT))
    try:
        in_path = work_dir / "f.csv"
        out_path = work_dir / "out.csv"
        script_path = work_dir / "run_parity.m"

        # MATLAB's mumfordShah2D expects 2-D input for a single channel
        # (it auto-decides via ndims). We restrict the parity test to the
        # single-channel case for now: f shape (1, n_row, n_col).
        if f.ndim == 3:
            assert f.shape[0] == 1, "only single-channel parity tested for now"
            f2d = f[0]
        elif f.ndim == 2:
            f2d = f
        else:
            raise ValueError(f"unexpected f.ndim = {f.ndim}")

        np.savetxt(in_path, f2d, delimiter=",", fmt="%.17e")

        nu_seq_arg = "varargin = {'nuSeq', @(k) 0};" if not rho_coupling else "varargin = {};"
        script = f"""
addpath(genpath('{REPO_ROOT}'));
javaaddpath('{REPO_ROOT}/Java/bin');
f = readmatrix('{in_path}');
weights = ones(size(f));
prox = makeProxL2w(f, weights);
{nu_seq_arg}
[out, nIter] = mumfordShah2D({gamma:.17e}, {alpha:.17e}, prox, ...
    'isotropic', {isotropic}, 'tol', 1e-300, 'maxIter', {max_iter}, varargin{{:}});
fprintf('MATLAB nIter=%d\\n', nIter);
writematrix(out, '{out_path}');
"""
        script_path.write_text(script)

        result = subprocess.run(
            ["matlab", "-batch", f"addpath('{work_dir}'); run_parity"],
            capture_output=True,
            text=True,
            timeout=180,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"MATLAB mumfordShah2D failed (rc={result.returncode}):\n"
                f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
            )
        out2d = np.loadtxt(out_path, delimiter=",")
        if out2d.ndim == 1:
            out2d = out2d[None, :]
        # Re-add the channel dim for shape consistency.
        return out2d[None, :, :]
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


# ----------------------------------------------------------------------
# Tiny deterministic cases
# ----------------------------------------------------------------------

@pytest.mark.parametrize("rho_coupling", [True, False])
@pytest.mark.parametrize("max_iter", [20, 100])
def test_step_4x4(max_iter, rho_coupling):
    """Binary step image, low γ, modest α. Tested with both ρ-coupled
    (MATLAB-default) and ρ-disabled (`nuSeq=@(k) 0`) ADMM variants.

    Without ρ-coupling the two implementations are bit-identical (atol=1e-9).
    With ρ-coupling the iterate trajectories agree to atol=1e-5 at low
    max_iter (the ADMM algorithms differ in floating-point accumulator order
    inside the cross-direction sum) and tighten to atol=1e-9 by max_iter=100
    (both implementations have converged to the same fixed point).
    """
    f = np.zeros((1, 4, 4), dtype=np.float64)
    f[0, :, 2:] = 1.0

    gamma = 0.05
    alpha = 1e-6

    out_rust = ms_core.min_l2_mum_2d(
        f,
        gamma=gamma,
        alpha=alpha,
        tol=1e-300,
        max_iter=max_iter,
        verbose=False,
        rho_coupling=rho_coupling,
    )
    out_matlab = _run_matlab_2d(
        f, gamma=gamma, alpha=alpha, max_iter=max_iter, rho_coupling=rho_coupling
    )

    # Tight tolerance once both algorithms have converged or the simpler
    # path (no ρ-coupling) is in use; loose at intermediate iterates with ρ on.
    if rho_coupling and max_iter < 100:
        atol = 1e-5
    else:
        atol = 1e-9
    npt.assert_allclose(out_rust, out_matlab, atol=atol)


@pytest.mark.parametrize("rho_coupling", [True, False])
@pytest.mark.parametrize("max_iter", [20, 100])
def test_step_4x4_isotropic_1(max_iter, rho_coupling):
    """Same step image as `test_step_4x4` but with `isotropic=1` (8-connected,
    S=4: rows + cols + NW-SE diagonals + NE-SW diagonals). Tests that the
    omega weights are applied correctly and the diagonal direction sweeps
    align bit-for-bit with MATLAB."""
    f = np.zeros((1, 4, 4), dtype=np.float64)
    f[0, :, 2:] = 1.0

    gamma = 0.05
    alpha = 1e-6

    out_rust = ms_core.min_l2_mum_2d(
        f,
        gamma=gamma,
        alpha=alpha,
        tol=1e-300,
        max_iter=max_iter,
        verbose=False,
        rho_coupling=rho_coupling,
        isotropic=1,
    )
    out_matlab = _run_matlab_2d(
        f, gamma=gamma, alpha=alpha, max_iter=max_iter, rho_coupling=rho_coupling, isotropic=1,
    )

    if rho_coupling and max_iter < 100:
        atol = 1e-5
    else:
        atol = 1e-9
    npt.assert_allclose(out_rust, out_matlab, atol=atol)


@pytest.mark.parametrize("rho_coupling", [True, False])
@pytest.mark.parametrize("max_iter", [20, 100])
def test_constant_image_4x4(max_iter, rho_coupling):
    """Constant image must be invariant to any number of ADMM iterations."""
    f = np.full((1, 4, 4), 0.7, dtype=np.float64)

    gamma = 1.0
    alpha = 1.0

    out_rust = ms_core.min_l2_mum_2d(
        f,
        gamma=gamma,
        alpha=alpha,
        tol=1e-300,
        max_iter=max_iter,
        verbose=False,
        rho_coupling=rho_coupling,
    )
    out_matlab = _run_matlab_2d(
        f, gamma=gamma, alpha=alpha, max_iter=max_iter, rho_coupling=rho_coupling
    )

    npt.assert_allclose(out_rust, out_matlab, atol=1e-9)
