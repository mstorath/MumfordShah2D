// Ported from Java/src/mumfordShah/ (Hohm, Storath & Weinmann)
// by Claude Sonnet coding agent, Anthropic, 2026.
//
// PyO3 module: exposes Phase 1 Rust primitives to Python.

mod image;
mod gauss_elim;
mod mumfordshah_1d;
mod direction_processor;
mod admm;

use numpy::{IntoPyArray, PyArray1, PyArray3, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3};
use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::admm::{
    admm_4connected_l2_ms, admm_8connected_l2_ms, admm_knight_l2_ms, default_mu_seq,
    default_nu_seq, no_rho_coupling, L2DataProx, Prox,
};
use crate::gauss_elim::{GaussElim, GaussL2Mum};
use crate::image::Image;

/// Adapter that lets a Python callable act as a Rust `Prox`.
///
/// On each call to `apply` the GIL is acquired, the consensus image `z`
/// is exposed as a `(channels, rows, cols)` numpy array, and the caller
/// is invoked as `callable(z, lambda)` exactly like the Python prox
/// factories return (see `mumfordshah2d/prox.py`). The returned numpy
/// array is copied into `out`.
///
/// One GIL crossing per outer ADMM iteration. Acceptable for the four
/// MATLAB-canonical proxes (closed-form NumPy ops); a Rust-native
/// dispatch is available for performance via `weights=` (which builds
/// `L2DataProx` directly without going through Python).
struct PyProx {
    callable: Py<PyAny>,
}

impl Prox for PyProx {
    fn apply(&self, z: &Image, lambda: f64, out: &mut Image) {
        Python::with_gil(|py| {
            let z_owned = z.data.clone();
            let z_arr = z_owned.into_pyarray(py);
            let result = self
                .callable
                .call1(py, (z_arr, lambda))
                .expect("python prox raised");
            let bound = result.bind(py);
            let pyarr = bound
                .downcast::<PyArray3<f64>>()
                .expect("python prox must return a float64 ndarray of shape (channels, rows, cols)");
            let readonly = pyarr.readonly();
            let view = readonly.as_array();
            assert_eq!(
                view.shape(),
                out.data.shape(),
                "python prox returned shape {:?}, expected {:?}",
                view.shape(),
                out.data.shape()
            );
            out.data.assign(&view);
        });
    }
}

/// Dispatch on `(isotropic, rho_coupling)` for any concrete `Prox` impl.
/// The four MATLAB combinations are listed explicitly; Rust monomorphises
/// them per-prox-type.
fn run_2d_admm<P, Mu, Nu>(
    img: Image,
    gamma: f64,
    alpha: f64,
    prox: &P,
    mu_seq: Mu,
    nu_seq: Nu,
    isotropic: u8,
    rho_coupling: bool,
    tol: f64,
    max_iter: usize,
    verbose: bool,
) -> Image
where
    P: Prox,
    Mu: Fn(usize) -> f64,
    Nu: Fn(usize) -> f64,
{
    match (isotropic, rho_coupling) {
        (0, true) => admm_4connected_l2_ms(img, gamma, alpha, prox, mu_seq, nu_seq, tol, max_iter, verbose),
        (0, false) => admm_4connected_l2_ms(img, gamma, alpha, prox, mu_seq, no_rho_coupling, tol, max_iter, verbose),
        (1, true) => admm_8connected_l2_ms(img, gamma, alpha, prox, mu_seq, nu_seq, tol, max_iter, verbose),
        (1, false) => admm_8connected_l2_ms(img, gamma, alpha, prox, mu_seq, no_rho_coupling, tol, max_iter, verbose),
        (2, true) => admm_knight_l2_ms(img, gamma, alpha, prox, mu_seq, nu_seq, tol, max_iter, verbose),
        (2, false) => admm_knight_l2_ms(img, gamma, alpha, prox, mu_seq, no_rho_coupling, tol, max_iter, verbose),
        _ => unreachable!(),
    }
}

/// Solve the L2-Mumford-Shah within-segment smoothing problem
/// `min_μ Σ(y_i - μ_i)² + α Σ(μ_{i+1} - μ_i)²` with cached LU.
///
/// Builds a fresh `GaussL2Mum(n=len(y), alpha=alpha)`, calls `compute_mu`,
/// and returns the result as a 1-D numpy array. For repeated solves with
/// the same `n`, prefer the cached primitive (exposed in later phases).
#[pyfunction]
fn gauss_l2_mum_solve<'py>(
    py: Python<'py>,
    y: PyReadonlyArray1<f64>,
    alpha: f64,
) -> Bound<'py, PyArray1<f64>> {
    let y_slice = y.as_slice().unwrap();
    let n = y_slice.len().max(1);
    let g = GaussL2Mum::new(n, alpha);
    let mu = g.compute_mu(y_slice);
    mu.into_pyarray(py)
}

/// Compute the within-segment cost
/// `Σ(y_i - μ_i)² + α Σ(μ_{i+1} - μ_i)²` at the L2-MS optimum.
#[pyfunction]
fn gauss_l2_mum_cost(y: PyReadonlyArray1<f64>, alpha: f64) -> f64 {
    let y_slice = y.as_slice().unwrap();
    let n = y_slice.len().max(1);
    let g = GaussL2Mum::new(n, alpha);
    g.compute_dlr(y_slice)
}

/// Solve the 2-D L2 Mumford-Shah problem on a `(channels, rows, cols)` image
/// using the 4-connected (anisotropic) ADMM driver. Returns the smoothed
/// image with discontinuities, same shape as input.
///
/// Parameters mirror `mumfordShah2D.m`:
///   `gamma`       Potts (jump) penalty.
///   `alpha`       within-segment smoothness weight.
///   `tol`         relative-discrepancy stopping tolerance. Set 0 to force
///                 `max_iter` iterations (useful for parity testing).
///   `max_iter`    iteration cap.
///   `verbose`     print per-iteration diagnostics to stderr.
#[pyfunction]
#[pyo3(signature = (f, gamma, alpha, tol = 1e-3, max_iter = 50000, verbose = false, rho_coupling = true, isotropic = 0, weights = None, mu_schedule = None, nu_schedule = None, prox = None))]
fn min_l2_mum_2d<'py>(
    py: Python<'py>,
    f: PyReadonlyArray3<f64>,
    gamma: f64,
    alpha: f64,
    tol: f64,
    max_iter: usize,
    verbose: bool,
    rho_coupling: bool,
    isotropic: u8,
    weights: Option<PyReadonlyArray2<f64>>,
    mu_schedule: Option<PyReadonlyArray1<f64>>,
    nu_schedule: Option<PyReadonlyArray1<f64>>,
    prox: Option<Py<PyAny>>,
) -> PyResult<Bound<'py, PyArray3<f64>>> {
    if prox.is_some() && weights.is_some() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "pass either `prox` (an arbitrary callable) or `weights` (for the built-in L2 prox), not both",
        ));
    }
    let arr = f.as_array().to_owned();
    let img = Image::from_array(arr);
    let s_count: usize = match isotropic {
        0 => 2,
        1 => 4,
        2 => 8,
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unsupported isotropic mode {other}; expected 0, 1, or 2"
            )))
        }
    };

    // Closures over schedule data: when the user passes pre-computed arrays
    // we index them; otherwise fall back to the MATLAB defaults. Schedule
    // arrays are 1-indexed at the algorithm level (k = 1..=max_iter), so we
    // subtract 1 when indexing into the 0-indexed Vec. Bounds are clamped
    // to the last entry — caller may pass a shorter array for diagnostics.
    let mu_vec: Option<Vec<f64>> = mu_schedule.map(|s| s.as_slice().unwrap().to_vec());
    let nu_vec: Option<Vec<f64>> = nu_schedule.map(|s| s.as_slice().unwrap().to_vec());

    let mu_seq = move |k: usize| -> f64 {
        match &mu_vec {
            None => default_mu_seq(k),
            Some(v) => v[v.len().min(k).saturating_sub(1).max(0)],
        }
    };
    let nu_seq_default = move |k: usize| -> f64 {
        match &nu_vec {
            None => default_nu_seq(k, s_count),
            Some(v) => v[v.len().min(k).saturating_sub(1).max(0)],
        }
    };

    let result = match prox {
        Some(callable) => {
            let p = PyProx { callable };
            run_2d_admm(img, gamma, alpha, &p, mu_seq, nu_seq_default, isotropic, rho_coupling, tol, max_iter, verbose)
        }
        None => {
            let p = match weights {
                None => L2DataProx::new(img.clone()),
                Some(w) => L2DataProx::with_weights(img.clone(), w.as_array().to_owned()),
            };
            run_2d_admm(img, gamma, alpha, &p, mu_seq, nu_seq_default, isotropic, rho_coupling, tol, max_iter, verbose)
        }
    };
    Ok(result.data.into_pyarray(py))
}

#[pymodule]
#[pyo3(name = "_core")]
fn mumfordshah2d_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gauss_l2_mum_solve, m)?)?;
    m.add_function(wrap_pyfunction!(gauss_l2_mum_cost, m)?)?;
    m.add_function(wrap_pyfunction!(min_l2_mum_2d, m)?)?;
    Ok(())
}
