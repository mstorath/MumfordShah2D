// Ported from Java/src/mumfordShah/ (Hohm, Storath & Weinmann)
// by Claude Sonnet coding agent, Anthropic, 2026.
//
// PyO3 module: exposes Phase 1 Rust primitives to Python.

mod image;
mod gauss_elim;
mod mumfordshah_1d;
mod direction_processor;
mod admm;

use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;

use crate::gauss_elim::{GaussElim, GaussL2Mum};

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
    mu.into_pyarray_bound(py)
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

#[pymodule]
#[pyo3(name = "_core")]
fn mumfordshah2d_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gauss_l2_mum_solve, m)?)?;
    m.add_function(wrap_pyfunction!(gauss_l2_mum_cost, m)?)?;
    Ok(())
}
