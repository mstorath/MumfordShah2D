//! 2-D L2 Mumford-Shah ADMM driver. Direct port of `mumfordShah2D.m`'s
//! consensus splitting strategy with the L2 within-direction inner solve.
//!
//! Algorithm (S = 2 for the 4-connected case, isotropic=0):
//!   z = 0
//!   for s in 1..=S:
//!       w_s    = (μ v + λ_s) / μ                  (no rho coupling for nu=0)
//!       γ', α' = (2/μ) ω_s · {γ, α}
//!       if s is even:
//!           u_s = rotate90⁻¹( DirectionProcessor( rotate90(w_s), [1,0], γ', α' ) )
//!       else:
//!           u_s = DirectionProcessor( w_s, [1,0], γ', α' )
//!       z += u_s − λ_s/μ
//!   z /= S
//!   v = prox(z, μS)
//!   λ_s = λ_s + μ(v − u_s)
//!
//! Stopping criterion (verbatim from `mumfordShah2D.m:192`):
//!     err = ‖u₁ − u₂‖₂ / (‖u₁‖₂ + ‖u₂‖₂) < tol
//!
//! Default `prox`: L2 data fidelity, `prox(z, λ) = (f + λ z) / (1 + λ)`.
//! Other variants (Tikhonov via linear operator A, inpainting, etc.) can
//! be plugged in via `Prox` later — kept minimal for v0.1.

use crate::direction_processor::run_direction;
use crate::image::{directions, Image};

/// Trait for the data-fidelity proximal operator.
/// `prox(z, lambda)` solves
///     min_v  data_fidelity(v)  +  (lambda / 2) ||v - z||²
/// in-place on `out`.
pub trait Prox {
    fn apply(&self, z: &Image, lambda: f64, out: &mut Image);
}

/// Default L2 data fidelity: `min_v ||v - f||² + (lambda/2)||v - z||²`
/// → `v = (f + (lambda/2) z) / (1 + lambda/2)`. With the convention used
/// by `mumfordShah2D.m` where the data term is `||u - f||²` (no half),
/// the prox simplifies to `v = (2 f + lambda z) / (2 + lambda)`. We use
/// the half-factor convention to match Pottslab's prox style.
pub struct L2DataProx {
    pub f: Image,
}

impl Prox for L2DataProx {
    fn apply(&self, z: &Image, lambda: f64, out: &mut Image) {
        assert_eq!(out.data.shape(), self.f.data.shape());
        let denom = 1.0 + lambda;
        for ((slot, fv), zv) in out
            .data
            .iter_mut()
            .zip(self.f.data.iter())
            .zip(z.data.iter())
        {
            *slot = (fv + lambda * zv) / denom;
        }
    }
}

/// 4-connected (anisotropic) ADMM for the 2-D L2 Mumford-Shah problem.
/// Two directions (S=2): rows + columns. Single ω = 1, no inter-direction
/// rho coupling.
///
/// `mu_seq[k]` is the dual step size at iteration `k` (1-indexed, like the
/// MATLAB driver). Default in MATLAB: `mu_seq(k) = k².⁰¹ · 10⁻⁶`. The
/// caller supplies the sequence (closure) so we don't bake in the schedule.
pub fn admm_4connected_l2_ms<P: Prox>(
    f_data: Image,
    gamma: f64,
    alpha: f64,
    prox: &P,
    mu_seq: impl Fn(usize) -> f64,
    tol: f64,
    max_iter: usize,
    verbose: bool,
) -> Image {
    let (channels, n_row, n_col) = (f_data.channels(), f_data.n_row(), f_data.n_col());
    let s_count: usize = 2; // 4-conn → S=2

    // u_1 = vertical (rows-as-stripes via rotation), u_2 = horizontal — but the
    // MATLAB driver enumerates them the other way: s=1 (odd) processes
    // direction [1,0] = column-stripes, s=2 (even) rotates first → row-stripes.
    let mut u: Vec<Image> = (0..s_count)
        .map(|_| Image::zeros(channels, n_row, n_col))
        .collect();
    let mut lam: Vec<Image> = (0..s_count)
        .map(|_| Image::zeros(channels, n_row, n_col))
        .collect();
    let mut v = f_data.clone();
    let mut z = Image::zeros(channels, n_row, n_col);

    for k in 1..=max_iter {
        let mu = mu_seq(k);
        let denom_w = mu; // no rho coupling for S=2 / nu=0 → denominator is mu
        let scale = 2.0 / denom_w;
        let gamma_p = scale * gamma;
        let alpha_p = scale * alpha;

        // Build z by accumulating u_s − λ_s/μ.
        z.data.fill(0.0);

        for s in 0..s_count {
            // w = (μ v + λ_s) / μ = v + λ_s / μ
            let mut w = Image::zeros(channels, n_row, n_col);
            for ((slot, vv), lv) in w.data.iter_mut().zip(v.data.iter()).zip(lam[s].data.iter()) {
                *slot = vv + lv / mu;
            }

            // Even s: rotate, run, rotate back.
            if s % 2 == 1 {
                let mut w_rot = w.rotate90_cw();
                run_direction(&mut w_rot, directions::ROW_STEP, gamma_p, alpha_p);
                u[s] = w_rot.rotate90_ccw();
            } else {
                run_direction(&mut w, directions::ROW_STEP, gamma_p, alpha_p);
                u[s] = w;
            }

            // z += u_s − λ_s / μ
            for ((slot, uv), lv) in z
                .data
                .iter_mut()
                .zip(u[s].data.iter())
                .zip(lam[s].data.iter())
            {
                *slot += uv - lv / mu;
            }
        }

        // z /= S
        for v_ in z.data.iter_mut() {
            *v_ /= s_count as f64;
        }

        // v = prox(z, μ S)
        prox.apply(&z, mu * s_count as f64, &mut v);

        // λ_s = λ_s + μ (v − u_s)
        for s in 0..s_count {
            for ((slot, vv), uv) in lam[s]
                .data
                .iter_mut()
                .zip(v.data.iter())
                .zip(u[s].data.iter())
            {
                *slot += mu * (vv - uv);
            }
        }

        // Stopping criterion: relative discrepancy between the two u_s.
        let mut diff_sq = 0.0;
        let mut u0_sq = 0.0;
        let mut u1_sq = 0.0;
        for ((u0, u1), _) in u[0].data.iter().zip(u[1].data.iter()).zip(0..) {
            diff_sq += (u0 - u1) * (u0 - u1);
            u0_sq += u0 * u0;
            u1_sq += u1 * u1;
        }
        let denom = u0_sq.sqrt() + u1_sq.sqrt();
        let err = if denom > 0.0 {
            diff_sq.sqrt() / denom
        } else {
            f64::INFINITY
        };

        if verbose {
            eprintln!(
                "iter {k}: mu={mu:.3e}, err={err:.3e}, gamma'={gamma_p:.3e}, alpha'={alpha_p:.3e}"
            );
        }

        if err < tol {
            return v;
        }
    }

    v
}

/// Default mu schedule from `mumfordShah2D.m`:  μ(k) = k².⁰¹ · 10⁻⁶.
pub fn default_mu_seq(k: usize) -> f64 {
    (k as f64).powf(2.01) * 1e-6
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array3;

    fn make_image(channels: usize, n_row: usize, n_col: usize, fill: impl Fn(usize, usize, usize) -> f64) -> Image {
        let mut a = Array3::<f64>::zeros((channels, n_row, n_col));
        for c in 0..channels {
            for i in 0..n_row {
                for j in 0..n_col {
                    a[[c, i, j]] = fill(c, i, j);
                }
            }
        }
        Image::from_array(a)
    }

    #[test]
    fn constant_image_converges_to_itself() {
        let f = make_image(1, 4, 4, |_, _, _| 7.0);
        let prox = L2DataProx { f: f.clone() };
        let result = admm_4connected_l2_ms(
            f.clone(),
            1.0,
            1.0,
            &prox,
            default_mu_seq,
            1e-3,
            500,
            false,
        );
        for &v in result.data.iter() {
            assert!((v - 7.0).abs() < 1e-3, "got {v}, expected 7.0");
        }
    }

    #[test]
    fn binary_step_image_recovers_two_segments_at_small_gamma() {
        // 4×4 image, left half 0, right half 1, no noise.
        let f = make_image(1, 4, 4, |_, _, j| if j < 2 { 0.0 } else { 1.0 });
        let prox = L2DataProx { f: f.clone() };
        let result = admm_4connected_l2_ms(
            f.clone(),
            0.05,
            1e-6,
            &prox,
            default_mu_seq,
            1e-3,
            2000,
            false,
        );
        // Left half should be near 0, right half near 1.
        for i in 0..4 {
            for j in 0..2 {
                assert!(
                    result.data[[0, i, j]].abs() < 0.05,
                    "left ({i},{j}) = {}",
                    result.data[[0, i, j]]
                );
            }
            for j in 2..4 {
                assert!(
                    (result.data[[0, i, j]] - 1.0).abs() < 0.05,
                    "right ({i},{j}) = {}",
                    result.data[[0, i, j]]
                );
            }
        }
    }
}
