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

/// Weighted L2 data-fidelity prox. Mirrors `Auxiliary/makeProxL2w.m`:
///
///     prox(z, λ) = (w · f + λ · z) / (w + λ)
///
/// where `w` are per-pixel weights (broadcast across channels). Pixels
/// with `w = 0` are treated as unobserved (inpainting); the prox at those
/// pixels is `z`, completely driven by the consensus.
///
/// Construct with `L2DataProx::new(f)` for uniform `weights = 1` (the
/// default denoising case) or `L2DataProx::with_weights(f, weights)` for
/// the heteroscedastic / inpainting cases.
pub struct L2DataProx {
    pub f: Image,
    /// Per-pixel weights, shape `(rows, cols)`. `None` ≡ all ones.
    pub weights: Option<ndarray::Array2<f64>>,
}

impl L2DataProx {
    pub fn new(f: Image) -> Self {
        Self { f, weights: None }
    }
    pub fn with_weights(f: Image, weights: ndarray::Array2<f64>) -> Self {
        debug_assert_eq!(
            weights.shape(),
            &[f.n_row(), f.n_col()],
            "weights shape must match (rows, cols)"
        );
        Self { f, weights: Some(weights) }
    }
}

impl Prox for L2DataProx {
    fn apply(&self, z: &Image, lambda: f64, out: &mut Image) {
        assert_eq!(out.data.shape(), self.f.data.shape());
        match &self.weights {
            None => {
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
            Some(w) => {
                let (channels, n_row, n_col) = (out.channels(), out.n_row(), out.n_col());
                for c in 0..channels {
                    for i in 0..n_row {
                        for j in 0..n_col {
                            let wv = w[[i, j]];
                            let denom = wv + lambda;
                            out.data[[c, i, j]] = (wv * self.f.data[[c, i, j]]
                                + lambda * z.data[[c, i, j]])
                                / denom;
                        }
                    }
                }
            }
        }
    }
}

/// Generalised consensus ADMM for the 2-D L2 Mumford-Shah problem.
///
/// `nhood` lists the primary direction stencils (each `[dx, dy]`) and their
/// per-direction ω weight. Each entry contributes 2 directions to the
/// ADMM (un-rotated and 90°-rotated), so `S = 2 · nhood.len()`. Common
/// choices match `mumfordShah2D.m`:
///
///   isotropic = 0 (4-connected, anisotropic): &[([1,0], 1.0)]                S = 2
///   isotropic = 1 (8-connected, near-isotropic): &[([1,0], √2-1), ([1,1], 1-√2/2)]    S = 4
///   isotropic = 2 (knight-move): &[([1,0], …), ([1,1], …), ([2,1], …), ([1,2], …)]   S = 8
///
/// Optional inter-direction ρ coupling: when `nu_seq(k) ≠ 0`, ρ dual
/// variables for each `r < t` pair (`C(S, 2)` of them) strengthen
/// convergence — this matches MATLAB's default `rhoFlag=true`. Pass
/// `nu_seq = |_| 0.0` to disable.
pub fn admm_l2_ms<P, MuSeq, NuSeq>(
    f_data: Image,
    gamma: f64,
    alpha: f64,
    prox: &P,
    nhood: &[([i64; 2], f64)],
    mu_seq: MuSeq,
    nu_seq: NuSeq,
    tol: f64,
    max_iter: usize,
    verbose: bool,
) -> Image
where
    P: Prox,
    MuSeq: Fn(usize) -> f64,
    NuSeq: Fn(usize) -> f64,
{
    assert!(!nhood.is_empty(), "nhood must contain at least one direction");
    let (channels, n_row, n_col) = (f_data.channels(), f_data.n_row(), f_data.n_col());
    let s_count: usize = 2 * nhood.len();

    // u_s for s in 0..S; ρ for the C(S, 2) directed-pair (r, t) with r < t.
    // MATLAB initialises u{s} = backproj = prox(0, mu(1)) ≈ f (rather than 0)
    // so that the iter-1 cross-direction `nu * u{t}` term sees a sensible
    // starting point. Mirroring that here.
    let mut u: Vec<Image> = (0..s_count).map(|_| f_data.clone()).collect();
    {
        // Apply prox(zero_image, mu(1)) once to seed u_s exactly the same way.
        let mu1 = mu_seq(1);
        let zero = Image::zeros(channels, n_row, n_col);
        let mut backproj = Image::zeros(channels, n_row, n_col);
        prox.apply(&zero, mu1, &mut backproj);
        for s in 0..s_count {
            u[s] = backproj.clone();
        }
    }
    let mut lam: Vec<Image> = (0..s_count)
        .map(|_| Image::zeros(channels, n_row, n_col))
        .collect();
    // rho_pair[(r, t)] for r < t. For S=2: one pair (0, 1). For S=4: six.
    let mut rho_pair: Vec<Image> = (0..s_count * (s_count - 1) / 2)
        .map(|_| Image::zeros(channels, n_row, n_col))
        .collect();
    let mut v = f_data.clone();
    let mut z = Image::zeros(channels, n_row, n_col);

    let pair_index = |r: usize, t: usize| -> usize {
        // Triangular indexing for r < t in 0..s_count: 0 for (0,1), then
        // (0,2), (1,2), (0,3), (1,3), (2,3), … (column-major upper triangle).
        debug_assert!(r < t);
        t * (t - 1) / 2 + r
    };

    let rho_flag = nu_seq(1) != 0.0;

    for k in 1..=max_iter {
        let mu = mu_seq(k);
        let nu = if rho_flag { nu_seq(k) } else { 0.0 };
        let denom_w = mu + nu * (s_count as f64 - 1.0);
        let scale = 2.0 / denom_w;

        z.data.fill(0.0);

        for s in 0..s_count {
            let nhood_idx = s / 2;
            let (direction, omega) = nhood[nhood_idx];
            let gamma_p = scale * omega * gamma;
            let alpha_p = scale * omega * alpha;
            // w = mu * v + lams[s]  (then optionally + Σ_{r<s} ν u[r] - ρ[r,s]
            //                                    + Σ_{t>s} ν u[t] + ρ[s,t])
            //   then divide by (mu + ν (S-1)).
            let mut w = Image::zeros(channels, n_row, n_col);
            for ((slot, vv), lv) in w.data.iter_mut().zip(v.data.iter()).zip(lam[s].data.iter()) {
                *slot = mu * vv + lv;
            }
            // Sign convention (verbatim from `mumfordShah2D.m:139–146`):
            //   r < s  →  w += ν u_r + ρ_{r,s}
            //   t > s  →  w += ν u_t − ρ_{s,t}
            if rho_flag {
                for r in 0..s {
                    let rho = &rho_pair[pair_index(r, s)];
                    for ((slot, ur), rh) in
                        w.data.iter_mut().zip(u[r].data.iter()).zip(rho.data.iter())
                    {
                        *slot += nu * ur + rh;
                    }
                }
                for t in (s + 1)..s_count {
                    let rho = &rho_pair[pair_index(s, t)];
                    for ((slot, ut), rh) in
                        w.data.iter_mut().zip(u[t].data.iter()).zip(rho.data.iter())
                    {
                        *slot += nu * ut - rh;
                    }
                }
            }
            for slot in w.data.iter_mut() {
                *slot /= denom_w;
            }

            // MATLAB rotates the *odd*-indexed (1-indexed) directions (s=1, 3, …),
            // which in 0-indexed Rust are s=0, 2, …. Rotation direction must
            // match MATLAB's `rotate90(w, 1)` (CCW before, CW back).
            if s % 2 == 0 {
                let mut w_rot = w.rotate90_ccw();
                run_direction(&mut w_rot, direction, gamma_p, alpha_p);
                u[s] = w_rot.rotate90_cw();
            } else {
                run_direction(&mut w, direction, gamma_p, alpha_p);
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

        // ρ_{r, t} = ρ_{r, t} + ν (u_r − u_t)
        if rho_flag {
            for r in 0..s_count {
                for t in (r + 1)..s_count {
                    let idx = pair_index(r, t);
                    // Borrow u[r] and u[t] immutably; rho_pair[idx] mutably.
                    let (ur, ut) = (&u[r], &u[t]);
                    let rho = &mut rho_pair[idx];
                    for ((slot, urv), utv) in rho.data.iter_mut().zip(ur.data.iter()).zip(ut.data.iter()) {
                        *slot += nu * (urv - utv);
                    }
                }
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
                "iter {k}: mu={mu:.3e}, nu={nu:.3e}, err={err:.3e}, S={s_count}"
            );
        }

        if err < tol {
            return v;
        }
    }

    v
}

/// Default ν schedule from `mumfordShah2D.m` for `S` directions:
///   ν(k) = μ(k) · S / nchoosek(S, 2)
/// For S=2: ν(k) = 2 μ(k); for S=4: ν(k) = (4/6) μ(k) = (2/3) μ(k).
pub fn default_nu_seq(k: usize, s_count: usize) -> f64 {
    let nck = (s_count * (s_count - 1) / 2) as f64;
    default_mu_seq(k) * s_count as f64 / nck
}

/// `nu_seq` constructor that returns 0.0 always — equivalent to passing
/// `'nuSeq', @(k) 0` in MATLAB. Disables ρ coupling.
pub fn no_rho_coupling(_k: usize) -> f64 {
    0.0
}

/// 4-connected (anisotropic) neighbourhood — single direction `[1, 0]`,
/// ω = 1. Maps to `mumfordShah2D.m` `'isotropic', 0`. S = 2.
pub const NHOOD_4_CONNECTED: &[([i64; 2], f64)] = &[([1, 0], 1.0)];

/// 8-connected (near-isotropic) neighbourhood — `[1,0]` and `[1,1]` with
/// the MATLAB-default isotropy weights. Maps to `'isotropic', 1`. S = 4.
/// `omega_1 = √2 − 1 ≈ 0.4142`, `omega_2 = 1 − √2/2 ≈ 0.2929`.
pub fn nhood_8_connected() -> [([i64; 2], f64); 2] {
    [
        ([1, 0], std::f64::consts::SQRT_2 - 1.0),
        ([1, 1], 1.0 - std::f64::consts::SQRT_2 / 2.0),
    ]
}

/// Knight-move neighbourhood — `[1,0]`, `[1,1]`, `[2,1]`, `[1,2]` with
/// the MATLAB-default near-isotropic weights. Maps to `'isotropic', 2`. S = 8.
///   ω_1 = √5 − 2 ≈ 0.2361
///   ω_2 = √5 − (3/2)√2 ≈ 0.1146
///   ω_3 = ω_4 = ½(1 + √2 − √5) ≈ 0.0827
pub fn nhood_knight_move() -> [([i64; 2], f64); 4] {
    let s5 = 5.0_f64.sqrt();
    let s2 = std::f64::consts::SQRT_2;
    let omega_3 = 0.5 * (1.0 + s2 - s5);
    [
        ([1, 0], s5 - 2.0),
        ([1, 1], s5 - 1.5 * s2),
        ([2, 1], omega_3),
        ([1, 2], omega_3),
    ]
}

/// Backward-compatible wrapper: 4-connected ADMM (S=2). Equivalent to
/// `admm_l2_ms(..., NHOOD_4_CONNECTED, ...)`.
pub fn admm_4connected_l2_ms<P, MuSeq, NuSeq>(
    f_data: Image,
    gamma: f64,
    alpha: f64,
    prox: &P,
    mu_seq: MuSeq,
    nu_seq: NuSeq,
    tol: f64,
    max_iter: usize,
    verbose: bool,
) -> Image
where
    P: Prox,
    MuSeq: Fn(usize) -> f64,
    NuSeq: Fn(usize) -> f64,
{
    admm_l2_ms(
        f_data,
        gamma,
        alpha,
        prox,
        NHOOD_4_CONNECTED,
        mu_seq,
        nu_seq,
        tol,
        max_iter,
        verbose,
    )
}

/// 8-connected (near-isotropic) ADMM (S=4). Maps to MATLAB's
/// `mumfordShah2D(..., 'isotropic', 1)`.
pub fn admm_8connected_l2_ms<P, MuSeq, NuSeq>(
    f_data: Image,
    gamma: f64,
    alpha: f64,
    prox: &P,
    mu_seq: MuSeq,
    nu_seq: NuSeq,
    tol: f64,
    max_iter: usize,
    verbose: bool,
) -> Image
where
    P: Prox,
    MuSeq: Fn(usize) -> f64,
    NuSeq: Fn(usize) -> f64,
{
    let nhood = nhood_8_connected();
    admm_l2_ms(
        f_data, gamma, alpha, prox, &nhood, mu_seq, nu_seq, tol, max_iter, verbose,
    )
}

/// Knight-move ADMM (S=8). Maps to MATLAB's `mumfordShah2D(..., 'isotropic', 2)`
/// — the closest-to-isotropic discretisation in the family.
pub fn admm_knight_l2_ms<P, MuSeq, NuSeq>(
    f_data: Image,
    gamma: f64,
    alpha: f64,
    prox: &P,
    mu_seq: MuSeq,
    nu_seq: NuSeq,
    tol: f64,
    max_iter: usize,
    verbose: bool,
) -> Image
where
    P: Prox,
    MuSeq: Fn(usize) -> f64,
    NuSeq: Fn(usize) -> f64,
{
    let nhood = nhood_knight_move();
    admm_l2_ms(
        f_data, gamma, alpha, prox, &nhood, mu_seq, nu_seq, tol, max_iter, verbose,
    )
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
        let prox = L2DataProx::new(f.clone());
        let result = admm_4connected_l2_ms(
            f.clone(),
            1.0,
            1.0,
            &prox,
            default_mu_seq,
            no_rho_coupling,
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
        let prox = L2DataProx::new(f.clone());
        let result = admm_4connected_l2_ms(
            f.clone(),
            0.05,
            1e-6,
            &prox,
            default_mu_seq,
            no_rho_coupling,
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
