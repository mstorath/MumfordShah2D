//! 1-D L2 Mumford-Shah solver: dynamic programming over Potts partitions
//! with within-segment smoothing via the cached tridiagonal LU
//! (`GaussL2Mum`).
//!
//! Solves
//!
//!     min_{u, J}  γ |J|  +  α Σ_{i ∉ J} (u_{i+1} - u_i)²  +  Σ_i ||u_i - y_i||²
//!
//! where J is the set of jump indices (Potts cardinality penalty) and the
//! within-segment smoothing energy is the L2-MS objective solved exactly
//! by the cached `GaussL2Mum` for any segment length.
//!
//! Direct port of `Java/src/mumfordShah/{AbstractPotts1D,L2MumfordShah1DFaster}.java`.
//! The cost-recurrence state machine in `compute_dlr` is preserved one-for-one;
//! the AbstractPotts1D DP loop drives it in the call order
//!
//!     for r in 1..=n:
//!         compute_dlr(1, r)
//!         for l in (2..=r).rev():
//!             compute_dlr(l, r)   // exits early via the gamma-bounded break
//!
//! which the four-state transition logic in `compute_dlr` exploits to
//! compute h(l, r) in O(1) per step using scriptF.

use crate::gauss_elim::{GaussElim, GaussL2Mum};

/// 1-D L2 Mumford-Shah solver. Holds owned input `f` and a borrowed
/// pre-built `GaussL2Mum` of capacity ≥ `f[0].len()`. Borrowing the cached
/// LU avoids re-factorisation across many stripe solves in the 2-D driver.
pub struct L2MumfordShah1D<'g> {
    n: usize,
    channels: usize,
    gamma: f64,
    f: Vec<Vec<f64>>, // f[channel][index]
    gauss: &'g GaussL2Mum,
    script_f: Vec<f64>,

    // recurrence state — preserved across compute_dlr() calls in the
    // specific order used by the AbstractPotts1D DP loop.
    cur_l: usize,
    cur_r: usize,
    prev_d1r: Vec<f64>,
    d1r: Vec<f64>,
    prev_u1r: Vec<f64>,
    u1r: Vec<f64>,
    prev_dlr: Vec<f64>,
    dlr: Vec<f64>,
    prev_ulr: Vec<f64>,
    ulr: Vec<f64>,
    temp_f: Vec<Vec<f64>>,

    // DP results
    partition: Vec<usize>, // partition[r] = best left-bound (0 if r belongs to first segment)
    result: Vec<Vec<f64>>,
}

impl<'g> L2MumfordShah1D<'g> {
    /// `f` has shape `(channels, n)`. `gauss` must have capacity ≥ n.
    pub fn new(f: Vec<Vec<f64>>, gamma: f64, alpha: f64, gauss: &'g GaussL2Mum) -> Self {
        let channels = f.len();
        assert!(channels > 0, "f must have at least one channel");
        let n = f[0].len();
        for ch in &f {
            assert_eq!(ch.len(), n, "all channels must have the same length");
        }
        assert!(gauss.n >= n, "GaussL2Mum capacity {} < n {}", gauss.n, n);

        // scriptF[k] for k=0..n-1: scriptF[0] = 1; scriptF[k] = α scriptF[k-1] / (α + scriptF[k-1]) + 1
        let mut script_f = vec![0.0; n];
        if n > 0 {
            script_f[0] = 1.0;
            for k in 1..n {
                let prev = script_f[k - 1];
                script_f[k] = alpha * prev / (alpha + prev) + 1.0;
            }
        }

        let result = vec![vec![0.0; n]; channels];

        Self {
            n,
            channels,
            gamma,
            f,
            gauss,
            script_f,
            cur_l: 0,
            cur_r: 0,
            prev_d1r: vec![],
            d1r: vec![],
            prev_u1r: vec![],
            u1r: vec![],
            prev_dlr: vec![],
            dlr: vec![],
            prev_ulr: vec![],
            ulr: vec![],
            temp_f: vec![],
            partition: vec![0; n + 1],
            result,
        }
    }

    /// Run the DP partition search, then reconstruct `u` per segment.
    /// Returns `(channels, n)` Mumford-Shah optimal output.
    pub fn solve(mut self) -> Vec<Vec<f64>> {
        self.find_best_partition();
        self.segmentation_from_partition();
        self.result
    }

    /// Forward DP: ports `AbstractPotts1D.findBestPartition`. Computes
    /// `partition[r]` = best left-bound (1-indexed origin in Java terms;
    /// stored 0-indexed here so `partition[r] = 0` means the first segment
    /// extends to `r`).
    fn find_best_partition(&mut self) {
        // B[r] is the optimal cost using indices [1..=r] (1-indexed in Java).
        // Initialise B[0] = -gamma so that B[r] = γ + h(1, r) for a single segment.
        let mut b_arr = vec![0.0; self.n + 1];
        b_arr[0] = -self.gamma;

        for r in 1..=self.n {
            b_arr[r] = self.compute_dlr(1, r);

            for l in (2..=r).rev() {
                let dlr = self.compute_dlr(l, r);
                // Acceleration: cost h(l, r) is non-decreasing as l shrinks (more data points
                // accumulated into the segment). When B[r] < γ + h(l, r), no smaller-l
                // candidate can improve B[r], so we stop.
                if b_arr[r] < self.gamma + dlr {
                    break;
                }
                let b = b_arr[l - 1] + self.gamma + dlr;
                if b <= b_arr[r] {
                    b_arr[r] = b;
                    self.partition[r] = l - 1;
                }
            }
        }
    }

    /// Backward reconstruction: ports `AbstractPotts1D.segmentationFromPartition`.
    fn segmentation_from_partition(&mut self) {
        let mut r = self.n;
        let mut l = self.partition[r];
        while r > 0 {
            // Java: computeMu_LtoR(l+1, r) takes 1-indexed bounds; this Rust port
            // takes 0-indexed half-open [l, r) so we map directly.
            let mu = self.compute_mu_ltr(l, r);
            for ch in 0..self.channels {
                for t in 0..(r - l) {
                    self.result[ch][l + t] = mu[ch][t];
                }
            }
            r = l;
            l = self.partition[r];
        }
    }

    /// Solve the within-segment L2-MS on indices `[l, r)` via the cached LU.
    /// Returns shape `(channels, r - l)`.
    fn compute_mu_ltr(&self, l: usize, r: usize) -> Vec<Vec<f64>> {
        let len = r - l;
        let mut out = vec![vec![0.0; len]; self.channels];
        let mut buf = vec![0.0; len];
        for ch in 0..self.channels {
            buf.copy_from_slice(&self.f[ch][l..r]);
            let mu = self.gauss.compute_mu(&buf);
            out[ch].copy_from_slice(&mu);
        }
        out
    }

    /// Stateful within-segment L2-MS cost recurrence. Mirrors
    /// `L2MumfordShah1DFaster.computeDlr` exactly.
    ///
    /// Both `l` and `r` are 1-indexed (matching Java's convention) for direct
    /// transcription; index arithmetic on `f` and `temp_f` uses `r-l` and
    /// `r-1` to recover 0-indexed positions.
    fn compute_dlr(&mut self, l: usize, r: usize) -> f64 {
        if r == 1 && l == 1 {
            // Initial state: starting the [1, r] series.
            self.prev_d1r = vec![0.0; self.channels];
            self.d1r = vec![0.0; self.channels];
            self.prev_u1r = vec![0.0; self.channels];
            self.u1r = vec![0.0; self.channels];
            for ch in 0..self.channels {
                self.u1r[ch] = self.f[ch][r - 1];
            }
            self.cur_r = r;
            self.cur_l = l;
        } else if r == self.cur_r + 1 && l == 1 {
            // Advance r in the [1, r] series.
            self.prev_d1r = std::mem::take(&mut self.d1r);
            self.d1r = vec![0.0; self.channels];
            self.prev_u1r = std::mem::take(&mut self.u1r);
            self.u1r = vec![0.0; self.channels];
            self.cur_r = r;
            self.cur_l = l;
        } else if r == self.cur_r && l == self.cur_r {
            // Start the backwards [l, r] recurrence at the current right bound.
            self.prev_dlr = vec![0.0; self.channels];
            self.dlr = vec![0.0; self.channels];
            self.prev_ulr = vec![0.0; self.channels];
            self.ulr = vec![0.0; self.channels];
            // Reverse f[0..r-1] into temp_f (Java prepares only positions
            // 0..r-1 because the recurrence walks from r downward and uses
            // tempF[r-l] as the "next" data point).
            let len = r - 1;
            self.temp_f = vec![vec![0.0; len]; self.channels];
            for ch in 0..self.channels {
                for iter in 0..len {
                    self.temp_f[ch][iter] = self.f[ch][r - 1 - iter];
                }
            }
            for ch in 0..self.channels {
                self.ulr[ch] = self.temp_f[ch][0];
            }
            self.cur_r = r;
            self.cur_l = l;
        } else if r == self.cur_r && l == self.cur_l - 1 {
            // Extend [l, r] leftward by one.
            self.prev_dlr = std::mem::take(&mut self.dlr);
            self.dlr = vec![0.0; self.channels];
            self.prev_ulr = std::mem::take(&mut self.ulr);
            self.ulr = vec![0.0; self.channels];
            self.cur_r = r;
            self.cur_l = l;
        } else {
            panic!(
                "compute_dlr called with unexpected (l={}, r={}); cur=({}, {})",
                l, r, self.cur_l, self.cur_r
            );
        }

        // Recurrence body. r-l ≥ 0 here.
        let k = r - l;
        let sf = self.script_f[k];

        if l == 1 {
            for ch in 0..self.channels {
                let prev_u = self.prev_u1r[ch];
                let f_val = self.f[ch][r - l]; // Java: f[dim][r-l] — 0-indexed access at position r-l
                let new_u = ((sf - 1.0) * prev_u + f_val) / sf;
                self.u1r[ch] = new_u;
                let du = new_u - prev_u;
                let df = new_u - f_val;
                self.d1r[ch] = self.prev_d1r[ch] + (sf - 1.0) * du * du + df * df;
            }
            self.d1r.iter().sum()
        } else {
            for ch in 0..self.channels {
                let prev_u = self.prev_ulr[ch];
                let f_val = self.temp_f[ch][r - l];
                let new_u = ((sf - 1.0) * prev_u + f_val) / sf;
                self.ulr[ch] = new_u;
                let du = new_u - prev_u;
                let df = new_u - f_val;
                self.dlr[ch] = self.prev_dlr[ch] + (sf - 1.0) * du * du + df * df;
            }
            self.dlr.iter().sum()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn solve_scalar(f: &[f64], gamma: f64, alpha: f64) -> Vec<f64> {
        let n = f.len();
        let gauss = GaussL2Mum::new(n, alpha);
        let solver = L2MumfordShah1D::new(vec![f.to_vec()], gamma, alpha, &gauss);
        let out = solver.solve();
        out.into_iter().next().unwrap()
    }

    #[test]
    fn constant_signal_unchanged_at_any_gamma() {
        let f = vec![3.14; 8];
        for &gamma in &[0.01, 1.0, 100.0] {
            let u = solve_scalar(&f, gamma, 1.0);
            for &v in &u {
                assert!((v - 3.14).abs() < 1e-12, "gamma={gamma}: got {v}");
            }
        }
    }

    #[test]
    fn sharp_step_split_at_small_gamma_no_smoothing_at_alpha_zero() {
        // alpha = 0 means within-segment cost is just sum of (y - mu)^2,
        // optimal mu per segment is the mean. Step signal with γ small enough
        // to admit a jump should split into two flat segments at the mean of
        // each half.
        let f = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let u = solve_scalar(&f, 0.01, 0.0);
        assert!(
            u.iter().take(4).all(|&v| v.abs() < 1e-12),
            "left half: {u:?}"
        );
        assert!(
            u.iter().skip(4).all(|&v| (v - 1.0).abs() < 1e-12),
            "right half: {u:?}"
        );
    }

    #[test]
    fn large_gamma_collapses_to_single_segment() {
        let f = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        // γ=1e6 prevents any jump; α=1.0 smooths the single segment.
        let u = solve_scalar(&f, 1e6, 1.0);
        // Reference: same input through GaussL2Mum directly with α=1.0.
        let gauss = GaussL2Mum::new(f.len(), 1.0);
        let mu = gauss.compute_mu(&f);
        for (a, b) in u.iter().zip(mu.iter()) {
            assert!((a - b).abs() < 1e-12, "got {a}, expected {b}");
        }
    }

    #[test]
    fn alpha_smoothing_at_inf_gamma() {
        // gamma = ∞, alpha > 0: a single segment of GaussL2Mum smoothing.
        let f = vec![0.0, 1.0, 0.0, 1.0, 0.0];
        let u = solve_scalar(&f, 1e9, 1.0);
        // Reference: GaussL2Mum directly solves the same problem.
        let gauss = GaussL2Mum::new(f.len(), 1.0);
        let mu = gauss.compute_mu(&f);
        for (a, b) in u.iter().zip(mu.iter()) {
            assert!((a - b).abs() < 1e-12, "got {a}, expected {b}");
        }
    }

    #[test]
    fn multichannel_keeps_segment_structure() {
        // Two-channel input where both channels have a step at the same place.
        let f0 = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let f1 = vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0];
        let n = f0.len();
        let gauss = GaussL2Mum::new(n, 0.0);
        let solver = L2MumfordShah1D::new(vec![f0.clone(), f1.clone()], 0.01, 0.0, &gauss);
        let out = solver.solve();
        // Each channel should split at the same place.
        for ch in 0..2 {
            let exp = if ch == 0 { (0.0, 1.0) } else { (0.0, 2.0) };
            for i in 0..3 {
                assert!((out[ch][i] - exp.0).abs() < 1e-12);
            }
            for i in 3..6 {
                assert!((out[ch][i] - exp.1).abs() < 1e-12);
            }
        }
    }
}
