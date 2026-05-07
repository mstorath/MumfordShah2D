// Ported from Java/src/mumfordShah/GaussElim.java and GaussL2Mum.java
// (Hohm, Storath & Weinmann) by Claude Sonnet coding agent, Anthropic, 2026.
//
// Solves the within-segment L2-Mumford-Shah smoothing problem
//
//     min_μ  Σ_i (y_i - μ_i)²  +  α Σ_i (μ_{i+1} - μ_i)²
//
// via a pre-computed tridiagonal LU factorisation of the matrix
//
//     [1+α   -α                       ]
//     [-α    1+2α  -α                 ]
//     [      -α    1+2α  -α           ]
//     [            ...                ]
//     [                  -α    1+α    ]
//
// The Java implementation caches LU factors for every length 1..n at
// construction so that every per-segment solve later is O(len). We mirror
// that caching strategy.
//
// Index convention: matches the Java source one-for-one. The diagonal
// vector for length-`k` systems is `diagonals[k]` of length `k+1` and the
// off-diagonal multiplier vector is `factors[k]` of length `k`.

/// Trait for the abstract `GaussElim` base class. The Java `computeDlr`
/// helper is provided for free against any concrete implementation.
pub trait GaussElim {
    fn alpha(&self) -> f64;
    fn compute_mu(&self, y: &[f64]) -> Vec<f64>;

    /// Within-segment L2-Mumford-Shah cost
    /// `Σ_i (y_i - μ_i)² + α Σ_i (μ_{i+1} - μ_i)²`. Direct port of
    /// `GaussElim.computeDlr` in the Java source.
    fn compute_dlr(&self, y: &[f64]) -> f64 {
        let mu = self.compute_mu(y);
        let alpha = self.alpha();
        let mut result = 0.0;
        for i in 0..mu.len() - 1 {
            let dy = y[i] - mu[i];
            let dm = mu[i + 1] - mu[i];
            result += dy * dy + alpha * dm * dm;
        }
        let last = mu.len() - 1;
        let dy = y[last] - mu[last];
        result += dy * dy;
        result
    }
}

/// Cached tridiagonal LU for the L2-Mumford-Shah within-segment solve.
pub struct GaussL2Mum {
    pub n: usize,
    pub alpha: f64,
    /// Diagonal entries `D[k][i]` after LU; `diagonals[k]` has length `k+1`.
    pub diagonals: Vec<Vec<f64>>,
    /// Off-diagonal factors `C[k][i]`; `factors[k]` has length `k` (empty
    /// for `k == 0`).
    pub factors: Vec<Vec<f64>>,
    /// Constant off-diagonal `b = -alpha`. Cached because the back-substitution
    /// uses it on every step.
    pub b: f64,
}

impl GaussL2Mum {
    pub fn new(n: usize, alpha: f64) -> Self {
        let mut me = Self {
            n,
            alpha,
            diagonals: Vec::with_capacity(n),
            factors: Vec::with_capacity(n),
            b: -alpha,
        };
        me.initialize();
        me
    }

    fn initialize(&mut self) {
        let alpha = self.alpha;
        let a = alpha + 1.0;       // boundary diagonal
        let b = -alpha;            // off-diagonal
        let c = 2.0 * alpha + 1.0; // interior diagonal

        self.diagonals.clear();
        self.factors.clear();

        // Length-1 system: trivial.
        self.diagonals.push(vec![1.0]);
        self.factors.push(Vec::new());

        if self.n == 1 {
            return;
        }

        // Length-2 system: just the two boundary rows.
        self.diagonals.push(vec![a, a - b * b / a]);
        self.factors.push(vec![-b / a]);

        // Length-`idx+1` systems for idx in 2..n.
        for idx in 2..self.n {
            let mut d = vec![0.0_f64; idx + 1];
            let mut c_arr = vec![0.0_f64; idx];
            d[0] = a;
            for cur in 1..idx {
                let multi = -b / d[cur - 1];
                c_arr[cur - 1] = multi;
                d[cur] = c + b * multi;
            }
            // Boundary at the bottom.
            let multi = -b / d[idx - 1];
            c_arr[idx - 1] = multi;
            d[idx] = a + b * multi;
            self.diagonals.push(d);
            self.factors.push(c_arr);
        }
    }
}

impl GaussElim for GaussL2Mum {
    fn alpha(&self) -> f64 { self.alpha }

    /// Solve `A μ = y` where `A` is the L2-Mumford-Shah tridiagonal of
    /// length `y.len()` using the cached LU factors. Direct port of
    /// `GaussL2Mum.computeMu`.
    fn compute_mu(&self, y: &[f64]) -> Vec<f64> {
        let len = y.len();
        if len == 0 {
            return Vec::new();
        }
        if len == 1 {
            return vec![y[0]];
        }
        assert!(
            len <= self.n,
            "GaussL2Mum: y.len()={} exceeds cached n={}",
            len,
            self.n
        );
        let diag = &self.diagonals[len - 1];
        let factor = &self.factors[len - 1];
        debug_assert_eq!(diag.len(), len);
        debug_assert_eq!(factor.len(), len - 1);

        // Forward substitution.
        let mut bvec = vec![0.0_f64; len];
        bvec[0] = y[0];
        for i in 1..len {
            bvec[i] = y[i] + factor[i - 1] * bvec[i - 1];
        }

        // Back substitution.
        let mut result = vec![0.0_f64; len];
        result[len - 1] = bvec[len - 1] / diag[len - 1];
        for i in (0..len - 1).rev() {
            // Java: result[i] = (b[i] - this.b * result[i+1]) / diagonal[i]
            // Here `self.b == -alpha`, so subtracting it adds `alpha`.
            result[i] = (bvec[i] - self.b * result[i + 1]) / diag[i];
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference: build the dense tridiagonal matrix and solve with naive
    /// Gauss elimination. Used to cross-check the cached LU.
    fn solve_dense(y: &[f64], alpha: f64) -> Vec<f64> {
        let n = y.len();
        if n == 0 {
            return Vec::new();
        }
        if n == 1 {
            return vec![y[0]];
        }
        let a_d = alpha + 1.0;
        let c_d = 2.0 * alpha + 1.0;
        let b_d = -alpha;
        // Build a dense (n x n+1) augmented matrix and solve.
        let mut m = vec![vec![0.0_f64; n + 1]; n];
        for i in 0..n {
            let diag = if i == 0 || i == n - 1 { a_d } else { c_d };
            m[i][i] = diag;
            if i > 0 { m[i][i - 1] = b_d; }
            if i < n - 1 { m[i][i + 1] = b_d; }
            m[i][n] = y[i];
        }
        // Forward elimination.
        for i in 0..n {
            let pivot = m[i][i];
            for k in i + 1..n {
                let factor = m[k][i] / pivot;
                for j in i..=n {
                    m[k][j] -= factor * m[i][j];
                }
            }
        }
        // Back substitution.
        let mut x = vec![0.0_f64; n];
        for i in (0..n).rev() {
            let mut s = m[i][n];
            for j in i + 1..n {
                s -= m[i][j] * x[j];
            }
            x[i] = s / m[i][i];
        }
        x
    }

    fn assert_close(a: &[f64], b: &[f64], tol: f64) {
        assert_eq!(a.len(), b.len());
        for i in 0..a.len() {
            let diff = (a[i] - b[i]).abs();
            assert!(diff <= tol, "idx {}: {} vs {} (diff {})", i, a[i], b[i], diff);
        }
    }

    #[test]
    fn alpha_zero_is_identity() {
        // With alpha=0 the matrix is the identity → mu == y.
        let g = GaussL2Mum::new(8, 0.0);
        for n in 1..=8 {
            let y: Vec<f64> = (0..n).map(|i| (i as f64) * 0.5 - 1.0).collect();
            let mu = g.compute_mu(&y);
            assert_close(&mu, &y, 1e-12);
        }
    }

    #[test]
    fn matches_dense_gauss_for_small_alpha() {
        let alpha = 0.7;
        let g = GaussL2Mum::new(20, alpha);
        for &y in [
            &[1.0_f64, 2.0, 3.0, 4.0][..],
            &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            &[1.5, -1.5, 0.5, -0.25, 2.0],
            &[3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0],
        ].iter() {
            let mu = g.compute_mu(y);
            let want = solve_dense(y, alpha);
            assert_close(&mu, &want, 1e-10);
        }
    }

    #[test]
    fn compute_dlr_alpha_zero_equals_zero_when_y_const() {
        let g = GaussL2Mum::new(6, 0.0);
        let y = vec![2.0; 6];
        let dlr = g.compute_dlr(&y);
        // alpha=0 means the optimum is mu=y; data error is zero.
        assert!(dlr.abs() < 1e-12, "dlr={}", dlr);
    }

    #[test]
    fn compute_dlr_positive_for_non_const() {
        let g = GaussL2Mum::new(6, 5.0);
        let y = vec![0.0_f64, 1.0, 2.0, 3.0, 4.0, 5.0];
        let dlr = g.compute_dlr(&y);
        // With heavy smoothing, mu won't equal y; cost must be > 0.
        assert!(dlr > 0.0, "dlr={}", dlr);
    }

    /// Cached LU is the same for any sub-length k <= n. A solver built
    /// with capacity n=20 must give identical answers for length-5 inputs
    /// to a solver built with capacity n=5 — the LU is monotone in n.
    #[test]
    fn cached_lu_is_consistent_across_capacities() {
        let alpha = 1.5;
        let g_small = GaussL2Mum::new(5, alpha);
        let g_large = GaussL2Mum::new(20, alpha);
        let y = vec![0.0, 0.5, -0.25, 1.0, 0.75];
        let a = g_small.compute_mu(&y);
        let b = g_large.compute_mu(&y);
        assert_close(&a, &b, 1e-14);
    }

    /// Translation equivariance: solve(y + c) = solve(y) + c. Highly
    /// sensitive to indexing bugs because shifting changes ALL entries.
    #[test]
    fn translation_equivariance() {
        let g = GaussL2Mum::new(8, 0.7);
        let y = vec![0.0, 1.0, -2.0, 0.5, 3.0, -1.0, 0.25, 1.75];
        let c = 7.5;
        let mu = g.compute_mu(&y);
        let y_shifted: Vec<f64> = y.iter().map(|v| v + c).collect();
        let mu_shifted = g.compute_mu(&y_shifted);
        let want: Vec<f64> = mu.iter().map(|v| v + c).collect();
        assert_close(&mu_shifted, &want, 1e-12);
    }

    /// Scale equivariance: solve(s y) = s solve(y) for any s. Both data
    /// and smoothness terms are quadratic in mu, so the optimum is
    /// homogeneous of degree 1.
    #[test]
    fn scale_equivariance() {
        let g = GaussL2Mum::new(8, 1.5);
        let y = vec![0.0, 1.0, -2.0, 0.5, 3.0, -1.0, 0.25, 1.75];
        let s = -3.5;
        let mu = g.compute_mu(&y);
        let y_scaled: Vec<f64> = y.iter().map(|v| s * v).collect();
        let mu_scaled = g.compute_mu(&y_scaled);
        let want: Vec<f64> = mu.iter().map(|v| s * v).collect();
        assert_close(&mu_scaled, &want, 1e-10);
    }

    /// Optimality check: cost at mu* should be <= cost at any small
    /// perturbation. Empirical proof that the solver lands in a minimum.
    #[test]
    fn solution_is_a_local_minimum() {
        let alpha = 1.0;
        let g = GaussL2Mum::new(10, alpha);
        let y = vec![0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 0.5, 1.0, 0.5, 0.0];
        let mu = g.compute_mu(&y);
        let cost = g.compute_dlr(&y);

        // Try 50 random perturbations; cost must always be >= optimum.
        let mut state: u64 = 0xC0FFEE;
        for _ in 0..50 {
            // Trivial LCG to avoid pulling in a PRNG dep.
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let mut perturb = vec![0.0_f64; mu.len()];
            for v in perturb.iter_mut() {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                let r = ((state >> 33) as f64 / (1u64 << 31) as f64) - 1.0;
                *v = r * 1e-3;
            }
            let mu_p: Vec<f64> = mu.iter().zip(perturb.iter()).map(|(m, e)| m + e).collect();
            let mut cost_p = 0.0_f64;
            for i in 0..mu.len() - 1 {
                let dy = y[i] - mu_p[i];
                let dm = mu_p[i + 1] - mu_p[i];
                cost_p += dy * dy + alpha * dm * dm;
            }
            let last = mu.len() - 1;
            cost_p += (y[last] - mu_p[last]).powi(2);
            assert!(
                cost_p >= cost - 1e-9,
                "perturbed cost {} below optimum {}",
                cost_p,
                cost
            );
        }
    }

    /// The LU factor pattern for a length-1 system must be `[1.0]` with
    /// no off-diagonal factors, irrespective of alpha.
    #[test]
    fn length_one_lu_is_identity_for_any_alpha() {
        for &alpha in &[0.0_f64, 0.5, 5.0, 1e6] {
            let g = GaussL2Mum::new(5, alpha);
            assert_eq!(g.diagonals[0], vec![1.0]);
            assert_eq!(g.factors[0].len(), 0);
        }
    }

    /// Diagonal storage layout matches the Java source: diagonals[k] has
    /// length k+1 and factors[k] has length k.
    #[test]
    fn lu_storage_layout_matches_java_convention() {
        let n = 6;
        let g = GaussL2Mum::new(n, 0.4);
        for k in 0..n {
            assert_eq!(g.diagonals[k].len(), k + 1, "diagonals[{}]", k);
            assert_eq!(g.factors[k].len(), k, "factors[{}]", k);
        }
    }
}
