//! Direction processor: drives `L2MumfordShah1D` over every stripe of a
//! given direction in a 2-D image. Mirrors
//! `Java/src/mumfordShah/L2L2MumfordShahDirectionProcessorFaster` and the
//! parallel-stripe loop in `AbstractDirectionProcessor.run`.
//!
//! Parallelism: stripes are independent (different start positions, no
//! overlapping pixels for the four primary stencils we support — verified
//! by `image::tests::set_direction_modifies_only_visited_cells`). Compute
//! runs in parallel via rayon; write-back is serialised because
//! `Image::set_direction` needs `&mut self`.

use ndarray::Array2;
use rayon::prelude::*;

use crate::gauss_elim::GaussL2Mum;
use crate::image::Image;
use crate::mumfordshah_1d::L2MumfordShah1D;

/// Run the L2 Mumford-Shah 1-D solver on every stripe of `image` in
/// `direction`, modifying `image` in place. `gamma` is the Potts penalty,
/// `alpha` the within-segment smoothness weight.
///
/// Allocates a single shared `GaussL2Mum` of capacity = max(n_row, n_col)
/// so every stripe length 1..=capacity reuses the cached LU factors.
pub fn run_direction(image: &mut Image, direction: [i64; 2], gamma: f64, alpha: f64) {
    let max_size = image.n_row().max(image.n_col());
    if max_size == 0 {
        return;
    }
    let gauss = GaussL2Mum::new(max_size, alpha);
    let n_stripes = image.n_stripes(direction);
    let channels = image.channels();

    // Compute phase — parallel across stripe indices.
    let results: Vec<(i64, Vec<Vec<f64>>)> = (0..n_stripes)
        .into_par_iter()
        .map(|idx| {
            let stripe = image.get_direction(direction, idx);
            // Convert (channels, len) ndarray into Vec<Vec<f64>>.
            let f: Vec<Vec<f64>> = stripe.outer_iter().map(|r| r.to_vec()).collect();
            let solver = L2MumfordShah1D::new(f, gamma, alpha, &gauss);
            (idx, solver.solve())
        })
        .collect();

    // Write-back phase — sequential because Image::set_direction takes &mut self.
    for (idx, channel_buf) in results {
        let len = channel_buf.first().map(|c| c.len()).unwrap_or(0);
        let mut flat = Vec::with_capacity(channels * len);
        for ch in 0..channels {
            flat.extend_from_slice(&channel_buf[ch]);
        }
        let arr = Array2::from_shape_vec((channels, len), flat).expect("rectangular stripe");
        image.set_direction(&arr, direction, idx);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::image::{directions, Image};
    use ndarray::Array3;

    fn make_image(
        channels: usize,
        n_row: usize,
        n_col: usize,
        fill: impl Fn(usize, usize, usize) -> f64,
    ) -> Image {
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
    fn constant_image_unchanged_for_any_direction() {
        for dir in [
            directions::ROW_STEP,
            directions::NW_SE,
            directions::KNIGHT_2_1,
            directions::KNIGHT_1_2,
        ] {
            let mut img = make_image(1, 5, 5, |_, _, _| 7.0);
            run_direction(&mut img, dir, 1.0, 1.0);
            for i in 0..5 {
                for j in 0..5 {
                    assert!(
                        (img.data[[0, i, j]] - 7.0).abs() < 1e-12,
                        "dir={dir:?}, i={i}, j={j}, got {}",
                        img.data[[0, i, j]]
                    );
                }
            }
        }
    }

    #[test]
    fn axis_step_with_alpha_zero_and_large_gamma_smooths_columns_to_data() {
        // Direction [1,0] = vertical stripes (one per column). With α=0 and
        // γ=∞, the within-segment cost is just sum of squared differences
        // and a single segment is forced. The L2-MS minimum at α=0 is the
        // data itself (no smoothing), so the image is unchanged.
        let mut img = make_image(1, 4, 4, |_, i, j| (i + j) as f64);
        let original = img.clone();
        run_direction(&mut img, directions::ROW_STEP, 1e9, 0.0);
        assert_eq!(img.data, original.data);
    }

    #[test]
    fn axis_step_with_alpha_inf_gamma_inf_collapses_columns_to_their_means() {
        // α=very-large, γ=∞: each column smooths to a single segment with
        // very high smoothness, which is the column's mean.
        let mut img = make_image(1, 4, 4, |_, i, j| (i + 4 * j) as f64);
        run_direction(&mut img, directions::ROW_STEP, 1e9, 1e9);
        // Each column's mean: column j has values (j, 1+j, 2+j, 3+j)·... wait — fill is i+4*j.
        // Column j: values are 0+4j, 1+4j, 2+4j, 3+4j. Mean = 1.5 + 4j.
        for j in 0..4 {
            let expected = 1.5 + 4.0 * j as f64;
            for i in 0..4 {
                assert!(
                    (img.data[[0, i, j]] - expected).abs() < 1e-3,
                    "col {j} row {i}: got {}, expected ~{expected}",
                    img.data[[0, i, j]]
                );
            }
        }
    }

    #[test]
    fn step_in_each_column_recovered_at_small_gamma_zero_alpha() {
        // Same step pattern in every column, axis-aligned direction.
        // α=0, γ small: each column splits at the step; output equals input.
        let mut img = make_image(1, 6, 3, |_, i, _| if i < 3 { 0.0 } else { 1.0 });
        let original = img.clone();
        run_direction(&mut img, directions::ROW_STEP, 0.01, 0.0);
        assert_eq!(img.data, original.data);
    }
}
