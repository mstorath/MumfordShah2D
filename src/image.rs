// Ported from Java/src/mumfordShah/Point.java and Image.java
// (Hohm, Storath & Weinmann) by Claude Sonnet coding agent, Anthropic, 2026.
//
// The Java code uses `(channels, rows, cols)` indexing throughout. This
// Rust port keeps the same convention internally so direction-slice math
// transfers literally; the PyO3 boundary in `lib.rs` handles conversion
// to/from the numpy `(rows, cols, channels)` layout.
//
// Several methods on `Image` are unused in Phase 1 — they're consumed by
// the 2D direction processor and outer ADMM in Phase 4. Silence the
// dead-code warnings here so they don't drown out real ones in later
// phases.

#![allow(dead_code)]

use ndarray::{Array2, Array3, ArrayView3, ArrayViewMut3};

/// A point with a left-sided gradient. Used by `TautString` (Phase 2) for
/// convex-hull book-keeping.
///
/// Direct port of `Java/src/mumfordShah/Point.java` (Hohm).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub position: i64,
    pub value: f64,
    pub gradient: f64,
}

impl Point {
    pub fn new(position: i64, value: f64, gradient: f64) -> Self {
        Self { position, value, gradient }
    }

    pub fn position(&self) -> i64 { self.position }
    pub fn value(&self) -> f64 { self.value }
    pub fn gradient(&self) -> f64 { self.gradient }
    pub fn set_gradient(&mut self, g: f64) { self.gradient = g; }
}

/// 3-D image stored as `(channels, rows, cols)`. Mirrors the Java `Image`
/// class — provides row, column, and direction-slice accessors used by the
/// 2-D ADMM (Phase 4).
///
/// Direct port of `Java/src/mumfordShah/Image.java` (Hohm/Storath/Weinmann).
#[derive(Debug, Clone)]
pub struct Image {
    /// Shape `(channels, rows, cols)`.
    pub data: Array3<f64>,
}

impl Image {
    /// Allocate a zero image of the given shape.
    pub fn zeros(channels: usize, n_row: usize, n_col: usize) -> Self {
        Self { data: Array3::zeros((channels, n_row, n_col)) }
    }

    /// Wrap an existing `(channels, rows, cols)` array.
    pub fn from_array(data: Array3<f64>) -> Self { Self { data } }

    pub fn channels(&self) -> usize { self.data.shape()[0] }
    pub fn n_row(&self) -> usize { self.data.shape()[1] }
    pub fn n_col(&self) -> usize { self.data.shape()[2] }

    /// Get a single row across all channels: shape `(channels, n_col)`.
    pub fn get_row(&self, index: usize) -> Array2<f64> {
        let (c, _, n_col) = (self.channels(), self.n_row(), self.n_col());
        let mut out = Array2::<f64>::zeros((c, n_col));
        for d in 0..c {
            for j in 0..n_col {
                out[[d, j]] = self.data[[d, index, j]];
            }
        }
        out
    }

    /// Get a single column across all channels: shape `(channels, n_row)`.
    pub fn get_col(&self, index: usize) -> Array2<f64> {
        let (c, n_row, _) = (self.channels(), self.n_row(), self.n_col());
        let mut out = Array2::<f64>::zeros((c, n_row));
        for d in 0..c {
            for j in 0..n_row {
                out[[d, j]] = self.data[[d, j, index]];
            }
        }
        out
    }

    /// Overwrite a row from a `(channels, n_col)` array.
    pub fn set_row(&mut self, row: &Array2<f64>, index: usize) {
        let (c, n_col) = (self.channels(), self.n_col());
        debug_assert_eq!(row.shape(), &[c, n_col]);
        for d in 0..c {
            for j in 0..n_col {
                self.data[[d, index, j]] = row[[d, j]];
            }
        }
    }

    /// Overwrite a column from a `(channels, n_row)` array.
    pub fn set_col(&mut self, col: &Array2<f64>, index: usize) {
        let (c, n_row) = (self.channels(), self.n_row());
        debug_assert_eq!(col.shape(), &[c, n_row]);
        for d in 0..c {
            for j in 0..n_row {
                self.data[[d, j, index]] = col[[d, j]];
            }
        }
    }

    /// Squared Euclidean error between two images.
    pub fn compute_error(u: &Image, v: &Image) -> f64 {
        assert_eq!(u.data.shape(), v.data.shape(), "Images must have the same shape");
        let mut err = 0.0;
        for ((a, b), _) in u.data.iter().zip(v.data.iter()).zip(0..) {
            let d = a - b;
            err += d * d;
        }
        err
    }

    /// In-place dual-variable update: `self += mu * (u - v)`.
    /// Mirrors `Image.updateDualParam` in the Java source.
    pub fn update_dual_param(&mut self, mu: f64, u: &Image, v: &Image) {
        let shape = self.data.shape().to_vec();
        debug_assert_eq!(u.data.shape(), &shape[..]);
        debug_assert_eq!(v.data.shape(), &shape[..]);
        let s = self.data.as_slice_mut().unwrap();
        let us = u.data.as_slice().unwrap();
        let vs = v.data.as_slice().unwrap();
        for i in 0..s.len() {
            s[i] += mu * (us[i] - vs[i]);
        }
    }

    /// Slice the image along a directional stencil.
    ///
    /// The four supported stencils correspond to the MATLAB `nhood` cell
    /// array in `mumfordShah2D.m`: `[1,0]` (vertical), `[1,1]` (main diag),
    /// `[2,1]`, `[1,2]` (knight moves). The `idx` parameter enumerates
    /// every starting line for the chosen stencil — see
    /// `compute_start_coordinates` for the indexing convention. Returns
    /// shape `(channels, line_length)`.
    pub fn get_direction(&self, direction: [i64; 2], idx: i64) -> Array2<f64> {
        let (x0, y0) = self.compute_start_coordinates(direction, idx);
        let size = self.compute_size(direction, x0, y0);
        let step_x = direction[0];
        let step_y = direction[1];
        let c = self.channels();
        let mut result = Array2::<f64>::zeros((c, size as usize));
        for iter in 0..size {
            for d in 0..c {
                let r = (x0 + iter * step_x) as usize;
                let q = (y0 + iter * step_y) as usize;
                result[[d, iter as usize]] = self.data[[d, r, q]];
            }
        }
        result
    }

    /// Write a slice back along a directional stencil. Mirror of
    /// `get_direction`.
    pub fn set_direction(&mut self, data_slice: &Array2<f64>, direction: [i64; 2], idx: i64) {
        let (x0, y0) = self.compute_start_coordinates(direction, idx);
        let size = self.compute_size(direction, x0, y0);
        let step_x = direction[0];
        let step_y = direction[1];
        let c = self.channels();
        debug_assert_eq!(data_slice.shape()[0], c);
        debug_assert!(data_slice.shape()[1] >= size as usize);
        for iter in 0..size {
            for d in 0..c {
                let r = (x0 + iter * step_x) as usize;
                let q = (y0 + iter * step_y) as usize;
                self.data[[d, r, q]] = data_slice[[d, iter as usize]];
            }
        }
    }

    /// Length of a directional slice starting at (x0, y0).
    pub fn compute_size(&self, direction: [i64; 2], x0: i64, y0: i64) -> i64 {
        let n_row = self.n_row() as i64;
        let n_col = self.n_col() as i64;
        let max_x = ((n_row - x0) as f64 / direction[0] as f64).ceil() as i64;
        let max_y = if direction[1] == 0 {
            i64::MAX
        } else {
            ((n_col - y0) as f64 / direction[1] as f64).ceil() as i64
        };
        max_x.min(max_y)
    }

    /// Starting-coordinate enumeration for the four MATLAB-side stencils.
    ///
    /// **Index convention (verbatim from the Java source):** `idx` enumerates
    /// every distinct starting line for the chosen stencil. For idx in
    /// `[0, n_col)` the start is on the top row at column `idx`. Beyond that,
    /// for `[1,1]` the start moves down the leftmost column; for `[2,1]` and
    /// `[1,2]` see the original `Image.computeStartCoordinates` for the
    /// hand-derived offsets.
    pub fn compute_start_coordinates(&self, direction: [i64; 2], idx: i64) -> (i64, i64) {
        let n_col = self.n_col() as i64;
        if idx < n_col {
            (0, idx)
        } else if direction == [1, 1] {
            (idx - n_col, 0)
        } else if direction == [2, 1] {
            if idx < 2 * n_col {
                (1, idx - n_col)
            } else {
                (2 + (idx - 2 * n_col), 0)
            }
        } else if direction == [1, 2] {
            (1 + (idx - n_col) / 2, (idx - n_col) % 2)
        } else {
            // direction [1,0] only ever uses idx in [0, n_col).
            (0, idx)
        }
    }

    /// Borrow the underlying ndarray.
    pub fn view(&self) -> ArrayView3<'_, f64> { self.data.view() }
    pub fn view_mut(&mut self) -> ArrayViewMut3<'_, f64> { self.data.view_mut() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_basic() {
        let p = Point::new(7, 1.5, -0.25);
        assert_eq!(p.position(), 7);
        assert_eq!(p.value(), 1.5);
        assert_eq!(p.gradient(), -0.25);
        let mut q = p;
        q.set_gradient(0.5);
        assert_eq!(q.gradient(), 0.5);
    }

    #[test]
    fn image_row_col_roundtrip() {
        let mut img = Image::zeros(2, 3, 4);
        for d in 0..2 {
            for r in 0..3 {
                for c in 0..4 {
                    img.data[[d, r, c]] = (d * 100 + r * 10 + c) as f64;
                }
            }
        }
        // Read row 1 then write it back — image unchanged.
        let row1 = img.get_row(1);
        let snapshot = img.data.clone();
        img.set_row(&row1, 1);
        assert_eq!(img.data, snapshot);

        // Read col 2 then write it back — image unchanged.
        let col2 = img.get_col(2);
        img.set_col(&col2, 2);
        assert_eq!(img.data, snapshot);
    }

    #[test]
    fn image_direction_horizontal_idx_in_first_band() {
        // [1,0] stencil with idx<n_col goes straight down at column = idx.
        let mut img = Image::zeros(1, 5, 4);
        for r in 0..5 {
            for c in 0..4 {
                img.data[[0, r, c]] = (r * 10 + c) as f64;
            }
        }
        let line = img.get_direction([1, 0], 2);
        // Start at (0, 2), step (1, 0), size = ceil((5-0)/1)=5.
        assert_eq!(line.shape(), &[1, 5]);
        for r in 0..5 {
            assert_eq!(line[[0, r]], (r * 10 + 2) as f64);
        }
    }

    /// Helper: fill a 1-channel image with cell values `r * 10 + c`.
    fn checker_image(n_row: usize, n_col: usize) -> Image {
        let mut img = Image::zeros(1, n_row, n_col);
        for r in 0..n_row {
            for c in 0..n_col {
                img.data[[0, r, c]] = (r * 10 + c) as f64;
            }
        }
        img
    }

    #[test]
    fn direction_main_diagonal_top_band() {
        // [1,1] stencil, idx in [0, n_col): start at top row, column = idx.
        let img = checker_image(4, 6);

        // idx=0 → start (0, 0), step (1, 1). Length = min(4, 6) = 4.
        let line = img.get_direction([1, 1], 0);
        assert_eq!(line.shape(), &[1, 4]);
        for k in 0..4 {
            assert_eq!(line[[0, k]], (k * 10 + k) as f64);
        }

        // idx=2 → start (0, 2), step (1, 1). Length = min(4, ceil((6-2)/1)) = 4.
        let line = img.get_direction([1, 1], 2);
        assert_eq!(line.shape(), &[1, 4]);
        for k in 0..4 {
            assert_eq!(line[[0, k]], (k * 10 + (k + 2)) as f64);
        }

        // idx=5 → start (0, 5), step (1, 1). Length = min(4, ceil((6-5)/1)) = 1.
        let line = img.get_direction([1, 1], 5);
        assert_eq!(line.shape(), &[1, 1]);
        assert_eq!(line[[0, 0]], 5.0);
    }

    #[test]
    fn direction_main_diagonal_left_band() {
        // [1,1], idx >= n_col: start on left column (idx - n_col, 0).
        // 4×6 image, n_col=6.
        let img = checker_image(4, 6);

        // idx=7 → start (1, 0), length = min(ceil((4-1)/1), 6) = 3.
        let line = img.get_direction([1, 1], 7);
        assert_eq!(line.shape(), &[1, 3]);
        for k in 0..3 {
            // cells (1, 0), (2, 1), (3, 2) → 10, 21, 32
            assert_eq!(line[[0, k]], ((k + 1) * 10 + k) as f64);
        }
    }

    #[test]
    fn direction_knight_2_1_top_band() {
        // [2,1] stencil with idx < n_col: start (0, idx), step (2, 1).
        let img = checker_image(7, 5);  // n_row=7, n_col=5

        // idx=1 → start (0, 1). Length = min(ceil((7-0)/2), ceil((5-1)/1)) = min(4, 4) = 4.
        let line = img.get_direction([2, 1], 1);
        assert_eq!(line.shape(), &[1, 4]);
        // Cells: (0, 1), (2, 2), (4, 3), (6, 4) → 1, 22, 43, 64
        let expected = [1.0, 22.0, 43.0, 64.0];
        for k in 0..4 {
            assert_eq!(line[[0, k]], expected[k]);
        }
    }

    #[test]
    fn direction_knight_2_1_middle_band() {
        // [2,1], idx in [n_col, 2*n_col): start (1, idx - n_col).
        let img = checker_image(7, 5);

        // idx=5 (= n_col) → start (1, 0), step (2, 1). Length = min(ceil((7-1)/2), ceil((5-0)/1)) = min(3, 5) = 3.
        let line = img.get_direction([2, 1], 5);
        assert_eq!(line.shape(), &[1, 3]);
        // Cells: (1, 0), (3, 1), (5, 2) → 10, 31, 52
        let expected = [10.0, 31.0, 52.0];
        for k in 0..3 {
            assert_eq!(line[[0, k]], expected[k]);
        }
    }

    #[test]
    fn direction_knight_1_2_top_band() {
        // [1,2] stencil with idx < n_col: start (0, idx), step (1, 2).
        let img = checker_image(5, 7);  // n_row=5, n_col=7

        // idx=1 → start (0, 1). Length = min(ceil(5/1), ceil((7-1)/2)) = min(5, 3) = 3.
        let line = img.get_direction([1, 2], 1);
        assert_eq!(line.shape(), &[1, 3]);
        // Cells: (0, 1), (1, 3), (2, 5) → 1, 13, 25
        let expected = [1.0, 13.0, 25.0];
        for k in 0..3 {
            assert_eq!(line[[0, k]], expected[k]);
        }
    }

    #[test]
    fn set_direction_round_trips_for_all_stencils() {
        // Read a slice along each stencil and write it back unchanged;
        // the image must be byte-identical afterwards.
        let img = checker_image(6, 6);
        for direction in [[1_i64, 0], [1, 1], [2, 1], [1, 2]] {
            // For each starting column 0..n_col, do a get/set round-trip.
            for idx in 0..6_i64 {
                let mut copy = img.clone();
                let slice = copy.get_direction(direction, idx);
                copy.set_direction(&slice, direction, idx);
                assert_eq!(
                    copy.data, img.data,
                    "round-trip changed image for direction {:?}, idx {}",
                    direction, idx
                );
            }
        }
    }

    #[test]
    fn set_direction_modifies_only_visited_cells() {
        // Writing -1 along a [1,0] vertical line at column 2 must change
        // exactly column 2 and leave all other columns unchanged.
        let img = checker_image(4, 5);
        let mut modified = img.clone();
        let mut new_col = ndarray::Array2::<f64>::from_elem((1, 4), -1.0);
        modified.set_direction(&new_col, [1, 0], 2);
        for r in 0..4 {
            for c in 0..5 {
                let want = if c == 2 { -1.0 } else { img.data[[0, r, c]] };
                assert_eq!(modified.data[[0, r, c]], want, "cell ({}, {})", r, c);
            }
        }
        // Silence unused-mut warning in stable Rust if `new_col` was
        // logically immutable at this point.
        new_col[[0, 0]] = -1.0;
    }

    #[test]
    fn compute_size_handles_horizontal_only_stencil() {
        // [1, 0] never decreases column; size is bounded by n_row only.
        let img = checker_image(7, 3);
        assert_eq!(img.compute_size([1, 0], 0, 0), 7);
        assert_eq!(img.compute_size([1, 0], 3, 0), 4);
    }

    #[test]
    fn update_dual_param_arithmetic() {
        // Image::update_dual_param implements λ += μ (u - v). Verify on a
        // tiny 1x1x1 case where the result is hand-checkable.
        let mut lambda = Image::from_array(ndarray::arr3(&[[[1.0]]]));
        let u = Image::from_array(ndarray::arr3(&[[[5.0]]]));
        let v = Image::from_array(ndarray::arr3(&[[[3.0]]]));
        lambda.update_dual_param(0.5, &u, &v);
        // 1 + 0.5 * (5 - 3) = 2
        assert_eq!(lambda.data[[0, 0, 0]], 2.0);
    }

    #[test]
    fn compute_error_is_squared_l2_distance() {
        let u = Image::from_array(ndarray::arr3(&[[[0.0, 1.0], [2.0, 3.0]]]));
        let v = Image::from_array(ndarray::arr3(&[[[0.0, 1.0], [2.0, 3.0]]]));
        assert_eq!(Image::compute_error(&u, &v), 0.0);

        let v2 = Image::from_array(ndarray::arr3(&[[[1.0, 2.0], [3.0, 4.0]]]));
        // Each entry differs by 1; 4 entries; sum of squares = 4.
        assert_eq!(Image::compute_error(&u, &v2), 4.0);
    }
}
