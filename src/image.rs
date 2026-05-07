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
}
