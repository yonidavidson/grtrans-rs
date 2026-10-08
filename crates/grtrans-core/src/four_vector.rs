//! Four-vector algebra with an attached metric.
//!
//! Direct translation of `class_four_vector.f90` (upstream GRTRANS).
//!
//! Upstream stores the metric as the upper triangle of the symmetric 4x4
//! matrix in the order `[g00, g01, g02, g03, g11, g12, g13, g22, g23, g33]`.
//! A [`Metric`] is that same array; [`FourVector`] carries its own metric so
//! that inner products can be taken anywhere.
//!
//! Note on semantics: upstream's arithmetic results (`a + b`, `a - b`,
//! `r * a`) copy only `.data` and leave `.metric` undefined. The Rust port
//! copies the metric from the left operand (or the only operand) so that
//! arithmetic is always well-defined; scientific code paths always call
//! `assign_metric` before taking inner products of arithmetic results, so
//! this is behavior-preserving.

/// The upper triangle of a symmetric 4x4 metric, upstream ordering:
/// `[g00, g01, g02, g03, g11, g12, g13, g22, g23, g33]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Metric(pub [f64; 10]);

/// Minkowski metric in the upstream component order.
pub const MINKOWSKI: Metric = Metric([-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0]);

impl Metric {
    pub const fn new(m: [f64; 10]) -> Self {
        Metric(m)
    }

    /// Full 4x4 matrix (row-major) reconstructed from the upper triangle.
    pub fn matrix(&self) -> [[f64; 4]; 4] {
        let m = &self.0;
        [
            [m[0], m[1], m[2], m[3]],
            [m[1], m[4], m[5], m[6]],
            [m[2], m[5], m[7], m[8]],
            [m[3], m[6], m[8], m[9]],
        ]
    }
}

/// A four-vector with its metric, upstream `four_Vector`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FourVector {
    pub data: [f64; 4],
    pub metric: Metric,
}

impl FourVector {
    /// Construct with an explicit metric (upstream `four_Vector(values, metric)`).
    pub const fn new(data: [f64; 4], metric: Metric) -> Self {
        FourVector { data, metric }
    }

    /// Construct with the Minkowski default metric (upstream `make_four_Vector`
    /// default and the `four_Vector(values)` constructor).
    pub const fn flat(data: [f64; 4]) -> Self {
        FourVector {
            data,
            metric: MINKOWSKI,
        }
    }

    /// Lower the index: `b_mu = g_{mu nu} a^nu` (upstream `lower`).
    pub fn lower(&self) -> [f64; 4] {
        let m = &self.metric.0;
        let a = &self.data;
        [
            m[0] * a[0] + m[1] * a[1] + m[2] * a[2] + m[3] * a[3],
            m[1] * a[0] + m[4] * a[1] + m[5] * a[2] + m[6] * a[3],
            m[2] * a[0] + m[5] * a[1] + m[7] * a[2] + m[8] * a[3],
            m[3] * a[0] + m[6] * a[1] + m[8] * a[2] + m[9] * a[3],
        ]
    }

    /// Inner product `a * b` using `self`'s metric (upstream `dot_four_Vector`).
    pub fn dot(&self, other: &FourVector) -> f64 {
        let l = self.lower();
        l[0] * other.data[0] + l[1] * other.data[1] + l[2] * other.data[2] + l[3] * other.data[3]
    }

    /// Squared norm `a * a`.
    pub fn norm_sq(&self) -> f64 {
        self.dot(self)
    }

    /// Normalize to unit norm (upstream `normalize_four_Vector`); returns a
    /// zero vector if the norm is below machine epsilon.
    pub fn normalized(&self) -> FourVector {
        let total = self.norm_sq().sqrt();
        let mut out = *self;
        if total < f64::EPSILON {
            out.data = [0.0; 4];
        } else {
            for i in 0..4 {
                out.data[i] = self.data[i] / total;
            }
        }
        out
    }

    /// Assign a new metric (upstream `assign_metric`).
    pub fn assign_metric(&mut self, metric: Metric) {
        self.metric = metric;
    }

    pub fn plus(&self, other: &FourVector) -> FourVector {
        let mut out = *self;
        for i in 0..4 {
            out.data[i] = self.data[i] + other.data[i];
        }
        out
    }

    pub fn minus(&self, other: &FourVector) -> FourVector {
        let mut out = *self;
        for i in 0..4 {
            out.data[i] = self.data[i] - other.data[i];
        }
        out
    }

    pub fn scale(&self, r: f64) -> FourVector {
        let mut out = *self;
        for i in 0..4 {
            out.data[i] = r * self.data[i];
        }
        out
    }

    /// Linear combination `(1 - w) * a + w * b` (Fortran
    /// `(1d0-weight)*a + weight*b`), used for interpolating rays. The metric
    /// is taken from `a`.
    pub fn lerp(a: &FourVector, b: &FourVector, w: f64) -> FourVector {
        let mut out = *a;
        for i in 0..4 {
            out.data[i] = (1.0 - w) * a.data[i] + w * b.data[i];
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minkowski_dot() {
        let a = FourVector::flat([2.0, 1.0, 0.0, 0.0]);
        assert_eq!(a.norm_sq(), -3.0);
        let b = FourVector::flat([1.0, 1.0, 0.0, 0.0]);
        assert_eq!(a.dot(&b), -1.0);
    }

    #[test]
    fn lowering_symmetric_metric() {
        let metric = Metric([-2.0, 0.0, 0.0, 0.5, 3.0, 0.0, 0.0, 4.0, 0.0, 5.0]);
        let v = FourVector::new([1.0, 2.0, 3.0, 4.0], metric);
        let l = v.lower();
        assert_eq!(l[0], -2.0 * 1.0 + 0.5 * 4.0);
        assert_eq!(l[1], 3.0 * 2.0);
        assert_eq!(l[2], 4.0 * 3.0);
        assert_eq!(l[3], 0.5 * 1.0 + 5.0 * 4.0);
    }
}
