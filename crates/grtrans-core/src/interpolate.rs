//! Interpolation and array bracketing helpers.
//!
//! Direct translation of `interpolate.f90`, `hunt.f` and `locate.f`
//! (upstream GRTRANS). The bracketing routines reproduce the Numerical
//! Recipes `hunt`/`locate` logic including its exact index conventions;
//! `jlo` is a **1-based** index into the array, exactly as in Fortran, so
//! comparisons with upstream debug output line up.

/// `locate(xx, x)`: find `j` such that `xx[j] <= x < xx[j+1]` for ascending
/// arrays (1-based `j`). Uses pure bisection (no initial guess).
pub fn locate(xx: &[f64], x: f64) -> isize {
    let n = xx.len() as isize;
    let ascnd = xx[xx.len() - 1] >= xx[0];
    let mut jl: isize = 0;
    let mut ju: isize = n + 1;
    while ju - jl > 1 {
        let jm = (ju + jl) / 2;
        if ascnd == (x >= xx[(jm - 1) as usize]) {
            jl = jm;
        } else {
            ju = jm;
        }
    }
    if x == xx[0] {
        1
    } else if x == xx[xx.len() - 1] {
        n - 1
    } else {
        jl
    }
}

/// `hunt(xx, x, jlo)`: bracketing with a starting guess. On entry `jlo`
/// should be the previous bracket (1-based) or any value outside `1..=n` to
/// fall back to bisection. Returns the updated 1-based lower bracket index.
pub fn hunt(xx: &[f64], x: f64, jlo: isize) -> isize {
    let n = xx.len() as isize;
    let ascnd = xx[xx.len() - 1] >= xx[0];
    let mut jlo = jlo;
    let mut jhi;
    if jlo <= 0 || jlo > n {
        jlo = 0;
        jhi = n + 1;
    } else {
        let mut inc: isize = 1;
        if (x >= xx[(jlo - 1) as usize]) == ascnd {
            loop {
                jhi = jlo + inc;
                if jhi > n {
                    jhi = n + 1;
                    break;
                } else if (x < xx[(jhi - 1) as usize]) == ascnd {
                    break;
                } else {
                    jlo = jhi;
                    inc += inc;
                }
            }
        } else {
            jhi = jlo;
            loop {
                jlo = jhi - inc;
                if jlo < 1 {
                    jlo = 0;
                    break;
                } else if (x >= xx[(jlo - 1) as usize]) == ascnd {
                    break;
                } else {
                    jhi = jlo;
                    inc += inc;
                }
            }
        }
    }
    loop {
        if jhi - jlo <= 1 {
            if x == xx[xx.len() - 1] {
                jlo = n - 1;
            }
            if x == xx[0] {
                jlo = 1;
            }
            break;
        } else {
            let jm = (jhi + jlo) / 2;
            if (x >= xx[(jm - 1) as usize]) == ascnd {
                jlo = jm;
            } else {
                jhi = jm;
            }
        }
    }
    jlo
}

/// `get_weight(xx, x, jlo)`: interpolation weight for the bracket
/// containing `x`, updating `jlo` as `hunt` does.
///
/// Returns `(weight, jlo)` with `weight = (x - xx[jlo]) / (xx[jlo+1] - xx[jlo])`
/// (1-based). Panics if the resulting `jlo` is outside the interpolatable
/// range (upstream would read out of bounds).
pub fn get_weight(xx: &[f64], x: f64, jlo: isize) -> (f64, isize) {
    let j = hunt(xx, x, jlo);
    assert!(
        j >= 1 && (j as usize) < xx.len(),
        "get_weight: bracket index {j} out of range for len {} (x={x})",
        xx.len()
    );
    let lo = xx[(j - 1) as usize];
    let hi = xx[j as usize];
    let weight = (x - lo) / (hi - lo);
    (weight, j)
}

/// Vectorized `get_weight_arr`: for a sorted array `x` of query points
/// returns `(indices, weights)` (1-based indices).
pub fn get_weight_arr(xx: &[f64], x: &[f64]) -> (Vec<isize>, Vec<f64>) {
    let mut indx = Vec::with_capacity(x.len());
    let mut jlo = locate(xx, x[0]);
    indx.push(jlo);
    for &xi in &x[1..] {
        jlo = hunt(xx, xi, jlo);
        indx.push(jlo);
    }
    let weight: Vec<f64> = x
        .iter()
        .zip(indx.iter())
        .map(|(&xi, &j)| {
            let lo = xx[(j - 1) as usize];
            let hi = xx[j as usize];
            (xi - lo) / (hi - lo)
        })
        .collect();
    (indx, weight)
}

/// Bilinear interpolation of a value at `(xd, yd)` in `[0,1]^2` given the
/// four corner values in upstream order `v00, v01, v10, v11`
/// (`vv(:,1..4)` in Fortran).
#[inline]
pub fn bilin(v00: f64, v01: f64, v10: f64, v11: f64, xd: f64, yd: f64) -> f64 {
    // w1 = v(:,1)*(1-yd) + v(:,2)*yd ; w2 = v(:,3)*(1-yd) + v(:,4)*yd
    let w1 = v00 * (1.0 - yd) + v01 * yd;
    let w2 = v10 * (1.0 - yd) + v11 * yd;
    w1 * (1.0 - xd) + w2 * xd
}

/// Single-precision bilinear interpolation (upstream `bilininterp` with
/// default `real` arguments, selected by gfortran for FFJET's mixed-kind
/// calls).
#[inline]
pub fn bilin_f32(v00: f32, v01: f32, v10: f32, v11: f32, xd: f32, yd: f32) -> f32 {
    let w1 = v00 * (1.0 - yd) + v01 * yd;
    let w2 = v10 * (1.0 - yd) + v11 * yd;
    w1 * (1.0 - xd) + w2 * xd
}

/// Trilinear interpolation with upstream corner order
/// `v000, v001, v010, v011, v100, v101, v110, v111`.
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn trilin(v: &[f64; 8], xd: f64, yd: f64, zd: f64) -> f64 {
    // first along z
    let i1 = v[0] * (1.0 - zd) + v[1] * zd;
    let i2 = v[2] * (1.0 - zd) + v[3] * zd;
    let j1 = v[4] * (1.0 - zd) + v[5] * zd;
    let j2 = v[6] * (1.0 - zd) + v[7] * zd;
    // then along y
    let w1 = i1 * (1.0 - yd) + i2 * yd;
    let w2 = j1 * (1.0 - yd) + j2 * yd;
    // finally along x
    w1 * (1.0 - xd) + w2 * xd
}

/// Quadrilinear interpolation with upstream corner order (16 corners,
/// Fortran `vv(:,1..16)`).
#[inline]
pub fn quadlin(v: &[f64; 16], td: f64, xd: f64, yd: f64, zd: f64) -> f64 {
    let i1 = v[0] * (1.0 - zd) + v[1] * zd;
    let i2 = v[2] * (1.0 - zd) + v[3] * zd;
    let j1 = v[4] * (1.0 - zd) + v[5] * zd;
    let j2 = v[6] * (1.0 - zd) + v[7] * zd;
    let k1 = v[8] * (1.0 - zd) + v[9] * zd;
    let k2 = v[10] * (1.0 - zd) + v[11] * zd;
    let l1 = v[12] * (1.0 - zd) + v[13] * zd;
    let l2 = v[14] * (1.0 - zd) + v[15] * zd;
    let m1 = i1 * (1.0 - yd) + i2 * yd;
    let m2 = j1 * (1.0 - yd) + j2 * yd;
    let n1 = k1 * (1.0 - yd) + k2 * yd;
    let n2 = l1 * (1.0 - yd) + l2 * yd;
    let w1 = m1 * (1.0 - xd) + m2 * xd;
    let w2 = n1 * (1.0 - xd) + n2 * xd;
    w1 * (1.0 - td) + w2 * td
}

/// Single-precision `hunt` (upstream `hunt` with default `real` arrays).
pub fn hunt_f32(xx: &[f32], x: f32, jlo: isize) -> isize {
    let n = xx.len() as isize;
    let ascnd = xx[xx.len() - 1] >= xx[0];
    let mut jlo = jlo;
    let mut jhi;
    if jlo <= 0 || jlo > n {
        jlo = 0;
        jhi = n + 1;
    } else {
        let mut inc: isize = 1;
        if (x >= xx[(jlo - 1) as usize]) == ascnd {
            loop {
                jhi = jlo + inc;
                if jhi > n {
                    jhi = n + 1;
                    break;
                } else if (x < xx[(jhi - 1) as usize]) == ascnd {
                    break;
                } else {
                    jlo = jhi;
                    inc += inc;
                }
            }
        } else {
            jhi = jlo;
            loop {
                jlo = jhi - inc;
                if jlo < 1 {
                    jlo = 0;
                    break;
                } else if (x >= xx[(jlo - 1) as usize]) == ascnd {
                    break;
                } else {
                    jhi = jlo;
                    inc += inc;
                }
            }
        }
    }
    loop {
        if jhi - jlo <= 1 {
            if x == xx[xx.len() - 1] {
                jlo = n - 1;
            }
            if x == xx[0] {
                jlo = 1;
            }
            break;
        } else {
            let jm = (jhi + jlo) / 2;
            if (x >= xx[(jm - 1) as usize]) == ascnd {
                jlo = jm;
            } else {
                jhi = jm;
            }
        }
    }
    jlo
}

/// Single-precision `get_weight`.
///
/// Note: upstream `hunt` returns a bracket index of 0 for `x` below the
/// first element and `n` for `x` above the last, in which case upstream
/// `get_weight` reads out of bounds. The port clamps to the edge interval
/// (deterministic extrapolation), which only differs for out-of-table
/// inputs; this affects the innermost SPHACC points (documented in
/// docs/VALIDATION_PLAN.md).
pub fn get_weight_f32(xx: &[f32], x: f32, jlo: isize) -> (f32, isize) {
    let n = xx.len() as isize;
    let mut j = hunt_f32(xx, x, jlo);
    if j < 1 {
        j = 1;
    }
    if j > n - 1 {
        j = n - 1;
    }
    let lo = xx[(j - 1) as usize];
    let hi = xx[j as usize];
    let weight = (x - lo) / (hi - lo);
    (weight, j)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locate_ascending() {
        let xx = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(locate(&xx, 1.0), 1);
        assert_eq!(locate(&xx, 2.5), 2);
        assert_eq!(locate(&xx, 4.0), 3);
    }

    #[test]
    fn hunt_matches_locate() {
        let xx: Vec<f64> = (0..100).map(|i| i as f64 * 0.5).collect();
        let mut j = 1isize;
        for i in 0..100 {
            let x = i as f64 * 0.5 + 0.2;
            j = hunt(&xx, x, j);
            assert_eq!(j, locate(&xx, x));
        }
    }

    #[test]
    fn get_weight_interpolates() {
        let xx = [0.0, 1.0, 2.0];
        let (w, j) = get_weight(&xx, 0.25, 0);
        assert_eq!(j, 1);
        assert!((w - 0.25).abs() < 1e-15);
    }

    #[test]
    fn bilinear_corner() {
        assert_eq!(bilin(1.0, 2.0, 3.0, 4.0, 0.0, 0.0), 1.0);
        assert_eq!(bilin(1.0, 2.0, 3.0, 4.0, 1.0, 1.0), 4.0);
        assert_eq!(bilin(1.0, 2.0, 3.0, 4.0, 0.5, 0.0), 2.0);
    }
}
