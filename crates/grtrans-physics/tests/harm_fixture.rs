#![allow(clippy::needless_range_loop)]
//! HARM fluid model validated against the Fortran fixture.
use grtrans_core::four_vector::FourVector;
use grtrans_physics::models::harm::{harm_vals, initialize_harm_model};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_harm_fluid.txt"
);

#[test]
fn harm_fluid_matches_fortran() {
    let dump = format!(
        "{}/../../reference/fixtures/harm/dump040",
        env!("CARGO_MANIFEST_DIR")
    );
    let data = initialize_harm_model(std::path::Path::new(&dump));
    let text = std::fs::read_to_string(FIXTURE).expect("harm fluid fixture missing");
    let mut rows: Vec<[f64; 14]> = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() != 14 {
            continue;
        }
        let mut a = [0.0f64; 14];
        let mut ok = true;
        for (i, t) in toks.iter().enumerate() {
            match t.replace(['d', 'D'], "e").parse::<f64>() {
                Ok(v) => a[i] = v,
                Err(_) => {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            rows.push(a);
        }
    }
    assert!(rows.len() >= 100, "only {} rows", rows.len());
    let nr = 12usize;
    let nth = 12usize;
    let mut rgrid = [0.0f64; 12];
    let mut thgrid = [0.0f64; 12];
    for i in 0..nr {
        // descending (matches the fixture driver)
        rgrid[i] = 39.0f64 * (-(i as f64) / ((nr - 1) as f64) * (39.0f64 / 1.5f64).ln()).exp();
    }
    for k in 0..nth {
        thgrid[k] = 0.05f64 + (3.0f64 - 0.05f64) * (k as f64) / ((nth - 1) as f64);
    }
    let mut pts = Vec::with_capacity(nr * nth);
    for k in 0..nth {
        for i in 0..nr {
            pts.push(FourVector::flat([0.0, rgrid[i], thgrid[k], 0.0]));
        }
    }
    let (rho, p, bmag, u, b) = harm_vals(&data, &pts, 0.9375);
    for n in [76usize,77,78,79,80,90] {
        println!("row {n}: r={:.3} th={:.3} | rho {:.6e}/{:.6e} p {:.6e}/{:.6e} bmag {:.6e}/{:.6e}",
            rows[n][12], rows[n][13], rho[n], rows[n][0], p[n], rows[n][1], bmag[n], rows[n][2]);
    }
    let mut worst = 0.0f64;
    let mut desc = String::new();
    for (n, row) in rows.iter().enumerate() {
        let checks: [(&str, f64, f64); 12] = [
            ("rho", rho[n] as f64, row[0]),
            ("p", p[n] as f64, row[1]),
            ("bmag", bmag[n] as f64, row[2]),
            ("bb", b[n].dot(&b[n]), row[3]),
            ("u0", u[n].data[0], row[4]),
            ("u1", u[n].data[1], row[5]),
            ("u2", u[n].data[2], row[6]),
            ("u3", u[n].data[3], row[7]),
            ("b0", b[n].data[0], row[8]),
            ("b1", b[n].data[1], row[9]),
            ("b2", b[n].data[2], row[10]),
            ("b3", b[n].data[3], row[11]),
        ];
        for (name, a, bb) in checks {
            let scale = row[0].abs().max(1.0);
            let e = (a - bb).abs() / scale;
            if e > worst {
                worst = e;
                desc = format!("row {n} {name}: rust {a:e} fort {bb:e} (r={}, th={})", row[12], row[13]);
            }
        }
    }
    println!("HARM fluid worst abs error (scaled): {worst:e} ({desc})");
    assert!(worst < 1e-4, "worst {worst:e}: {desc}");
}
