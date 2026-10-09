//! FFJET fluid model validated against the Fortran fixture.
use grtrans_core::four_vector::FourVector;
use grtrans_physics::models::ffjet::{ffjet_vals, initialize_ffjet_model};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_ffjet_fluid.txt"
);

#[test]
fn ffjet_fluid_matches_fortran() {
    let dump = format!(
        "{}/../../reference/fixtures/ffjet/m87bl09rfp10xi5a998fluidvars.bin",
        env!("CARGO_MANIFEST_DIR")
    );
    let data = initialize_ffjet_model(std::path::Path::new(&dump));
    let text = std::fs::read_to_string(FIXTURE).expect("ffjet fluid fixture missing");
    let mut rows: Vec<[f64; 13]> = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() != 13 {
            continue;
        }
        let mut a = [0.0f64; 13];
        let mut ok = true;
        for (i, t) in toks.iter().enumerate() {
            match t.parse::<f64>() {
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
    let x0: Vec<FourVector> = rows
        .iter()
        .map(|r| {
            let rr = r[12];
            // z from the fixture's construction is implicit; we recompute the
            // same theta by using the fixture's spherical r only for the
            // radius. To keep this exact, the driver stores r = sqrt(rr^2+zz^2)
            // and theta = acos(zz/r); here we invert using the stored
            // four-vector: the fixture does not store theta, so this test
            // instead re-evaluates on the same (r, theta) grid as the driver.
            FourVector::flat([0.0, rr, 0.0, 0.0])
        })
        .collect();
    // The driver used a 12x12 (r,z) grid; rebuild theta from the same grid.
    let nr = 12usize;
    let nz = 12usize;
    let mut rgrid = [0.0f64; 12];
    let mut zgrid = [0.0f64; 12];
    for i in 0..nr {
        rgrid[i] = 1.2f64 * ((i as f64) / ((nr - 1) as f64) * (120.0f64 / 1.2f64).ln()).exp();
    }
    for k in 0..nz {
        zgrid[k] = 0.3f64 * ((k as f64) / ((nz - 1) as f64) * (120.0f64 / 0.3f64).ln()).exp();
    }
    let mut pts = Vec::with_capacity(nr * nz);
    for k in 0..nz {
        for i in 0..nr {
            let rr = rgrid[i];
            let zz = zgrid[k];
            let r = (rr * rr + zz * zz).sqrt();
            pts.push(FourVector::flat([0.0, r, (zz / r).acos(), 0.0]));
        }
    }
    let _ = x0;
    let (rho, p, bmag, u, b) = ffjet_vals(&data, &pts, 0.998);
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
                desc = format!("row {n} {name}: rust {a:e} fort {bb:e}");
            }
        }
    }
    println!("FFJET fluid worst abs error (scaled): {worst:e} ({desc})");
    // The FFJET model is single precision internally (upstream `real`), and
    // gfortran resolves the mixed-kind calls to the single-precision
    // overloads; the Rust port reproduces that. Agreement is therefore at
    // the f32-precision level (~1e-5 worst in a velocity component).
    assert!(worst < 1e-4, "worst {worst:e}: {desc}");
}
