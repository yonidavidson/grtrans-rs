#![allow(clippy::needless_range_loop)]
//! SARIAF fluid model validated against the Fortran fixture.
use grtrans_core::four_vector::FourVector;
use grtrans_physics::fluid::{initialize_fluid_model, LoadedFluid, ModelState};
use grtrans_physics::models::sariaf::{get_sariaf_fluidvars, init_sariaf, SariafState};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_sariaf_fluid.txt"
);

#[test]
fn sariaf_fluid_matches_fortran() {
    let text = std::fs::read_to_string(FIXTURE).expect("sariaf fluid fixture missing");
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
    assert!(rows.len() >= 100);
    let nr = 12usize;
    let nth = 12usize;
    let mut rgrid = [0.0f64; 12];
    let mut thgrid = [0.0f64; 12];
    for i in 0..nr {
        rgrid[i] = 1.5f64 * ((i as f64) / ((nr - 1) as f64) * (40.0f64 / 1.5f64).ln()).exp();
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
    let mut state = SariafState::default();
    init_sariaf(&mut state, 3.0e7, 1.7e11, 8.0e4, 2.9, 10.0, 0);
    let loaded = LoadedFluid {
        name: "SARIAF".to_string(),
        state: ModelState::Sariaf(state),
    };
    let mut f = initialize_fluid_model(&loaded, nr * nth);
    get_sariaf_fluidvars(&state, &pts, 0.9375, &mut f);
    let mut worst = 0.0f64;
    let mut desc = String::new();
    for (n, row) in rows.iter().enumerate() {
        let checks: [(&str, f64, f64); 12] = [
            ("rho", f.rho[n] as f64, row[0]),
            ("p", f.p[n] as f64, row[1]),
            ("bmag", f.bmag[n] as f64, row[2]),
            ("rho2", f.rho2[n] as f64, row[3]),
            ("u0", f.u[n].data[0], row[4]),
            ("u1", f.u[n].data[1], row[5]),
            ("u2", f.u[n].data[2], row[6]),
            ("u3", f.u[n].data[3], row[7]),
            ("b0", f.b[n].data[0], row[8]),
            ("b1", f.b[n].data[1], row[9]),
            ("b2", f.b[n].data[2], row[10]),
            ("b3", f.b[n].data[3], row[11]),
        ];
        for (name, a, bb) in checks {
            let scale = row[0].abs().max(1.0);
            let e = (a - bb).abs() / scale;
            if e > worst {
                worst = e;
                desc = format!(
                    "row {n} {name}: rust {a:e} fort {bb:e} (r={}, th={})",
                    row[12], row[13]
                );
            }
        }
    }
    println!("SARIAF fluid worst abs error (scaled): {worst:e} ({desc})");
    assert!(worst < 1e-5, "worst {worst:e}: {desc}");
}
