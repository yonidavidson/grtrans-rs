#![allow(clippy::approx_constant)]

//! SCHNITTMAN hot-spot fluid model validated against the upstream
//! `schnittman_hotspot_vals` / `get_schnittman_hotspot_fluidvars`.

use grtrans_core::four_vector::FourVector;
use grtrans_physics::fluid::{get_fluid_vars, initialize_fluid_model, LoadedFluid, ModelState};
use grtrans_physics::models::schnittman::{init_schnittman_hotspot, SchnittmanState};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_schnittman_fluid.txt"
);

#[test]
fn schnittman_matches_fortran() {
    let text = std::fs::read_to_string(FIXTURE).expect("schnittman fixture missing");
    let mut rows: Vec<[f64; 11]> = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let t: Vec<f64> = line
            .split_whitespace()
            .filter_map(|x| x.replace(['d', 'D'], "e").parse::<f64>().ok())
            .collect();
        if t.len() == 11 {
            let mut row = [0.0f64; 11];
            row.copy_from_slice(&t);
            rows.push(row);
        }
    }
    let r = [6.0, 6.2, 5.8, 6.0, 6.0, 6.5, 7.0, 8.0, 6.0, 6.0, 3.0, 12.0];
    let th = [
        1.5708, 1.5708, 1.5708, 1.4, 1.7, 1.5708, 1.5708, 1.5708, 0.8, 2.5, 1.5708, 1.5708,
    ];
    let phi = [0.0, 0.1, -0.1, 0.2, -0.2, 0.5, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0];
    let mut state = SchnittmanState::default();
    init_schnittman_hotspot(&mut state, 1.5, 6.0, 1.0e4);
    let loaded = LoadedFluid {
        name: "SCHNITTMAN".to_string(),
        state: ModelState::Schnittman(state),
    };
    let np = r.len();
    let mut f = initialize_fluid_model(&loaded, np);
    let x0: Vec<FourVector> = (0..np)
        .map(|i| FourVector::flat([0.0, r[i], th[i], phi[i]]))
        .collect();
    let k0: Vec<FourVector> = (0..np)
        .map(|_| FourVector::flat([1.0, 0.1, 0.05, 0.2]))
        .collect();
    let mut ll = loaded.clone();
    get_fluid_vars(&mut ll, &x0, &k0, 0.9375, &mut f);
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    let mut worst = 0.0f64;
    let mut desc = String::new();
    for (i, row) in rows.iter().enumerate() {
        let checks: [(&str, f64, f64); 8] = [
            ("rho", f.rho[i] as f64, row[3]),
            ("bmag", f.bmag[i] as f64, row[4]),
            ("u0", f.u[i].data[0], row[5]),
            ("u3", f.u[i].data[3], row[6]),
            ("b0", f.b[i].data[0], row[7]),
            ("b1", f.b[i].data[1], row[8]),
            ("b2", f.b[i].data[2], row[9]),
            ("b3", f.b[i].data[3], row[10]),
        ];
        for (name, a, bb) in checks {
            num += (a - bb).abs();
            den += bb.abs();
            let scale = row[4].abs().max(1.0);
            let e = (a - bb).abs() / scale;
            if e > worst {
                worst = e;
                desc = format!("pt {i} {name}: rust {a:e} fort {bb:e}");
            }
        }
    }
    let agg = num / den.max(1e-300);
    println!("SCHNITTMAN: aggregate={agg:e} maxscaled={worst:e} ({desc})");
    assert!(agg < 1e-5, "SCHNITTMAN aggregate {agg:e} too large");
}
