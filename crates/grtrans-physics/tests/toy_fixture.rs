//! TOY fluid model validated against the upstream `toy_vals` and the
//! contravariant four-velocity construction (the well-defined parts).
//! The upstream aleph/bb are garbage (stale f%u / unassigned bb), so the
//! magnetic field is not compared (we use documented deterministic choices).

use grtrans_core::four_vector::FourVector;
use grtrans_physics::fluid::{get_fluid_vars, initialize_fluid_model, LoadedFluid, ModelState};
use grtrans_physics::models::toy::{init_toy, toy_vals, ToyState};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_toy_fluid.txt"
);

#[test]
fn toy_matches_fortran() {
    let text = std::fs::read_to_string(FIXTURE).expect("toy fixture missing");
    let mut vals_rows: Vec<[f64; 4]> = Vec::new();
    let mut u_rows: Vec<[f64; 4]> = Vec::new();
    let mut section = 0;
    for line in text.lines() {
        match line {
            "# toy_vals" => {
                section = 1;
                continue;
            }
            "# toy_u" => {
                section = 2;
                continue;
            }
            _ => {}
        }
        let t: Vec<f64> = line
            .split_whitespace()
            .filter_map(|x| x.replace(['d', 'D'], "e").parse::<f64>().ok())
            .collect();
        if t.len() == 4 {
            let mut row = [0.0f64; 4];
            row.copy_from_slice(&t);
            if section == 1 {
                vals_rows.push(row);
            } else if section == 2 {
                u_rows.push(row);
            }
        }
    }
    let np = 12usize;
    let r: Vec<f64> = (0..np)
        .map(|i| 2.0 + 18.0 * (i as f64) / ((np - 1) as f64))
        .collect();
    let th: Vec<f64> = (0..np)
        .map(|i| 0.3 + 2.5 * (i as f64) / ((np - 1) as f64))
        .collect();
    let a = 0.9375;
    let mut state = ToyState::default();
    init_toy(&mut state, 3.0e7, 1.7, 2.0, 0);
    let mu: Vec<f64> = th.iter().map(|t| t.cos()).collect();
    let (n, l) = toy_vals(&state, a, &mu, &r);
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for (i, row) in vals_rows.iter().enumerate() {
        num += (n[i] - row[2]).abs();
        den += row[2].abs();
        num += (l[i] - row[3]).abs();
        den += row[3].abs();
    }
    println!("TOY vals: aggregate={:e}", num / den.max(1e-300));
    assert!(num / den < 1e-9, "toy_vals aggregate too large");

    let loaded = LoadedFluid {
        name: "TOY".to_string(),
        state: ModelState::Toy(state),
    };
    let mut f = initialize_fluid_model(&loaded, np);
    let x0: Vec<FourVector> = (0..np)
        .map(|i| FourVector::flat([0.0, r[i], th[i], 0.0]))
        .collect();
    let k0: Vec<FourVector> = (0..np)
        .map(|_| FourVector::flat([1.0, 0.1, 0.05, 0.2]))
        .collect();
    let mut ll = loaded.clone();
    get_fluid_vars(&mut ll, &x0, &k0, a, &mut f);
    let mut num2 = 0.0f64;
    let mut den2 = 0.0f64;
    for (i, row) in u_rows.iter().enumerate() {
        num2 += (f.u[i].data[0] - row[2]).abs();
        den2 += row[2].abs();
        num2 += (f.u[i].data[3] - row[3]).abs();
        den2 += row[3].abs();
    }
    println!("TOY u: aggregate={:e}", num2 / den2.max(1e-300));
    assert!(num2 / den2 < 1e-9, "toy u aggregate too large");
}
