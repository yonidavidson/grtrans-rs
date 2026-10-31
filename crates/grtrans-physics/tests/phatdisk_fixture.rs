//! PHATDISK tables, interpolation and Keplerian fluid velocity validated
//! against the upstream routines (`test_phatdisk_fluid`).

use grtrans_core::four_vector::FourVector;
use grtrans_physics::fluid::{
    get_fluid_vars, initialize_fluid_model, load_fluid_model, FluidArgs, LoadedFluid, ModelState,
};
use grtrans_physics::models::phatdisk::{init_phatdisk, phatdisk_vals, PhatdiskState};
use grtrans_physics::models::thindisk::ThindiskState;
use std::sync::Arc;

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_phatdisk_fluid.txt"
);

fn state() -> PhatdiskState {
    let thin = ThindiskState {
        mbh: 6.2e8,
        mdot: 1.57e15,
        rin: -1.0,
        rout: 1e8,
    };
    let mut s = PhatdiskState::new(thin);
    s.nw = 20;
    s.wmin = 1e-4;
    s.wmax = 1e4;
    s.nfreq_tab = 10;
    s.fmin = 1e10;
    s.fmax = 1e15;
    s.rmax = 1e4;
    s.nr = 30;
    s.sigt = 0.4;
    s.fcol = 1.7;
    s
}

fn parse() -> (Vec<[f64; 6]>, Vec<[f64; 4]>, Vec<[f64; 6]>) {
    let text = std::fs::read_to_string(FIXTURE).expect("phatdisk fixture missing");
    let mut tables = Vec::new();
    let mut vals = Vec::new();
    let mut fluid = Vec::new();
    let mut section = 0;
    for line in text.lines() {
        match line {
            "# tables" => {
                section = 1;
                continue;
            }
            "# vals" => {
                section = 2;
                continue;
            }
            "# fluid" => {
                section = 3;
                continue;
            }
            _ => {}
        }
        let t: Vec<f64> = line
            .split_whitespace()
            .filter_map(|x| x.replace(['d', 'D'], "e").parse::<f64>().ok())
            .collect();
        match section {
            1 if t.len() == 6 => tables.push([t[0], t[1], t[2], t[3], t[4], t[5]]),
            2 if t.len() == 4 => vals.push([t[0], t[1], t[2], t[3]]),
            3 if t.len() == 6 => fluid.push([t[0], t[1], t[2], t[3], t[4], t[5]]),
            _ => {}
        }
    }
    (tables, vals, fluid)
}

#[test]
fn phatdisk_matches_fortran() {
    let a = 0.9375f32;
    let mut st = state();
    init_phatdisk(&mut st, a);
    let (tables, vals, fluid) = parse();
    let nft = st.nfreq_tab as usize;

    let mut worst_t = 0.0f64;
    let mut sum_t = 0.0f64;
    let mut den_t = 0.0f64;
    for row in &tables {
        let (i, k) = (row[0] as usize - 1, row[1] as usize - 1);
        let pairs = [
            (st.r_tab[i] as f64, row[2]),
            (st.freq_tab[k] as f64, row[3]),
            (st.fnu_tab[i * nft + k] as f64, row[4]),
            (st.om_tab[i] as f64, row[5]),
        ];
        for (ours, fort) in pairs {
            let e = (ours - fort).abs();
            sum_t += e;
            den_t += fort.abs();
            let rel = e / fort.abs().max(1e-30);
            if rel > worst_t {
                worst_t = rel;
            }
        }
    }
    println!(
        "PHATDISK tables: aggregate={:e} maxrel={worst_t:e}",
        sum_t / den_t
    );

    // interpolated values at 12 radii
    let np = 12usize;
    let rgrid: Vec<f32> = (0..np)
        .map(|i| 2.0f32 * ((i as f32) / ((np - 1) as f32) * (50.0f32 / 2.0).ln()).exp())
        .collect();
    let (fnu, omega) = phatdisk_vals(&st, &rgrid, a);
    let mut sum_v = 0.0f64;
    let mut den_v = 0.0f64;
    for (i, row) in vals.iter().enumerate() {
        let pairs = [
            (omega[i] as f64, row[1]),
            (fnu[i * nft] as f64, row[2]),
            (fnu[i * nft + nft - 1] as f64, row[3]),
        ];
        for (ours, fort) in pairs {
            sum_v += (ours - fort).abs();
            den_v += fort.abs();
        }
    }
    println!("PHATDISK vals: aggregate={:e}", sum_v / den_v);

    // full fluid velocity through the dispatcher
    let loaded = LoadedFluid {
        name: "PHATDISK".to_string(),
        state: ModelState::Phatdisk(Arc::new(st.clone())),
    };
    let mut f = initialize_fluid_model(&loaded, np);
    let mut x0 = Vec::new();
    let mut k0 = Vec::new();
    for (i, &r) in rgrid.iter().enumerate() {
        let th = 0.3f64 + 2.5f64 * (i as f64) / ((np - 1) as f64);
        x0.push(FourVector::flat([0.0, r as f64, th, 0.0]));
        k0.push(FourVector::flat([1.0, 0.0, 0.0, 0.0]));
    }
    get_fluid_vars(&mut loaded.clone(), &x0, &k0, a as f64, &mut f);
    let mut sum_f = 0.0f64;
    let mut den_f = 0.0f64;
    for (i, row) in fluid.iter().enumerate() {
        let pairs = [
            (f.u[i].data[0], row[2]),
            (f.u[i].data[3], row[3]),
            (f.fnu[i * nft] as f64, row[4]),
            (f.freq_tab[0] as f64, row[5]),
        ];
        for (ours, fort) in pairs {
            sum_f += (ours - fort).abs();
            den_f += fort.abs();
        }
    }
    println!("PHATDISK fluid: aggregate={:e}", sum_f / den_f);

    assert!(sum_t / den_t < 1e-5, "tables aggregate too large");
    assert!(sum_v / den_v < 1e-5, "vals aggregate too large");
    assert!(sum_f / den_f < 1e-5, "fluid aggregate too large");
    let _ = FluidArgs::default();
    let _ = load_fluid_model;
}
