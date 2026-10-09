//! Polarized power-law synchrotron coefficients validated against the
//! Fortran fixture.
use grtrans_physics::polsynch::{polsynchpl, polsynchth, synchemis, synchpl};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_polsynch.txt"
);

fn inputs() -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let nu = vec![
        2.7555415006877594e11,
        2.683007249434544e11,
        2.731003087487894e11,
        3.45e11,
        3.45e11,
        3.45e11,
        1.0e11,
        1.0e12,
        2.0e11,
        5.0e11,
        3.45e11,
        3.45e11,
    ];
    let nnth = vec![
        4.107198715209961e-1,
        1.6513969898223877e0,
        4.701847076416016e0,
        1.0,
        1.0,
        1.0,
        0.1,
        10.0,
        1.5,
        0.5,
        1.0,
        1.0,
    ];
    let b = vec![
        2.221096083521843e0,
        3.065144419670105e0,
        4.255986772477627e0,
        3.0,
        3.0,
        3.0,
        0.5,
        20.0,
        2.0,
        5.0,
        1.0,
        10.0,
    ];
    let th = vec![
        3.368731396732495e-1,
        2.6194631579424915e-1,
        1.8918080719155692e-1,
        0.3,
        1.2,
        0.05,
        0.5,
        1.0,
        0.8,
        0.2,
        1.5,
        0.1,
    ];
    (nu, nnth, b, th)
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

#[test]
fn polsynchpl_matches_fortran() {
    let text = std::fs::read_to_string(FIXTURE).expect("test_polsynch fixture missing");
    let mut section = "pl";
    let (nu, nnth, b, th) = inputs();
    let p = vec![3.5f64; nu.len()];
    let gmin = vec![100.0f64; nu.len()];
    let ours_pl = polsynchpl(&nu, &nnth, &b, &th, &p, &gmin, 1e5);
    let ours_un = synchpl(&nu, &nnth, &b, &th, &p, &gmin, 1e5);
    let tt = vec![
        1.0e10, 3.0e10, 1.0e11, 1.0e9, 3.0e11, 5.0e9, 2.0e10, 1.0e12, 1.0e8, 4.0e10, 1.0e11, 7.0e9,
    ];
    let ours_th = polsynchth(&nu, &nnth, &b, &tt, &th);
    let ours_mh = synchemis(&nu, &nnth, &b, &tt);
    let mut i = 0usize;
    let mut worst = 0.0f64;
    let mut worst_desc = String::new();
    for line in text.lines() {
        if line.starts_with('#') {
            section = match line {
                "# synchpl" => "synch",
                "# polsynchth" => "th",
                "# synchemis" => "mh",
                _ => section,
            };
            i = 0;
            continue;
        }
        let v: Vec<f64> = line
            .split_whitespace()
            .map(|t| t.parse().unwrap())
            .collect();
        let row = match section {
            "pl" => ours_pl[i],
            "synch" => ours_un[i],
            "th" => ours_th[i],
            _ => ours_mh[i],
        };
        // upstream `synchpl`/`synchemis` only set columns 1 and 5 (jI,
        // alphaI); the remaining columns are left uninitialized upstream.
        let cols: Vec<usize> = if section == "pl" || section == "th" {
            (0..11).collect()
        } else {
            vec![0, 4]
        };
        for q in cols {
            let scale = v
                .iter()
                .cloned()
                .fold(0.0f64, |m, x| m.max(x.abs()))
                .max(1e-300);
            let e = (row[q] - v[q]).abs() / scale;
            if e > worst {
                worst = e;
                worst_desc = format!(
                    "{section} row {i} col {q}: rust {:e} fort {:e}",
                    row[q], v[q]
                );
            }
        }
        i += 1;
    }
    println!("worst relative-to-row-max error: {worst:e} ({worst_desc})");
    assert!(worst < 1e-6, "worst {worst:e}: {worst_desc}");
}

#[test]
fn polsynchpl_point0_matches_fortran_exactly() {
    let text = std::fs::read_to_string(FIXTURE).unwrap();
    let (nu, nnth, b, th) = inputs();
    let p = vec![3.5f64; nu.len()];
    let gmin = vec![100.0f64; nu.len()];
    let ours = polsynchpl(&nu, &nnth, &b, &th, &p, &gmin, 1e5);
    let first: Vec<f64> = text
        .lines()
        .next()
        .unwrap()
        .split_whitespace()
        .map(|t| t.parse().unwrap())
        .collect();
    for q in 0..11 {
        assert!(
            rel(ours[0][q], first[q]) < 1e-10,
            "col {q}: rust {:e} fort {:e}",
            ours[0][q],
            first[q]
        );
    }
}
