//! BREMS, SYNCHTHBREMS, SYNCHTHAVNOABS and HYBRID emissivities validated
//! against the upstream Fortran routines/compositions (`test_emis_extra`).

use grtrans_physics::emissivity::Emis;

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_emis_extra.txt"
);

fn load() -> Vec<(String, Vec<[f64; 11]>)> {
    let text = std::fs::read_to_string(FIXTURE).expect("emis_extra fixture missing");
    let mut cases: Vec<(String, Vec<[f64; 11]>)> = Vec::new();
    for line in text.lines() {
        if let Some(name) = line.strip_prefix("# ") {
            cases.push((name.trim().to_string(), Vec::new()));
            continue;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() != 11 {
            continue;
        }
        let mut row = [0.0f64; 11];
        for (i, t) in toks.iter().enumerate() {
            row[i] = t.replace(['d', 'D'], "e").parse().unwrap();
        }
        cases.last_mut().unwrap().1.push(row);
    }
    cases
}

#[test]
fn emis_extra_matches_fortran() {
    let n = 12usize;
    let nu: Vec<f64> = vec![
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
    let nnth: Vec<f64> = vec![
        4.107198715209961e-1,
        1.6513969898223877,
        4.701847076416016,
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
    let ne: Vec<f64> = nnth.iter().map(|v| v * 1.0e6).collect();
    let b: Vec<f64> = vec![
        2.221096083521843,
        3.065144419670105,
        4.255986772477627,
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
    let th: Vec<f64> = vec![
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
    let tt: Vec<f64> = vec![
        1.0e10, 3.0e10, 1.0e11, 1.0e9, 3.0e11, 5.0e9, 2.0e10, 1.0e12, 1.0e8, 4.0e10, 1.0e11, 7.0e9,
    ];
    let gmin = vec![100.0f64; n];
    let rshift = vec![1.0f64; n];
    let cosne = vec![0.0f64; n];

    for (name, rows) in load() {
        assert_eq!(rows.len(), n, "case {name} row count");
        let ename = match name.as_str() {
            "brems" => "BREMS",
            "synchthbrems" => "SYNCHTHBREMS",
            "synchthavnoabs" => "SYNCHTHAVNOABS",
            "hybridthpl" => "HYBRIDTHPL",
            "hybridth" => "HYBRIDTH",
            "hybridpl" => "HYBRIDPL",
            other => panic!("unknown case {other}"),
        };
        let mut e = Emis::select(ename);
        e.initialize(n, &rshift, &th, &cosne);
        e.assign_synch_params(&ne, &nnth, &b, &tt);
        e.set_synchpl_model(3.5, &gmin, 1.0e5);
        e.calc_emissivity(&nu);
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for (i, row) in rows.iter().enumerate() {
            if e.neq == 4 {
                for q in 0..4 {
                    num += (e.j[i * 4 + q] - row[q]).abs();
                    den += row[q].abs();
                }
                for q in 0..7 {
                    num += (e.kcoef[i * 7 + q] - row[4 + q]).abs();
                    den += row[4 + q].abs();
                }
            } else {
                num += (e.j[i] - row[0]).abs();
                den += row[0].abs();
                num += (e.kcoef[i] - row[4]).abs();
                den += row[4].abs();
            }
        }
        let err = num / den.max(1e-300);
        println!("{ename}: relative error vs Fortran = {err:e}");
        let tol = if ename.starts_with("HYBRID") {
            1e-9
        } else {
            1e-12
        };
        assert!(err < tol, "{ename}: error {err:e} exceeds {tol:e}");
    }
}
