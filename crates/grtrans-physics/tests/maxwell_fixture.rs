//! MAXJUTT and MAXCOMP emissivities validated against the upstream Fortran
//! `calc_maxjutt`/`calc_maxcomp` subroutines (`test_maxwell`).

use grtrans_physics::emissivity::Emis;

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_maxwell.txt"
);

fn load() -> Vec<(String, Vec<[f64; 11]>)> {
    let text = std::fs::read_to_string(FIXTURE).expect("maxwell fixture missing");
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
fn maxwell_matches_fortran() {
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
    let weights = [
        0.10714295, 0.4459917, 1.34219997, 2.17450422, 1.36515333, 0.52892829,
    ];
    let rshift = vec![1.0f64; n];
    let cosne = vec![0.0f64; n];

    for (name, rows) in load() {
        assert_eq!(rows.len(), n, "case {name} row count");
        let (ename, mut otherargs) = match name.as_str() {
            "maxjutt" => ("MAXJUTT", vec![2.0]),
            "maxcomp_sel6" => ("MAXCOMP", vec![2.0, 6.0]),
            "maxcomp_sel0" => ("MAXCOMP", vec![2.0, 0.0]),
            other => panic!("unknown case {other}"),
        };
        otherargs.extend_from_slice(&weights);
        let mut e = Emis::select(ename);
        e.otherargs = otherargs;
        e.initialize(n, &rshift, &th, &cosne);
        e.assign_synch_params(&ne, &ne, &b, &tt);
        e.calc_emissivity(&nu);
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for (i, row) in rows.iter().enumerate() {
            for q in 0..4 {
                num += (e.j[i * 4 + q] - row[q]).abs();
                den += row[q].abs();
            }
            for q in 0..7 {
                num += (e.kcoef[i * 7 + q] - row[4 + q]).abs();
                den += row[4 + q].abs();
            }
        }
        let err = num / den.max(1e-300);
        println!("{name}: relative error vs Fortran = {err:e}");
        assert!(err < 1e-9, "{name}: error {err:e} exceeds 1e-9");
    }
}
