//! Carlson elliptic integrals validated against the Fortran fixture.
use grtrans_geodesics::geokerr::elliptic::{rc, rd, rf, rj};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_elliptic.txt"
);

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

#[test]
fn carlson_integrals_match_fortran() {
    let text = std::fs::read_to_string(FIXTURE).expect("test_elliptic fixture missing");
    let mut section = "";
    let mut checked = 0usize;
    let mut worst = 0.0f64;
    let mut worst_where = String::new();
    for line in text.lines() {
        if line.starts_with('#') {
            section = match line {
                "# rf" => "rf",
                "# rc" => "rc",
                "# rd" => "rd",
                "# rj" => "rj",
                _ => section,
            };
            continue;
        }
        let v: Vec<f64> = line
            .split_whitespace()
            .map(|t| t.parse().unwrap())
            .collect();
        let (ours, theirs) = match section {
            "rf" => (rf(v[0], v[1], v[2]), v[3]),
            "rc" => (rc(v[0], v[1]), v[2]),
            "rd" => (rd(v[0], v[1], v[2]), v[3]),
            "rj" => (rj(v[0], v[1], v[2], v[3]), v[4]),
            _ => continue,
        };
        let e = rel(ours, theirs);
        if e > worst {
            worst = e;
            worst_where = format!("{section}({:?})", &v[..v.len() - 1]);
        }
        checked += 1;
    }
    assert!(checked >= 140, "only {checked} values checked");
    assert!(
        worst < 1e-14,
        "worst Carlson rel err {worst:e} at {worst_where}"
    );
}
