//! GEOMU validation against the Fortran fixture (pixel 11, ncase=5).
use grtrans_geodesics::geokerr::geomu::{geomu, GeomuInputs};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_geomu.txt"
);

#[test]
fn geomu_matches_fortran() {
    let text = std::fs::read_to_string(FIXTURE).expect("test_geomu fixture missing");
    // The fixture interleaves GEOMU records (# call N + 3 data lines) with
    // GEOPHITIME records (# phitime N + 1 data line).
    let mut lines = text.lines().peekable();
    let mut checked = 0;
    loop {
        // find the next GEOMU record header
        let mut found = false;
        while let Some(l) = lines.next() {
            if l.starts_with("# call") {
                found = true;
                break;
            }
            if l.starts_with("# phitime") {
                lines.next(); // skip the phitime data line
            }
        }
        if !found {
            break;
        }
        let l1 = lines.next().unwrap();
        let v1: Vec<f64> = l1.split_whitespace().map(|t| t.parse().unwrap()).collect();
        let l2 = lines.next().unwrap();
        let v2: Vec<f64> = l2.split_whitespace().map(|t| t.parse().unwrap()).collect();
        let l3 = lines.next().unwrap();
        let v3: Vec<f64> = l3.split_whitespace().map(|t| t.parse().unwrap()).collect();
        // v1: uf, iu, muf, u1, u2, u3, u4
        let uf_f = v1[0];
        let iu_f = v1[1];
        let muf_f = v1[2];
        let u1_f = v1[3];
        let u4_f = v1[6];
        // v2: rffu0, rffu1, rffmu1, rffmu2, rffmu3, iu0, i1mu, i3mu
        let rffu0_f = v2[0];
        let rffmu1_f = v2[2];
        let rffmu2_f = v2[3];
        let rffmu3_f = v2[4];
        let i1mu_f = v2[6];
        let i3mu_f = v2[7];
        // v3: ncase, tpr, h1
        let ncase_f = v3[0] as i32;

        // pixel 11 parameters
        let (u0, mu0, a, l, l2, q2, su, sm) = (
            1.7506389832288786e-7,
            0.6428,
            0.9375,
            0.0,
            0.0,
            1.019_934_310_937_5e1,
            1.0,
            -1.0,
        );
        let out = geomu(
            u0,
            uf_f,
            mu0,
            0.0,
            a,
            l,
            l2,
            q2,
            1,
            0,
            su,
            sm,
            &GeomuInputs::first(),
            true,
            true,
        );
        assert_eq!(out.ncase, ncase_f, "ncase");
        let rel = |x: f64, y: f64| (x - y).abs() / y.abs().max(1e-300);
        for (name, ours, theirs) in [
            ("iu", out.iu, iu_f),
            ("muf", out.muf, muf_f),
            ("u1", out.u1, u1_f),
            ("u4", out.u4, u4_f),
            ("rffu0", out.rffu0, rffu0_f),
            ("rffmu1", out.rffmu1, rffmu1_f),
            ("rffmu2", out.rffmu2, rffmu2_f),
            ("rffmu3", out.rffmu3, rffmu3_f),
            ("i1mu", out.i1mu, i1mu_f),
            ("i3mu", out.i3mu, i3mu_f),
        ] {
            assert!(
                rel(ours, theirs) < 1e-13,
                "{name}: rust {ours} vs fortran {theirs} (rel {})",
                rel(ours, theirs)
            );
        }
        checked += 1;
    }
    assert!(checked >= 5);
}
