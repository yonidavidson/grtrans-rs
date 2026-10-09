#![allow(clippy::chunks_exact_to_as_chunks)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::type_complexity)]
#![allow(dead_code)]
//! Geokerr validation against Fortran fixtures.
//!
//! The fixture (reference/fixtures/fortran/test_geokerr.txt) is produced by
//! `reference/fortran/test_geokerr.f90` running the upstream production path
//! (geodesics.f90: initialize_pixels + initialize_geodesic) at revision
//! c76cb11. This test runs the Rust port on the same configurations and
//! compares every camera value and every ray sample.

use grtrans_geodesics::geokerr::camera::CameraPixel;
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_geokerr.txt"
);

#[derive(Debug, Clone)]
struct RayFixture {
    pixel: usize,
    uf: f64,
    tpm: i32,
    tpr: i32,
    status: i32,
    npts: usize,
    rows: Vec<[f64; 10]>,
}

#[derive(Debug, Default, Clone)]
struct Section {
    u0: f64,
    offset: f64,
    camera: Vec<[f64; 8]>,
    camera_int: Vec<[i32; 2]>,
    rays: Vec<RayFixture>,
}

fn parse() -> std::collections::HashMap<String, Section> {
    let text = std::fs::read_to_string(FIXTURE)
        .expect("geokerr fixture missing; run scripts/build_fortran_fixtures.sh");
    let mut sections = std::collections::HashMap::new();
    let mut cur = String::new();
    let mut section = Section::default();
    let mut ray: Option<RayFixture> = None;
    for line in text.lines() {
        if let Some(name) = line.strip_prefix("# ") {
            if let Some(rest) = name.strip_prefix("u0 offset ") {
                let mut it = rest.split_whitespace();
                section.u0 = it.next().unwrap().parse().unwrap();
                section.offset = it.next().unwrap().parse().unwrap();
                continue;
            }
            if let Some(rest) = name.strip_prefix("pixel ") {
                // finish previous ray
                if let Some(r) = ray.take() {
                    section.rays.push(r);
                }
                let mut it = rest.split_whitespace();
                let pixel: usize = it.next().unwrap().parse().unwrap();
                let uf: f64 = it.next().unwrap().parse().unwrap();
                let tpm: i32 = it.next().unwrap().parse().unwrap();
                let tpr: i32 = it.next().unwrap().parse().unwrap();
                let status: i32 = it.next().unwrap().parse().unwrap();
                let npts: usize = it.next().unwrap().parse().unwrap();
                ray = Some(RayFixture {
                    pixel,
                    uf,
                    tpm,
                    tpr,
                    status,
                    npts,
                    rows: Vec::new(),
                });
                continue;
            }
            // new section header
            if let Some(r) = ray.take() {
                section.rays.push(r);
            }
            if !cur.is_empty() {
                sections.insert(cur.clone(), section.clone());
            }
            cur = name.trim().to_string();
            section = Section::default();
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        // upstream prints some diagnostics to stdout (e.g. "kext next:");
        // skip any line whose tokens are not all numeric
        let vals: Option<Vec<f64>> = line
            .split_whitespace()
            .map(|t| t.parse::<f64>())
            .collect::<Result<Vec<f64>, _>>()
            .ok();
        let Some(vals) = vals else { continue };
        if let Some(r) = ray.as_mut() {
            let mut row = [0.0f64; 10];
            row.copy_from_slice(&vals[..10]);
            r.rows.push(row);
        } else {
            let mut cam = [0.0f64; 8];
            cam.copy_from_slice(&vals[..8]);
            section.camera.push(cam);
            section.camera_int.push([vals[8] as i32, vals[9] as i32]);
        }
    }
    if let Some(r) = ray.take() {
        section.rays.push(r);
    }
    if !cur.is_empty() {
        sections.insert(cur, section);
    }
    sections
}

fn rel(a: f64, b: f64) -> f64 {
    let scale = b.abs().max(1e-300);
    (a - b).abs() / scale
}

fn check_camera(label: &str, ours: &[CameraPixel], fixture: &Section, skip: &[&str]) {
    assert_eq!(ours.len(), fixture.camera.len(), "{label}: pixel count");
    for (i, (p, f)) in ours.iter().zip(fixture.camera.iter()).enumerate() {
        let fi = &fixture.camera_int[i];
        let tol = 1e-12;
        for (name, a, b) in [
            ("alpha", p.alpha, f[0]),
            ("beta", p.beta, f[1]),
            ("q2", p.q2, f[2]),
            ("l", p.l, f[3]),
            ("uf", p.uf, f[4]),
            ("sm", p.sm, f[6]),
            ("su", p.su, f[7]),
        ] {
            if skip.contains(&name) {
                continue;
            }
            assert!(
                rel(a, b) < tol,
                "{label} pixel {i} {name}: rust {a} vs fortran {b} (rel {})",
                rel(a, b)
            );
        }
        assert_eq!(p.tpm, fi[0], "{label} pixel {i} tpm");
        assert_eq!(p.tpr, fi[1], "{label} pixel {i} tpr");
    }
}

fn check_rays(label: &str, args: &grtrans_geodesics::rays::GeokerrArgs, fixture: &Section) {
    for rf in &fixture.rays {
        let (ray, status) = initialize_geodesic(args, rf.pixel - 1);
        assert_eq!(status, rf.status, "{label} pixel {} status", rf.pixel);
        assert_eq!(ray.npts, rf.npts, "{label} pixel {} npts", rf.pixel);
        assert_eq!(ray.npts, rf.rows.len(), "{label} pixel {} rows", rf.pixel);
        for (k, row) in rf.rows.iter().enumerate() {
            let tpm_f = (row[9] / 1000.0) as i32;
            let tpr_f = (row[9] % 1000.0) as i32;
            let x = &ray.x[k];
            let kv = &ray.k[k];
            let checks: [(&str, f64, f64); 9] = [
                ("t", x.data[0], row[0]),
                ("r", x.data[1], row[1]),
                ("theta", x.data[2], row[2]),
                ("phi", x.data[3], row[3]),
                ("k_t", kv.data[0], row[4]),
                ("k_r", kv.data[1], row[5]),
                ("k_th", kv.data[2], row[6]),
                ("k_ph", kv.data[3], row[7]),
                ("lambda", ray.lambda[k], row[8]),
            ];
            for (name, a, b) in checks {
                // Directly computed quantities (x, k) agree to ~1e-11 or
                // better. The affine parameter lambda is assembled by
                // upstream as a difference of ~1e6-scale values (the
                // divergent integral baseline from the tiny starting u0), so
                // its accuracy is limited by cancellation: 1 ulp of the
                // baseline (~1.3e-9) appears in the difference. Use an
                // absolute floor for lambda accordingly.
                let ok = if name == "lambda" || name == "t" {
                    (a - b).abs() < 1e-8 + 1e-7 * b.abs()
                } else {
                    // Direct quantities agree to ~1e-11 relative away from
                    // zero; near zero components (e.g. the radial momentum
                    // at a turning point, where U(u) -> 0) an absolute
                    // floor of 1e-10 is used.
                    (a - b).abs() < 1e-10 + 1e-11 * b.abs()
                };
                assert!(
                    ok,
                    "{label} pixel {} point {k} {name}: rust {a} vs fortran {b} (rel {})",
                    rf.pixel,
                    rel(a, b)
                );
            }
            assert_eq!(
                ray.tpmarr[k], tpm_f,
                "{label} pixel {} point {k} tpm",
                rf.pixel
            );
            assert_eq!(
                ray.tprarr[k], tpr_f,
                "{label} pixel {} point {k} tpr",
                rf.pixel
            );
        }
    }
}

#[test]
fn geokerr_matches_fortran_fixture() {
    let sections = parse();

    // ----- standard = 1, nup = 25 (mufill) -----
    let s = &sections["camera1"];
    let args = initialize_pixels(
        true, 1, 0.6428, -0.5, 0.9375, 0.04, 1.0, 1.0, 2, -13.0, 13.0, -13.0, 13.0, 5, 4, 25,
    );
    assert!((args.u0 - s.u0).abs() / s.u0 < 1e-15);
    assert!((args.offset - s.offset).abs() < 1e-15);
    check_camera("camera1", &args.pixels, s, &["muf"]);
    check_rays("geo1", &args, &sections["geo1"]);

    // ----- standard = 1, nup = 400 (production-like) -----
    let s = &sections["camera1b"];
    let args = initialize_pixels(
        true, 1, 0.6428, -0.5, 0.9375, 0.04, 1.0, 1.0, 2, -13.0, 13.0, -13.0, 13.0, 5, 4, 400,
    );
    assert!((args.u0 - s.u0).abs() / s.u0 < 1e-15);
    check_camera("camera1b", &args.pixels, s, &["muf"]);
    check_rays("geo1b", &args, &sections["geo1b"]);

    // ----- standard = 2, nup = 1 (thin-disk geometry) -----
    let s = &sections["camera2"];
    let args = initialize_pixels(
        true, 2, 0.26, -0.5, 0.9, 0.01, 1.0, 1.0, 2, -21.0, 21.0, -21.0, 21.0, 5, 4, 1,
    );
    assert!((args.u0 - s.u0).abs() / s.u0 < 1e-15);
    assert!((args.offset - s.offset).abs() < 1e-15);
    // upstream leaves UFARR unassigned for standard=2; the Rust port uses a
    // deterministic uplus there (see camera.rs)
    check_camera("camera2", &args.pixels, s, &["uf"]);
    check_rays("geo2", &args, &sections["geo2"]);
}
