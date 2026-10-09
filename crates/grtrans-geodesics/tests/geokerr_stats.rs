//! Diagnostic (non-asserting): report max relative differences between the
//! Rust port and the Fortran geokerr fixture, per field and section.
//! Run with `cargo test -p grtrans-geodesics --test geokerr_stats -- --nocapture`.

use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_geokerr.txt"
);

struct Ray {
    pixel: usize,
    rows: Vec<[f64; 10]>,
}

fn parse() -> std::collections::HashMap<String, (Vec<[f64; 10]>, Vec<Ray>)> {
    let text = std::fs::read_to_string(FIXTURE).unwrap();
    let mut out: std::collections::HashMap<String, (Vec<[f64; 10]>, Vec<Ray>)> =
        std::collections::HashMap::new();
    let mut cur = String::new();
    let mut cam: Vec<[f64; 10]> = Vec::new();
    let mut rays: Vec<Ray> = Vec::new();
    let mut ray: Option<Ray> = None;
    let flush = |out: &mut std::collections::HashMap<String, (Vec<[f64; 10]>, Vec<Ray>)>,
                 cur: &String,
                 cam: &mut Vec<[f64; 10]>,
                 rays: &mut Vec<Ray>,
                 ray: &mut Option<Ray>| {
        if let Some(r) = ray.take() {
            rays.push(r);
        }
        if !cur.is_empty() {
            out.insert(cur.clone(), (std::mem::take(cam), std::mem::take(rays)));
        }
    };
    for line in text.lines() {
        if let Some(name) = line.strip_prefix("# ") {
            if name.starts_with("u0 offset") {
                continue;
            }
            if let Some(rest) = name.strip_prefix("pixel ") {
                if let Some(r) = ray.take() {
                    rays.push(r);
                }
                let mut it = rest.split_whitespace();
                let pixel: usize = it.next().unwrap().parse().unwrap();
                ray = Some(Ray {
                    pixel,
                    rows: Vec::new(),
                });
                continue;
            }
            flush(&mut out, &cur, &mut cam, &mut rays, &mut ray);
            cur = name.trim().to_string();
            continue;
        }
        let vals: Option<Vec<f64>> = line
            .split_whitespace()
            .map(|t| t.parse::<f64>())
            .collect::<Result<_, _>>()
            .ok();
        let Some(vals) = vals else { continue };
        if let Some(r) = ray.as_mut() {
            let mut row = [0.0; 10];
            row.copy_from_slice(&vals[..10]);
            r.rows.push(row);
        } else {
            let mut row = [0.0; 10];
            row.copy_from_slice(&vals[..10]);
            cam.push(row);
        }
    }
    flush(&mut out, &cur, &mut cam, &mut rays, &mut ray);
    out
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

#[test]
fn geokerr_stats() {
    let sections = parse();
    let configs = [
        (
            "camera1", "geo1", 25usize, 1i32, 0.6428, 0.9375, 0.04, 13.0, 5usize, 4usize,
        ),
        (
            "camera1b", "geo1b", 400, 1, 0.6428, 0.9375, 0.04, 13.0, 5, 4,
        ),
        ("camera2", "geo2", 1, 2, 0.26, 0.9, 0.01, 21.0, 5, 4),
    ];
    let names = [
        "t", "r", "theta", "phi", "k_t", "k_r", "k_th", "k_ph", "lambda",
    ];
    for (_camsec, geosec, nup, standard, mu0, spin, uout, half, nro, nphi) in configs {
        let (_, rays) = &sections[geosec];
        let args = initialize_pixels(
            true, standard, mu0, -0.5, spin, uout, 1.0, 1.0, 2, -half, half, -half, half, nro,
            nphi, nup,
        );
        let mut maxerr = [0.0f64; 9];
        let mut argmax = [0usize; 9];
        for rf in rays {
            let (ray, _status) = initialize_geodesic(&args, rf.pixel - 1);
            if rf.pixel == 11 && geosec == "geo1" {
                for (k, row) in rf.rows.iter().enumerate() {
                    let vals = [
                        ray.x[k].data[0],
                        ray.x[k].data[1],
                        ray.x[k].data[2],
                        ray.x[k].data[3],
                        ray.k[k].data[0],
                        ray.k[k].data[1],
                        ray.k[k].data[2],
                        ray.k[k].data[3],
                        ray.lambda[k],
                    ];
                    let mut worst = 0.0f64;
                    let mut wf = 0;
                    for f in 0..9 {
                        let e = rel(vals[f], row[f]);
                        if e > worst {
                            worst = e;
                            wf = f;
                        }
                    }
                    println!(
                        "  p11 k={k:2} worst {:e} f{} | r {:e}/{:e} t {:e}/{:e} lam {:e}/{:e}",
                        worst, wf, vals[1], row[1], vals[0], row[0], vals[8], row[8]
                    );
                }
            }
            for (k, row) in rf.rows.iter().enumerate() {
                let vals = [
                    ray.x[k].data[0],
                    ray.x[k].data[1],
                    ray.x[k].data[2],
                    ray.x[k].data[3],
                    ray.k[k].data[0],
                    ray.k[k].data[1],
                    ray.k[k].data[2],
                    ray.k[k].data[3],
                    ray.lambda[k],
                ];
                for f in 0..9 {
                    let e = rel(vals[f], row[f]);
                    if e > maxerr[f] {
                        maxerr[f] = e;
                        argmax[f] = rf.pixel;
                    }
                }
            }
        }
        println!("== {geosec} ==");
        for f in 0..9 {
            println!(
                "  {:8} max rel err {:e} (pixel {})",
                names[f], maxerr[f], argmax[f]
            );
        }
    }
}
