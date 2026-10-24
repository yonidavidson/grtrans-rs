//! FITS output (cflag=1): structure, positional keywords and big-endian data
//! validated with a minimal independent parser (the same layout assumptions
//! as upstream `grtrans_batch.read_grtrans_output`).

use grtrans::fits::write_fits;
use grtrans::inputs::Inputs;
use grtrans::run::Image;

struct Hdu {
    cards: Vec<String>,
    data: Vec<f32>,
}

fn parse_fits(bytes: &[u8]) -> Vec<Hdu> {
    assert_eq!(bytes.len() % 2880, 0, "file not a multiple of 2880");
    let mut hdus = Vec::new();
    let mut off = 0usize;
    while off < bytes.len() {
        // header: 80-char cards until END, padded to 2880 with spaces
        let mut cards = Vec::new();
        let mut naxis1 = 0usize;
        loop {
            let card = std::str::from_utf8(&bytes[off..off + 80]).unwrap();
            let card = card.to_string();
            if card.starts_with("NAXIS1") {
                naxis1 = card[10..30].trim().parse().unwrap();
            }
            let done = card.starts_with("END");
            cards.push(card);
            off += 80;
            if done {
                break;
            }
        }
        let hrem = off % 2880;
        if hrem != 0 {
            off += 2880 - hrem;
        }
        let nbytes = naxis1 * 4;
        let mut data = Vec::with_capacity(naxis1);
        for k in 0..naxis1 {
            let b = &bytes[off + 4 * k..off + 4 * k + 4];
            data.push(f32::from_be_bytes(b.try_into().unwrap()));
        }
        off += nbytes;
        let drem = off % 2880;
        if drem != 0 {
            off += 2880 - drem;
        }
        hdus.push(Hdu { cards, data });
    }
    hdus
}

#[test]
fn fits_structure_and_data() {
    let nx = 2usize;
    let ny = 3usize;
    let nvals = 4usize;
    let npix = nx * ny;
    let mut images = Vec::new();
    for (k, nu) in [1.0e11, 2.0e11].iter().enumerate() {
        let mut img = Image {
            nx,
            ny,
            nvals,
            nu: *nu,
            ab: vec![[0.0; 2]; npix],
            ivals: vec![0.0; npix * nvals],
        };
        for i in 0..npix {
            img.ab[i] = [i as f32, -(i as f32)];
            for q in 0..nvals {
                img.ivals[i * nvals + q] = (1000 * k + 10 * i + q) as f32;
            }
        }
        images.push(img);
    }
    let inp = Inputs {
        fname: "POWERLAW".to_string(),
        ename: "POLSYNCHTH".to_string(),
        iname: "lsoda".to_string(),
        ..Default::default()
    };
    let path = std::env::temp_dir().join("grtrans_fits_test.fits");
    write_fits(&images, &inp, &path).expect("write fits");
    let bytes = std::fs::read(&path).unwrap();
    let hdus = parse_fits(&bytes);
    assert_eq!(hdus.len(), 3, "primary + two image extensions");
    // primary pixel locations
    assert!(hdus[0].cards[0].starts_with("SIMPLE"));
    let ab = &hdus[0].data;
    assert_eq!(ab.len(), npix * 2);
    for i in 0..npix {
        assert_eq!(ab[2 * i], i as f32);
        assert_eq!(ab[2 * i + 1], -(i as f32));
    }
    for (k, h) in hdus[1..].iter().enumerate() {
        assert!(h.cards[0].starts_with("XTENSION"));
        // positional keywords, exactly as the upstream reader assumes
        assert!(h.cards[6].starts_with("NX"), "card 6 = NX");
        assert!(h.cards[7].starts_with("NY"), "card 7 = NY");
        assert!(h.cards[8].starts_with("NU"), "card 8 = NU");
        assert_eq!(h.cards[6][10..30].trim(), "2.000000E+00");
        assert_eq!(h.cards[7][10..30].trim(), "3.000000E+00");
        let expected_nu = if k == 0 {
            "1.000000E+11"
        } else {
            "2.000000E+11"
        };
        assert_eq!(h.cards[8][10..30].trim(), expected_nu);
        assert_eq!(h.data.len(), npix * nvals);
        for i in 0..npix {
            for q in 0..nvals {
                assert_eq!(h.data[i * nvals + q], (1000 * k + 10 * i + q) as f32);
            }
        }
    }
    let _ = std::fs::remove_file(&path);
}

/// Full CLI path on the THINDISK deck (cflag=1): the FITS file must contain
/// the 25 frequency extensions with data matching the upstream reference.
#[test]
fn cli_fits_thindisk_matches_reference() {
    let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture_dir = root.join("reference/fixtures/thindisk");
    let work = root.join("target/cli_fits_thindisk_test");
    let _ = std::fs::remove_dir_all(&work);
    std::fs::create_dir_all(&work).unwrap();
    std::fs::copy(fixture_dir.join("inputs.in"), work.join("inputs.in")).unwrap();
    std::fs::write(
        work.join("files.in"),
        "&files\n ifile=\"inputs.in\",\n ofile=\"grtrans_test.fits\",\n/\n",
    )
    .unwrap();
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_grtrans"))
        .arg(work.join("files.in"))
        .current_dir(&work)
        .status()
        .unwrap();
    assert!(status.success(), "grtrans CLI failed");
    let bytes = std::fs::read(work.join("grtrans_test.fits")).unwrap();
    let hdus = parse_fits(&bytes);
    let npix = 10000usize;
    let nvals = 4usize;
    let nfreq = 25usize;
    assert_eq!(hdus.len(), nfreq + 1);
    assert_eq!(hdus[0].data.len(), npix * 2);
    let reference: Vec<f64> = std::fs::read(fixture_dir.join("ivals.f64.bin"))
        .unwrap()
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
        .collect();
    assert_eq!(reference.len(), npix * nvals * nfreq);
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for k in 0..nfreq {
        let h = &hdus[k + 1];
        assert_eq!(h.data.len(), npix * nvals, "extension {k} size");
        assert!(h.cards[6].starts_with("NX"));
        assert!(h.cards[7].starts_with("NY"));
        assert!(h.cards[8].starts_with("NU"));
        for i in 0..npix {
            for q in 0..nvals {
                let a = h.data[i * nvals + q] as f64;
                // reference layout is (pixel, stokes, freq) with freq fastest
                let b = reference[(i * nvals + q) * nfreq + k];
                num += (a - b).abs();
                den += b.abs();
            }
        }
    }
    let err = num / den;
    println!("CLI FITS THINDISK relative error vs reference: {err:e}");
    assert!(err < 1e-4, "FITS data deviates from reference: {err:e}");
}
