//! End-to-end CLI test: run the `grtrans` binary on the THINDISK input deck
//! (binary output) and compare with the upstream reference fixture.

use std::io::Read;
use std::path::PathBuf;

fn read_record(f: &mut std::fs::File) -> Vec<u8> {
    let mut m = [0u8; 4];
    f.read_exact(&mut m).unwrap();
    let len = i32::from_le_bytes(m) as usize;
    let mut buf = vec![0u8; len];
    f.read_exact(&mut buf).unwrap();
    let mut m2 = [0u8; 4];
    f.read_exact(&mut m2).unwrap();
    assert_eq!(i32::from_le_bytes(m2) as usize, len);
    buf
}

#[test]
fn cli_thindisk_binary_output_matches_reference() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture_dir = root.join("reference/fixtures/thindisk");
    let work = root.join("target/cli_thindisk_test");
    let _ = std::fs::remove_dir_all(&work);
    std::fs::create_dir_all(&work).unwrap();

    // patch cflag to 0 (binary output) for this test
    let inputs = std::fs::read_to_string(fixture_dir.join("inputs.in")).unwrap();
    let inputs = inputs.replace("cflag=1,", "cflag=0,");
    std::fs::write(work.join("inputs.in"), inputs).unwrap();
    std::fs::write(
        work.join("files.in"),
        "&files\n ifile=\"inputs.in\",\n ofile=\"grtrans_test.bin\",\n/\n",
    )
    .unwrap();

    let status = std::process::Command::new(env!("CARGO_BIN_EXE_grtrans"))
        .arg(work.join("files.in"))
        .current_dir(&work)
        .status()
        .unwrap();
    assert!(status.success(), "grtrans CLI failed");

    // parse the binary output: 25 images of 10000x4
    let mut f = std::fs::File::open(work.join("grtrans_test.bin")).unwrap();
    let reference: Vec<f64> = std::fs::read(fixture_dir.join("ivals.f64.bin"))
        .unwrap()
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
        .collect();
    let npix = 10000usize;
    let nvals = 4usize;
    let nfreq = 25usize;
    let mut ours = vec![0.0f64; npix * nvals * nfreq];
    let mut seen = 0usize;
    while seen < nfreq {
        let hdr = read_record(&mut f);
        let dims: Vec<i32> = hdr
            .chunks_exact(4)
            .map(|c| i32::from_le_bytes(c.try_into().unwrap()))
            .collect();
        assert_eq!(dims, vec![100, 100, 4]);
        let nkey_rec = read_record(&mut f);
        let nkey = i32::from_le_bytes(nkey_rec[..4].try_into().unwrap());
        for _ in 0..nkey {
            let _ = read_record(&mut f);
        }
        let ab = read_record(&mut f);
        assert_eq!(ab.len(), 2 * npix * 4);
        let iv = read_record(&mut f);
        assert_eq!(iv.len(), npix * nvals * 4);
        let vals: Vec<f32> = iv
            .chunks_exact(4)
            .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
            .collect();
        for i in 0..npix {
            for q in 0..nvals {
                ours[(i * nvals + q) * nfreq + seen] = vals[i * nvals + q] as f64;
            }
        }
        seen += 1;
    }
    let num: f64 = ours
        .iter()
        .zip(reference.iter())
        .map(|(a, b)| (a - b).abs())
        .sum();
    let den: f64 = reference.iter().map(|b| b.abs()).sum();
    let err = num / den;
    println!("CLI THINDISK relative error vs reference: {err:e}");
    assert!(err < 1e-2, "CLI THINDISK error {err:e}");
}
