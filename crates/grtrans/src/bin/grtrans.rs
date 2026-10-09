//! grtrans-rs command-line driver.
//!
//! Mirrors the upstream `grtrans` program: reads `files.in` (ifile/ofile)
//! and the input deck, runs the ray tracing, and writes the output.
//! Binary output (`cflag=0`) is implemented; FITS output (`cflag=1`) is not
//! ported yet (use cflag=0 or the Python bindings).

use grtrans::{inputs, run};
use std::path::PathBuf;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let files_in = if args.len() > 1 {
        PathBuf::from(&args[1])
    } else {
        PathBuf::from("files.in")
    };
    let (ifile, ofile) = inputs::read_files_in(&files_in);
    println!("grtrans-rs: ifile={ifile} ofile={ofile}");
    let inp = inputs::read_inputs(&PathBuf::from(&ifile));
    if inp.cflag == 1 {
        eprintln!(
            "warning: FITS output (cflag=1) is not ported; writing the plain \
             binary format instead (see docs/PORTING_MATRIX.md)"
        );
    }
    let images = run::run(&inp);
    run::write_binary(&images, &PathBuf::from(&ofile)).expect("writing output");
    println!("wrote {} image(s) to {ofile}", images.len());
}
