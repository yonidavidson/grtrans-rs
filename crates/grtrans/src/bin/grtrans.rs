//! grtrans-rs command-line driver.
//!
//! Mirrors the upstream `grtrans` program: reads `files.in` (ifile/ofile)
//! and the input deck, runs the ray tracing, and writes the output.
//! Output is plain binary (`cflag=0`) or FITS (`cflag=1`).

use grtrans::{fits, inputs, run};
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
    let images = run::run(&inp);
    if inp.cflag == 1 {
        fits::write_fits(&images, &inp, &PathBuf::from(&ofile)).expect("writing FITS output");
    } else {
        run::write_binary(&images, &PathBuf::from(&ofile)).expect("writing output");
    }
    println!("wrote {} image(s) to {ofile}", images.len());
}
