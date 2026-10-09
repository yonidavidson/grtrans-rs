fn main() {
    let inp = grtrans::inputs::read_inputs(std::path::Path::new(
        "reference/fixtures/thindisk/inputs.in",
    ));
    println!(
        "standard={} mumin={} nmu={} spin={} uout={} gridvals={:?} nn={:?}",
        inp.standard, inp.mumin, inp.nmu, inp.spin, inp.uout, inp.gridvals, inp.nn
    );
    println!(
        "fname={} mdot={} mbh={} rin={} rout={} tscl={} rscl={}",
        inp.fname, inp.fmdot, inp.mbh, inp.frin, inp.frout, inp.ftscl, inp.frscl
    );
    println!(
        "ename={} nfreq={} fmin={} fmax={} nvals={} iname={} cflag={}",
        inp.ename, inp.nfreq, inp.fmin, inp.fmax, inp.nvals, inp.iname, inp.cflag
    );
    let fr = grtrans::inputs::freqs(&inp);
    println!("freqs[0..3]={:?} freqs[24]={}", &fr[..3], fr[24]);
    println!(
        "mus={:?} mdots={:?}",
        grtrans::inputs::mus(&inp),
        grtrans::inputs::mdots(&inp)
    );
}
