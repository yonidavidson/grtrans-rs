//! Minimal FITS writer matching upstream `kwrite_raytrace_camera` (cflag=1).
//!
//! Layout (validated against upstream `grtrans.out` / `sphacc_abs.out`):
//! * primary HDU: `SIMPLE`, `BITPIX=-32`, `NAXIS=1`, `NAXIS1=nx*ny*2`,
//!   `EXTEND`, then the pixel locations (alpha,beta) as big-endian f32.
//! * one `IMAGE` extension per image (frequency): `XTENSION`, `BITPIX`,
//!   `NAXIS`, `NAXIS1=nx*ny*nvals`, `PCOUNT`, `GCOUNT`, then `NX`, `NY`, `NU`
//!   at card positions 6,7,8 (the upstream Python reader
//!   `grtrans_batch.read_grtrans_output` reads them positionally), followed
//!   by the input parameters. Data are big-endian f32.
//!
//! Headers are 80-character cards in 2880-byte blocks; data blocks are padded
//! with zeros. No external FITS dependency.

use crate::inputs::Inputs;
use crate::run::Image;
use std::io::Write;
use std::path::Path;

const BLOCK: usize = 2880;
const CARD: usize = 80;

/// cfitsio-style scientific notation with six decimals and a signed,
/// at-least-two-digit exponent (e.g. `2.000000E+02`, `-5.000000E-01`).
fn fmt_e6(v: f64) -> String {
    if v == 0.0 {
        return "0.000000E+00".to_string();
    }
    let s = format!("{:.6E}", v); // e.g. "2.000000E2" or "-5.000000E-1"
    let (mant, exp) = s.split_once('E').expect("scientific format");
    let e: i32 = exp.parse().expect("exponent");
    format!("{}E{:+03}", mant, e)
}

fn pad_card(s: String) -> [u8; CARD] {
    let mut out = [b' '; CARD];
    let b = s.as_bytes();
    let n = b.len().min(CARD);
    out[..n].copy_from_slice(&b[..n]);
    out
}

/// numeric/bool keyword: value right-justified in the 20-character field.
fn card_num(name: &str, value: &str, comment: &str) -> [u8; CARD] {
    pad_card(format!("{:<8}= {:>20} / {}", name, value, comment))
}

/// string keyword: quoted value left-justified in the 20-character field.
fn card_str(name: &str, value: &str, comment: &str) -> [u8; CARD] {
    pad_card(format!(
        "{:<8}= {:<20} / {}",
        name,
        format!("'{}'", value),
        comment
    ))
}

fn card_hierarch_num(name: &str, value: &str, comment: &str) -> [u8; CARD] {
    pad_card(format!("HIERARCH {} = {:>20} / {}", name, value, comment))
}

fn card_hierarch_str(name: &str, value: &str, comment: &str) -> [u8; CARD] {
    pad_card(format!("HIERARCH {} = '{}' / {}", name, value, comment))
}

struct Header {
    cards: Vec<[u8; CARD]>,
}

impl Header {
    fn new() -> Self {
        Header { cards: Vec::new() }
    }
    fn push(&mut self, card: [u8; CARD]) {
        self.cards.push(card);
    }
    fn write<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        let mut bytes = Vec::with_capacity((self.cards.len() + 1) * CARD);
        for c in &self.cards {
            bytes.extend_from_slice(c);
        }
        bytes.extend_from_slice(&pad_card("END".to_string()));
        let rem = bytes.len() % BLOCK;
        if rem != 0 {
            bytes.resize(bytes.len() + (BLOCK - rem), b' ');
        }
        w.write_all(&bytes)
    }
}

fn write_data<W: Write>(w: &mut W, vals: &[f32]) -> std::io::Result<()> {
    let mut bytes = Vec::with_capacity(vals.len() * 4);
    for v in vals {
        bytes.extend_from_slice(&v.to_be_bytes());
    }
    let rem = bytes.len() % BLOCK;
    if rem != 0 {
        bytes.resize(bytes.len() + (BLOCK - rem), 0);
    }
    w.write_all(&bytes)
}

fn primary_header(npix2: usize) -> Header {
    let mut h = Header::new();
    h.push(card_num(
        "SIMPLE",
        "T",
        "file does conform to FITS standard",
    ));
    h.push(card_num("BITPIX", "-32", "number of bits per data pixel"));
    h.push(card_num("NAXIS", "1", "number of data axes"));
    h.push(card_num(
        "NAXIS1",
        &npix2.to_string(),
        "length of data axis 1",
    ));
    h.push(card_num(
        "EXTEND",
        "T",
        "FITS dataset may contain extensions",
    ));
    h.push(pad_card(
        "COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy"
            .to_string(),
    ));
    h.push(pad_card(
        "COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"
            .to_string(),
    ));
    h
}

fn extension_header(img: &Image, inp: &Inputs) -> Header {
    let mut h = Header::new();
    h.push(pad_card(format!(
        "{:<8}= {:<20} / {}",
        "XTENSION", "'IMAGE   '", "IMAGE extension"
    )));
    h.push(card_num("BITPIX", "-32", "number of bits per data pixel"));
    h.push(card_num("NAXIS", "1", "number of data axes"));
    h.push(card_num(
        "NAXIS1",
        &(img.nx * img.ny * img.nvals).to_string(),
        "length of data axis 1",
    ));
    h.push(card_num("PCOUNT", "0", "required keyword; must = 0"));
    h.push(card_num("GCOUNT", "1", "required keyword; must = 1"));
    // positions 6,7,8: NX, NY, NU (read positionally by the upstream reader).
    // upstream passes these as single-precision kvals, hence the f32 cast.
    h.push(card_num("NX", &fmt_e6(img.nx as f32 as f64), "# x pixels"));
    h.push(card_num("NY", &fmt_e6(img.ny as f32 as f64), "# y pixels"));
    h.push(card_num(
        "NU",
        &fmt_e6(img.nu as f32 as f64),
        "Frequency (Hz)",
    ));
    // input parameters (kwrite_raytrace_camera)
    h.push(card_num("STANDARD", &inp.standard.to_string(), "standard"));
    h.push(card_num("MUMIN", &fmt_e6(inp.mumin), "mumin"));
    h.push(card_num("MUMAX", &fmt_e6(inp.mumax), "mumax"));
    h.push(card_num("NMU", &inp.nmu.to_string(), "nmu"));
    h.push(card_num("PHI0", &fmt_e6(inp.phi0), "phi0"));
    h.push(card_num("SPIN", &fmt_e6(inp.spin), "spin"));
    h.push(card_num("UOUT", &fmt_e6(inp.uout), "uout"));
    h.push(card_num("UIN", &fmt_e6(inp.uin), "uin"));
    h.push(card_num("RCUT", &fmt_e6(inp.rcut), "rcut"));
    h.push(card_num("NROTYPE", &inp.nrotype.to_string(), "nrotype"));
    for k in 0..4 {
        h.push(card_hierarch_num(
            &format!("GRIDVALS({})", k + 1),
            &fmt_e6(inp.gridvals[k]),
            &format!("gridvals({})", k + 1),
        ));
    }
    for k in 0..3 {
        h.push(card_hierarch_num(
            &format!("NN({})", k + 1),
            &inp.nn[k].to_string(),
            &format!("nn({})", k + 1),
        ));
    }
    h.push(card_str("FNAME", inp.fname.trim(), "fname"));
    h.push(card_num("DT", &fmt_e6(inp.dt), "dt"));
    h.push(card_num("NT", &inp.nt.to_string(), "nt"));
    h.push(card_num("NLOAD", &inp.nload.to_string(), "nload"));
    h.push(card_num("NMDOT", &inp.nmdot.to_string(), "nmdot"));
    h.push(card_num("MDOTMIN", &fmt_e6(inp.mdotmin), "mdotmin"));
    h.push(card_num("MDOTMAX", &fmt_e6(inp.mdotmax), "mdotmax"));
    h.push(card_str("ENAME", inp.ename.trim(), "ename"));
    h.push(card_num("MBH", &fmt_e6(inp.mbh), " mbh"));
    h.push(card_num("NFREQ", &inp.nfreq.to_string(), "nfreq"));
    h.push(card_num("FMIN", &fmt_e6(inp.fmin), "fmin"));
    h.push(card_num("FMAX", &fmt_e6(inp.fmax), "fmax"));
    h.push(card_num("MUVAL", &fmt_e6(inp.muval), "muval"));
    h.push(card_num("GMIN", &fmt_e6(inp.gmin), "gmin"));
    h.push(card_num("GMAX", &fmt_e6(inp.gmax), "gmax"));
    h.push(card_num("P1", &fmt_e6(inp.p1), "p1"));
    h.push(card_num("P2", &fmt_e6(inp.p2), "p2"));
    h.push(card_num("JETALPHA", &fmt_e6(inp.jetalpha), "jetalpha"));
    h.push(card_str("STYPE", inp.stype.trim(), "stype"));
    h.push(card_hierarch_str(
        "USE_GEOKERR",
        if inp.use_geokerr { "T" } else { "F" },
        "use_geokerr",
    ));
    h.push(card_num("NVALS", &inp.nvals.to_string(), "nvals"));
    h.push(card_str("INAME", inp.iname.trim(), "iname"));
    h.push(card_num("CFLAG", &inp.cflag.to_string(), "cflag"));
    h.push(card_num("EXTRA", &inp.extra.to_string(), "extra"));
    h
}

/// Write all images to a FITS file (`cflag=1`). The first image supplies the
/// primary pixel-location HDU; every image becomes one `IMAGE` extension.
pub fn write_fits(images: &[Image], inp: &Inputs, path: &Path) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    let first = images.first().expect("at least one image");
    primary_header(first.ab.len() * 2).write(&mut f)?;
    let ab: Vec<f32> = first.ab.iter().flat_map(|p| [p[0], p[1]]).collect();
    write_data(&mut f, &ab)?;
    for img in images {
        extension_header(img, inp).write(&mut f)?;
        write_data(&mut f, &img.ivals)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::fmt_e6;

    #[test]
    fn e6_format() {
        assert_eq!(fmt_e6(200.0), "2.000000E+02");
        assert_eq!(fmt_e6(3.45e11), "3.450000E+11");
        assert_eq!(fmt_e6(-0.5), "-5.000000E-01");
        assert_eq!(fmt_e6(0.0), "0.000000E+00");
        assert_eq!(fmt_e6(1.0), "1.000000E+00");
    }
}
