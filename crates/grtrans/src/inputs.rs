//! Input-file handling: Fortran namelist parsing for `files.in` and the
//! GRTRANS input deck.
//!
//! Direct translation of the input handling in `grtrans_program.f90` +
//! `read_inputs.f90` (upstream GRTRANS). The namelist parser supports the
//! subset of Fortran namelist syntax used by GRTRANS: scalar/array values,
//! strings in single or double quotes, logicals (`T`/`F`/`.TRUE.`), and
//! `d`-exponent reals.

use std::collections::HashMap;
use std::path::Path;

/// One parsed namelist block: key -> raw string values.
pub type Block = HashMap<String, Vec<String>>;

/// Parse a Fortran namelist file into blocks keyed by namelist name.
pub fn parse_namelists(text: &str) -> HashMap<String, Block> {
    // strip comments
    let mut cleaned = String::new();
    let mut in_quote: Option<char> = None;
    for line in text.lines() {
        for ch in line.chars() {
            if let Some(q) = in_quote {
                cleaned.push(ch);
                if ch == q {
                    in_quote = None;
                }
            } else if ch == '\'' || ch == '"' {
                in_quote = Some(ch);
                cleaned.push(ch);
            } else if ch == '!' {
                break;
            } else {
                cleaned.push(ch);
            }
        }
        cleaned.push('\n');
    }
    let mut blocks: HashMap<String, Block> = HashMap::new();
    let bytes: Vec<char> = cleaned.chars().collect();
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == '&' {
            // read the namelist name
            let mut name = String::new();
            i += 1;
            while i < bytes.len() && bytes[i].is_alphanumeric() {
                name.push(bytes[i]);
                i += 1;
            }
            // find the terminating '/' outside quotes
            let start = i;
            let mut end = i;
            let mut q: Option<char> = None;
            while end < bytes.len() {
                let c = bytes[end];
                if let Some(qq) = q {
                    if c == qq {
                        q = None;
                    }
                } else if c == '\'' || c == '"' {
                    q = Some(c);
                } else if c == '/' {
                    break;
                }
                end += 1;
            }
            let body: String = bytes[start..end].iter().collect();
            let block = parse_block(&body);
            blocks.insert(name.to_lowercase(), block);
            i = end + 1;
        } else {
            i += 1;
        }
    }
    blocks
}

fn parse_block(body: &str) -> Block {
    let mut out: Block = HashMap::new();
    // split on commas outside quotes
    let mut items: Vec<String> = Vec::new();
    let mut cur = String::new();
    let mut q: Option<char> = None;
    for ch in body.chars() {
        if let Some(qq) = q {
            cur.push(ch);
            if ch == qq {
                q = None;
            }
        } else if ch == '\'' || ch == '"' {
            q = Some(ch);
            cur.push(ch);
        } else if ch == ',' || ch == '\n' {
            if !cur.trim().is_empty() {
                items.push(cur.clone());
            }
            cur.clear();
        } else {
            cur.push(ch);
        }
    }
    if !cur.trim().is_empty() {
        items.push(cur);
    }
    // array values are comma-separated and were split above; values without
    // '=' continue the most recently seen key (insertion order matters)
    let mut last_key: Option<String> = None;
    for item in items {
        if let Some(eq) = item.find('=') {
            let key = item[..eq].trim().to_lowercase();
            let val = item[eq + 1..].trim().to_string();
            if key.is_empty() {
                if let Some(last) = &last_key {
                    out.get_mut(last).unwrap().push(val);
                }
            } else {
                out.insert(key.clone(), vec![val]);
                last_key = Some(key);
            }
        } else if let Some(last) = &last_key {
            out.get_mut(last).unwrap().push(item.trim().to_string());
        }
    }
    out
}

fn as_f64(s: &str) -> f64 {
    let t = s.trim().trim_end_matches('_').replace(['d', 'D'], "e");
    t.parse::<f64>()
        .unwrap_or_else(|_| panic!("namelist: cannot parse real from {s:?}"))
}

fn as_i64(s: &str) -> i64 {
    let t = s.trim();
    // tolerate Fortran-style integer suffixes and "1.0" forms
    t.parse::<i64>()
        .or_else(|_| t.parse::<f64>().map(|v| v as i64))
        .unwrap_or_else(|_| panic!("namelist: cannot parse integer from {s:?}"))
}

fn as_bool(s: &str) -> bool {
    let t = s.trim().to_uppercase();
    matches!(t.as_str(), "T" | ".TRUE." | "TRUE")
}

fn as_string(s: &str) -> String {
    s.trim().trim_matches(|c| c == '\'' || c == '"').to_string()
}

/// Typed accessor for a namelist block.
pub struct Nml<'a> {
    pub block: &'a Block,
}

impl<'a> Nml<'a> {
    pub fn new(block: &'a Block) -> Self {
        Nml { block }
    }
    fn get(&self, key: &str) -> Option<&Vec<String>> {
        self.block.get(&key.to_lowercase())
    }
    pub fn f64(&self, key: &str, default: f64) -> f64 {
        self.get(key).map(|v| as_f64(&v[0])).unwrap_or(default)
    }
    pub fn i64(&self, key: &str, default: i64) -> i64 {
        self.get(key).map(|v| as_i64(&v[0])).unwrap_or(default)
    }
    pub fn bool(&self, key: &str, default: bool) -> bool {
        self.get(key).map(|v| as_bool(&v[0])).unwrap_or(default)
    }
    pub fn string(&self, key: &str, default: &str) -> String {
        self.get(key)
            .map(|v| as_string(&v[0]))
            .unwrap_or_else(|| default.to_string())
    }
    pub fn f64_array<const N: usize>(&self, key: &str, default: [f64; N]) -> [f64; N] {
        match self.get(key) {
            None => default,
            Some(v) => {
                let mut out = default;
                for (i, s) in v.iter().take(N).enumerate() {
                    out[i] = as_f64(s);
                }
                out
            }
        }
    }
    pub fn i64_array<const N: usize>(&self, key: &str, default: [i64; N]) -> [i64; N] {
        match self.get(key) {
            None => default,
            Some(v) => {
                let mut out = default;
                for (i, s) in v.iter().take(N).enumerate() {
                    out[i] = as_i64(s);
                }
                out
            }
        }
    }
}

/// Read `files.in` (upstream `grtrans_program.f90`).
pub fn read_files_in(path: &Path) -> (String, String) {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    let blocks = parse_namelists(&text);
    let b = blocks
        .get("files")
        .unwrap_or_else(|| panic!("{}: missing &files namelist", path.display()));
    let n = Nml::new(b);
    (
        n.string("ifile", "inputs.in"),
        n.string("ofile", "grtrans.out"),
    )
}

/// The complete GRTRANS input deck (upstream `grtrans_inputs`).
#[derive(Clone, Debug)]
pub struct Inputs {
    // geodata
    pub standard: i32,
    pub mumin: f64,
    pub mumax: f64,
    pub nmu: i64,
    pub phi0: f64,
    pub spin: f64,
    pub uout: f64,
    pub uin: f64,
    pub rcut: f64,
    pub nrotype: i32,
    pub gridvals: [f64; 4],
    pub nn: [i64; 3],
    pub i1: i64,
    pub i2: i64,
    pub extra: i32,
    pub debug: i32,
    // fluiddata
    pub fname: String,
    pub dt: f64,
    pub nt: i64,
    pub nload: i64,
    pub nmdot: i64,
    pub mdotmin: f64,
    pub mdotmax: f64,
    pub sigcut: f64,
    // emisdata
    pub ename: String,
    pub mbh: f64,
    pub nfreq: i64,
    pub fmin: f64,
    pub fmax: f64,
    pub muval: f64,
    pub gmin: f64,
    pub gmax: f64,
    pub p1: f64,
    pub p2: f64,
    pub jetalpha: f64,
    pub stype: String,
    pub delta: f64,
    pub nweights: i64,
    pub coefindx: [i64; 7],
    // general
    pub use_geokerr: bool,
    pub nvals: i64,
    pub iname: String,
    pub cflag: i32,
    // harm/fluid-file arguments
    pub fdfile: String,
    pub fhfile: String,
    pub fgfile: String,
    pub fsim: String,
    pub fnt: i64,
    pub findf: i64,
    pub fnfiles: i64,
    pub fjonfix: i64,
    pub foffset: i64,
    pub fdindf: i64,
    pub fmagcrit: i64,
    pub fscalefac: f64,
    // analytic arguments
    pub fnw: i64,
    pub fwmin: f64,
    pub fwmax: f64,
    pub fnfreq_tab: i64,
    pub ffmin: f64,
    pub ffmax: f64,
    pub frmax: f64,
    pub fnr: i64,
    pub fsigt: f64,
    pub ffcol: f64,
    pub frspot: f64,
    pub fr0spot: f64,
    pub fn0spot: f64,
    pub ftscl: f64,
    pub frscl: f64,
    pub fmdot: f64,
    pub fnscl: f64,
    pub fnnthscl: f64,
    pub fnnthp: f64,
    pub fbeta: f64,
    pub fbl06: i64,
    pub fnp: f64,
    pub ftp: f64,
    pub frin: f64,
    pub frout: f64,
    pub fthin: f64,
    pub fthout: f64,
    pub fphiin: f64,
    pub fphiout: f64,
}

impl Default for Inputs {
    fn default() -> Self {
        // defaults from grtrans_batch.py init()
        Inputs {
            standard: 2,
            mumin: 0.1,
            mumax: 1.0,
            nmu: 10,
            phi0: -0.5,
            spin: 0.998,
            uout: 0.04,
            uin: 1.0,
            rcut: 1.0,
            nrotype: 2,
            gridvals: [-25.0, 25.0, -25.0, 25.0],
            nn: [100, 100, 1],
            i1: 1,
            i2: -1,
            extra: 0,
            debug: 0,
            fname: "THINDISK".to_string(),
            dt: 5.0,
            nt: 1,
            nload: 1,
            nmdot: 1,
            mdotmin: 1.5e15,
            mdotmax: 1.5e15,
            sigcut: 1e10,
            ename: "BB".to_string(),
            mbh: 10.0,
            nfreq: 1,
            fmin: 1e17,
            fmax: 3e19,
            muval: 0.25,
            gmin: 100.0,
            gmax: 1e5,
            p1: 3.5,
            p2: 3.5,
            jetalpha: 0.02,
            stype: "const".to_string(),
            delta: 1.0,
            nweights: 1,
            coefindx: [1, 1, 1, 1, 1, 1, 1],
            use_geokerr: true,
            nvals: 1,
            iname: "lsoda".to_string(),
            cflag: 1,
            fdfile: String::new(),
            fhfile: String::new(),
            fgfile: String::new(),
            fsim: String::new(),
            fnt: 1,
            findf: 1,
            fnfiles: 1,
            fjonfix: 1,
            foffset: 0,
            fdindf: 1,
            fmagcrit: 0,
            fscalefac: 1.0,
            fnw: 500,
            fwmin: 1e-4,
            fwmax: 1e4,
            fnfreq_tab: 100,
            ffmin: 3.33e16,
            ffmax: 9e19,
            frmax: 1e4,
            fnr: 500,
            fsigt: 0.4,
            ffcol: 1.7,
            frspot: 1.5,
            fr0spot: 6.0,
            fn0spot: 1e4,
            ftscl: 1.0,
            frscl: 6.0,
            fmdot: 0.1,
            fnscl: 3e7,
            fnnthscl: 8e4,
            fnnthp: 2.9,
            fbeta: 10.0,
            fbl06: 0,
            fnp: 0.0,
            ftp: 0.0,
            frin: -1.0,
            frout: 1e8,
            fthin: -10.0,
            fthout: 10.0,
            fphiin: 0.0,
            fphiout: 1e4,
        }
    }
}

/// Read the input deck (upstream `read_inputs`).
pub fn read_inputs(path: &Path) -> Inputs {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    let blocks = parse_namelists(&text);
    let mut inp = Inputs::default();
    if let Some(b) = blocks.get("geodata") {
        let n = Nml::new(b);
        inp.standard = n.i64("standard", inp.standard as i64) as i32;
        inp.mumin = n.f64("mumin", inp.mumin);
        inp.mumax = n.f64("mumax", inp.mumax);
        inp.nmu = n.i64("nmu", inp.nmu);
        inp.phi0 = n.f64("phi0", inp.phi0);
        inp.spin = n.f64("spin", inp.spin);
        inp.uout = n.f64("uout", inp.uout);
        inp.uin = n.f64("uin", inp.uin);
        inp.rcut = n.f64("rcut", inp.rcut);
        inp.nrotype = n.i64("nrotype", inp.nrotype as i64) as i32;
        inp.gridvals = n.f64_array("gridvals", inp.gridvals);
        inp.nn = n.i64_array("nn", inp.nn);
        inp.i1 = n.i64("i1", inp.i1);
        inp.i2 = n.i64("i2", inp.i2);
        inp.extra = n.i64("extra", inp.extra as i64) as i32;
        inp.debug = n.i64("debug", inp.debug as i64) as i32;
    }
    if let Some(b) = blocks.get("fluiddata") {
        let n = Nml::new(b);
        inp.fname = n.string("fname", &inp.fname);
        inp.dt = n.f64("dt", inp.dt);
        inp.nt = n.i64("nt", inp.nt);
        inp.nload = n.i64("nload", inp.nload);
        inp.nmdot = n.i64("nmdot", inp.nmdot);
        inp.mdotmin = n.f64("mdotmin", inp.mdotmin);
        inp.mdotmax = n.f64("mdotmax", inp.mdotmax);
        inp.sigcut = n.f64("sigcut", inp.sigcut);
    }
    if let Some(b) = blocks.get("emisdata") {
        let n = Nml::new(b);
        inp.ename = n.string("ename", &inp.ename);
        inp.mbh = n.f64("mbh", inp.mbh);
        inp.nfreq = n.i64("nfreq", inp.nfreq);
        inp.fmin = n.f64("fmin", inp.fmin);
        inp.fmax = n.f64("fmax", inp.fmax);
        inp.muval = n.f64("muval", inp.muval);
        inp.gmin = n.f64("gmin", inp.gmin);
        inp.gmax = n.f64("gmax", inp.gmax);
        inp.p1 = n.f64("p1", inp.p1);
        inp.p2 = n.f64("p2", inp.p2);
        inp.jetalpha = n.f64("jetalpha", inp.jetalpha);
        inp.stype = n.string("stype", &inp.stype);
        inp.delta = n.f64("delta", inp.delta);
        inp.nweights = n.i64("nweights", inp.nweights);
        inp.coefindx = n.i64_array("coefindx", inp.coefindx);
    }
    if let Some(b) = blocks.get("general") {
        let n = Nml::new(b);
        inp.use_geokerr = n.bool("use_geokerr", inp.use_geokerr);
        inp.nvals = n.i64("nvals", inp.nvals);
        inp.iname = n.string("iname", &inp.iname);
        inp.cflag = n.i64("cflag", inp.cflag as i64) as i32;
    }
    if let Some(b) = blocks.get("harm") {
        let n = Nml::new(b);
        inp.fdfile = n.string("fdfile", &inp.fdfile);
        inp.fhfile = n.string("fhfile", &inp.fhfile);
        inp.fgfile = n.string("fgfile", &inp.fgfile);
        inp.fsim = n.string("fsim", &inp.fsim);
        inp.fnt = n.i64("fnt", inp.fnt);
        inp.findf = n.i64("findf", inp.findf);
        inp.fnfiles = n.i64("fnfiles", inp.fnfiles);
        inp.fjonfix = n.i64("fjonfix", inp.fjonfix);
        inp.foffset = n.i64("foffset", inp.foffset);
        inp.fdindf = n.i64("fdindf", inp.fdindf);
        inp.fmagcrit = n.i64("fmagcrit", inp.fmagcrit);
        inp.fscalefac = n.f64("fscalefac", inp.fscalefac);
    }
    if let Some(b) = blocks.get("analytic") {
        let n = Nml::new(b);
        inp.fnw = n.i64("fnw", inp.fnw);
        inp.fwmin = n.f64("fwmin", inp.fwmin);
        inp.fwmax = n.f64("fwmax", inp.fwmax);
        inp.fnfreq_tab = n.i64("fnfreq_tab", inp.fnfreq_tab);
        inp.ffmin = n.f64("ffmin", inp.ffmin);
        inp.ffmax = n.f64("ffmax", inp.ffmax);
        inp.frmax = n.f64("frmax", inp.frmax);
        inp.fnr = n.i64("fnr", inp.fnr);
        inp.fsigt = n.f64("fsigt", inp.fsigt);
        inp.ffcol = n.f64("ffcol", inp.ffcol);
        inp.frspot = n.f64("frspot", inp.frspot);
        inp.fr0spot = n.f64("fr0spot", inp.fr0spot);
        inp.fn0spot = n.f64("fn0spot", inp.fn0spot);
        inp.ftscl = n.f64("ftscl", inp.ftscl);
        inp.frscl = n.f64("frscl", inp.frscl);
        inp.fmdot = n.f64("fmdot", inp.fmdot);
        inp.fnscl = n.f64("fnscl", inp.fnscl);
        inp.fnnthscl = n.f64("fnnthscl", inp.fnnthscl);
        inp.fnnthp = n.f64("fnnthp", inp.fnnthp);
        inp.fbeta = n.f64("fbeta", inp.fbeta);
        inp.fbl06 = n.i64("fbl06", inp.fbl06);
        inp.fnp = n.f64("fnp", inp.fnp);
        inp.ftp = n.f64("ftp", inp.ftp);
        inp.frin = n.f64("frin", inp.frin);
        inp.frout = n.f64("frout", inp.frout);
        inp.fthin = n.f64("fthin", inp.fthin);
        inp.fthout = n.f64("fthout", inp.fthout);
        inp.fphiin = n.f64("fphiin", inp.fphiin);
        inp.fphiout = n.f64("fphiout", inp.fphiout);
    }
    // frequency and mdot grids (upstream read_inputs lines 71-86)
    inp
}

/// Frequency grid (upstream `read_inputs` log spacing).
pub fn freqs(inp: &Inputs) -> Vec<f64> {
    let n = inp.nfreq as usize;
    if n <= 1 {
        return vec![inp.fmin];
    }
    (0..n)
        .map(|i| inp.fmin * (((i as f64) * (inp.fmax / inp.fmin).ln()) / ((n - 1) as f64)).exp())
        .collect()
}

/// Mdot grid.
pub fn mdots(inp: &Inputs) -> Vec<f64> {
    let n = inp.nmdot as usize;
    if n <= 1 {
        return vec![inp.mdotmin];
    }
    (0..n)
        .map(|i| {
            inp.mdotmin * (((i as f64) * (inp.mdotmax / inp.mdotmin).ln()) / ((n - 1) as f64)).exp()
        })
        .collect()
}

/// Observer mu grid.
pub fn mus(inp: &Inputs) -> Vec<f64> {
    let n = inp.nmu as usize;
    if n <= 1 {
        return vec![inp.mumin];
    }
    (0..n)
        .map(|i| inp.mumin + (inp.mumax - inp.mumin) / ((n - 1) as f64) * (i as f64))
        .collect()
}
