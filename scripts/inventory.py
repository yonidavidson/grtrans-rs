#!/usr/bin/env python3
"""Generate docs/SOURCE_INVENTORY.md from an upstream grtrans checkout.

Usage:
    python3 scripts/inventory.py --upstream /path/to/upstream [--out docs/SOURCE_INVENTORY.md]

Line counts are computed live from the checkout so the document is always
consistent with the recorded revision (see docs/UPSTREAM.md).
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

# Descriptive metadata for every upstream source file.
# fields: responsibility | numerical methods / physical assumptions |
#         Rust target crate::module | dependencies
META = {
    # ---- core infrastructure -------------------------------------------------
    "phys_constants.f90": ("Physical constants (CGS)", "h,k,c,e,G,m,mp,msun,sigma_T,sigma_B", "grtrans-core::constants", "-"),
    "math.f90": ("Cumulative trapezoid sums, dot products, Brent root finding", "Brent (zbrent), cumulative trapezoid TSUM", "grtrans-core::math", "-"),
    "interpolate.f90": ("Multilinear interpolation and bracketing helpers", "bilinear/trilinear/quadrilinear interpolation; hunt/locate", "grtrans-core::interpolate", "hunt.f, locate.f, polyvl.f"),
    "class_four_vector.f90": ("Four-vector algebra with stored metric", "10-component symmetric metric storage; raise/lower; inner products", "grtrans-core::four_vector", "-"),
    "kerr.f90": ("Kerr metric, frames, geodesic constants, polarization transport", "BL covariant/contravariant metric, KS metric, LNRF transforms, uks2ubl, bl2ks, calcg, calc_nullp, Walker-Penrose transport, comoving orthonormal frame, Maxwell-Juttner utilities", "grtrans-core::kerr", "class_four_vector, math, phys_constants"),
    "bessel.f90": ("Incomplete modified Bessel functions", "Series/continued-fraction evaluation (NR style)", "grtrans-core::bessel", "-"),
    "chandra_tab24.f90": ("Chandrasekhar (1960) table 24 polarization loading/interpolation", "Table interpolation of degree of polarization vs viewing angle", "grtrans-core::chandra_tab24", "ch24_vals.txt"),
    "calcgmin.f90": ("Minimum Lorentz factor helper for power-law synchrotron", "Root solving", "grtrans-physics::emis", "polsynchemis"),
    "calc_maxjutt.f90": ("Bessel-function evaluation maxima (Maxwell-Juttner)", "Numerical maxima of Bessel combinations", "grtrans-physics::emis", "polsynchemis"),
    "calc_maxcomp.f90": ("Bessel-function evaluation maxima (power-law)", "Numerical maxima of Bessel combinations", "grtrans-physics::emis", "polsynchemis"),
    "geokerr_wrapper.f": ("geokerr Kerr geodesic solver (Dexter & Agol 2009)", "Semi-analytic geodesics via elliptic integrals and Jacobi elliptic functions; camera initialization; turning-point bookkeeping", "grtrans-geodesics::geokerr", "SNCNDN, ZROOTS, LAGUER (self-contained)"),
    "hunt.f": ("Monotonic-array bracketing (NR-style)", "Bisection from an initial guess", "grtrans-core::interpolate", "-"),
    "locate.f": ("Monotonic-array bracketing (NR-style)", "Bisection", "grtrans-core::interpolate", "-"),
    "polint.f": ("Polynomial interpolation", "Neville's algorithm", "grtrans-core::bessel", "-"),
    "polyvl.f": ("Polynomial value/derivative evaluation", "Horner scheme", "grtrans-core::bessel", "-"),
    # ---- radiative transfer --------------------------------------------------
    "rad_trans.f90": ("Radiation-transfer object and integrator dispatch", "iname -> iflag mapping", "grtrans-transfer::rad_trans", "grtrans_inputs"),
    "radtrans_integrate.f90": ("Polarized transfer equation integration", "Del Zanna & Bucciantini (2002) scheme; formal solution; RHS/Jacobian for LSODA; optical-depth quadrature; 4x4 matrix inversion", "grtrans-transfer::integrate", "odepack, interpolate, math"),
    "odepack.f90": ("LSODA wrapper", "calls ODEPACK LSODA (F77)", "grtrans-transfer::lsoda", "opkda1.f, opkda2.f, opkdmain.f"),
    "odepack_aux.f": ("LSODA auxiliary glue (F77)", "legacy interface", "grtrans-transfer::lsoda", "-"),
    "opkda1.f": ("ODEPACK driver/stepper (F77)", "LSODA/LSODE (Hindmarsh, Petzold)", "grtrans-transfer::lsoda", "-"),
    "opkda2.f": ("ODEPACK linear solvers (F77)", "dense/banded LU; nonlinear iteration", "grtrans-transfer::lsoda", "-"),
    "opkdmain.f": ("ODEPACK core solver (F77)", "LSODA automatic stiff/nonstiff switching", "grtrans-transfer::lsoda", "-"),
    # ---- emissivity ----------------------------------------------------------
    "emis.f90": ("Emissivity dispatch and blackbody/table models", "BB/FBB, INTERP, RHO, hybrid composition; relativistic invariants; Faraday rotation of coefficients", "grtrans-physics::emissivity", "polsynchemis, math, chandra_tab24, calc_maxjutt, calc_maxcomp"),
    "polsynchemis.f90": ("Polarized synchrotron emission/absorption", "Melrose (1983)/Dexter (2011) coefficients; thermal (Maxwell-Juttner) and power-law distributions; background subtraction (MAXCOMP/MAXJUTT)", "grtrans-physics::polsynch", "phys_constants, bessel"),
    "interpolate_aux.f": ("f2py glue", "-", "- (driver only)", "-"),
    # ---- fluid models --------------------------------------------------------
    "fluid.f90": ("Fluid-model dispatch and geometry of the emitting region", "model selection, coordinate transforms, cgs conversion, source parameters per model", "grtrans-physics::fluid", "class_four_vector, phys_constants, interpolate, kerr, all fluid_model_*"),
    "fluid_model_thindisk.f90": ("NT73/Page-Thorne thin disk", "analytic NT73 solution, Krolik flux function", "grtrans-physics::models::thindisk", "kerr"),
    "fluid_model_phatdisk.f90": ("Dexter & Agol (2011) inhomogeneous disk", "log-normal temperature table, blackbody convolution", "grtrans-physics::models::phatdisk", "fluid_model_thindisk, math, interpolate"),
    "fluid_model_numdisk.f90": ("Numerical thin-disk temperature table", "read/interpolate T_eff(r,phi)", "grtrans-physics::models::numdisk", "fluid_model_thindisk, interpolate"),
    "fluid_model_hotspot.f90": ("Orbiting hot spot (Broderick & Loeb 2005/2006)", "Keplerian spot, power-law emissivity", "grtrans-physics::models::hotspot", "kerr, fluid_model_thindisk"),
    "fluid_model_hotspot_schnittman.f90": ("Orbiting hot spot (Schnittman & Bertschinger 2004)", "Keplerian spot", "grtrans-physics::models::hotspot_schnittman", "kerr"),
    "fluid_model_sphacc.f90": ("Bondi/Michel spherical accretion", "semi-analytic GR accretion solution with numerical radial table", "grtrans-physics::models::sphacc", "interpolate, kerr"),
    "fluid_model_ffjet.f90": ("Broderick & Loeb (2009) jet", "semi-analytic jet on a conical grid", "grtrans-physics::models::ffjet", "interpolate, kerr"),
    "fluid_model_powerlaw.f90": ("Power-law test flow", "analytic density/temperature profiles", "-", "phys_constants"),
    "fluid_model_sariaf.f90": ("Sari & Esin/self-similar analytic flow", "analytic profiles", "-", "phys_constants"),
    "fluid_model_toy.f90": ("Toy test model", "analytic", "grtrans-physics::models::toy", "phys_constants"),
    "fluid_model_harm.f90": ("HARM GRMHD 2D data reader", "binary dump IO; KS<->BL; LNRF velocity reconstruction", "grtrans-physics::models::harm", "class_four_vector, interpolate, kerr, math"),
    "fluid_model_harm3d.f90": ("HARM3D GRMHD 3D data reader", "binary dump IO, trilinear interpolation", "-", "class_four_vector, interpolate, kerr, math"),
    "fluid_model_harmpi.f90": ("HARM-pi GRMHD data reader", "binary dump IO, primitive recovery", "-", "class_four_vector, interpolate, kerr, math"),
    "fluid_model_iharm.f90": ("iharm GRMHD data reader", "binary dump IO", "-", "class_four_vector, interpolate, kerr, math"),
    "fluid_model_koral.f90": ("KORAL GRMHD 2D data reader", "binary dump IO", "-", "class_four_vector, interpolate, kerr, math"),
    "fluid_model_koral3d.f90": ("KORAL3D GRMHD data reader", "binary dump IO", "-", "class_four_vector, interpolate, kerr, math"),
    "fluid_model_mb09.f90": ("McKinney & Blandford (2009) GRMHD data reader", "binary dump IO", "-", "class_four_vector, interpolate, kerr, mat"),
    "fluid_model_thickdisk.f90": ("Thick-disk GRMHD composite model", "binary dump IO; torus construction", "-", "class_four_vector, interpolate, kerr, math"),
    # ---- driver / IO / camera ------------------------------------------------
    "geodesics.f90": ("Geodesic sampling and geokerr orchestration", "per-pixel geokerr call; BL coordinate assembly; wave-vector reconstruction", "grtrans-geodesics::geodesics", "class_four_vector, interpolate, kerr, geokerr_wrapper.f"),
    "class_geokerr.f90": ("geokerr f2py wrapper object", "state for f2py driver", "grtrans-python (validation only)", "geokerr_wrapper.f"),
    "camera.f90": ("Camera/pixel geometry", "polar (nrotype=1) and Cartesian (nrotype=2) pixel grids", "grtrans-geodesics::camera", "fits.f90"),
    "read_inputs.f90": ("Input namelist parsing", "Fortran namelists geodata/fluiddata/emisdata/general/harm/analytic", "grtrans-io::inputs", "-"),
    "grtrans_driver.f90": ("Per-ray driver: fluid, emissivity, transfer, output", "invariant scaling j*g^2, K/g; optical depth; LSODA/delo/formal dispatch; debug quantities", "grtrans-driver (crate `grtrans`)", "all of the above"),
    "pgrtrans.f90": ("f2py entry point / run orchestration", "loops over mu, calls grtrans_driver with OpenMP", "grtrans-python", "grtrans_driver"),
    "grtrans.f90": ("Input-file entry point", "read inputs, run, save", "grtrans CLI", "pgrtrans, read_inputs"),
    "grtrans_program.f90": ("Main program reading files.in", "namelist files: ifile/ofile", "grtrans CLI", "grtrans.f90"),
    "fits.f90": ("FITS output via cfitsio", "images + header keys", "grtrans-io::fits", "cfitsio"),
    "fits_prealwin_07102014.f90": ("Historical FITS variant (unused)", "-", "-", "cfitsio"),
    "emis_old.f90": ("Historical emissivity variant (unused)", "-", "-", "-"),
    # ---- tests (Fortran) -----------------------------------------------------
    "test_kerr.f90": ("Kerr-metric unit checks", "prints metric values along rays", "grtrans-core tests (fixtures)", "-"),
    "test_harm.f90": ("HARM reader check", "prints fluid values", "grtrans-physics tests (fixtures)", "fluid_model_harm"),
    "test_ffjet.f90": ("FFJET reader check", "prints fluid values", "grtrans-physics tests (fixtures)", "fluid_model_ffjet"),
    "test_hotspot.f90": ("Hotspot model check", "prints fluid values", "grtrans-physics tests (fixtures)", "fluid_model_hotspot"),
    # ---- Python --------------------------------------------------------------
    "grtrans_batch.py": ("Python driver: builds inputs, runs binary, reads FITS output", "-", "scripts/gen_reference.py equivalent", "numpy, astropy, pgrtrans"),
    "run_grtrans_test_problems_public.py": ("Regression driver for the six shipped fixtures", "-", "tests/reference", "numpy"),
    "unit_tests_public.py": ("Small standalone checks", "-", "tests", "numpy"),
    "unit_tests_integration.py": ("Integration checks of transfers", "-", "tests", "numpy"),
    "simple_radtrans_integrate_tests.py": ("Transfer-equation checks", "-", "tests", "numpy"),
    "analytic_pol_rad_trans.py": ("Analytic comparison solution for polarized transfer", "Degl'Innocenti-style analytic solution", "tests/reference", "numpy"),
    "namelist.py": ("Namelist writer used by the Python driver", "-", "grtrans-io::inputs (writer)", "-"),
    "geokerr_interface.py": ("Python wrapper of geokerr via f2py", "-", "validation scripts", "geokerr module"),
    "read_geodebug_file.py": ("geodebug.out reader (debug mode)", "-", "debug tooling", "numpy"),
    "pgriter.py": ("Parallel job submission helper (cluster-specific)", "-", "- (not ported)", "python"),
    "pgrface.py": ("Parallel job submission helper (cluster-specific)", "-", "- (not ported)", "python"),
    "ppslave.py": ("Parallel job submission helper (cluster-specific)", "-", "- (not ported)", "python"),
    "ray_integrate.py": ("Ray-integration development script", "-", "not ported (development script)", "numpy"),
    "synch_integrate.py": ("Synchrotron-integration development script", "-", "not ported (development script)", "numpy"),
    "emitter_observer.py": ("Emitter-observer development script", "-", "not ported (development script)", "numpy"),
    "quick_radtrans_integrate_test.py": ("Quick transfer test", "-", "not ported (development script)", "numpy"),
    "radtrans_integrate_delo_harm_test.py": ("Delo integrator development test", "-", "not ported (development script)", "numpy"),
    "commit_grtrans.py": ("Upstream commit helper", "-", "not ported", "-"),
}

CATEGORIES = [
    ("Core infrastructure", ["phys_constants.f90", "math.f90", "interpolate.f90", "class_four_vector.f90", "kerr.f90", "bessel.f90", "chandra_tab24.f90", "calcgmin.f90", "calc_maxjutt.f90", "calc_maxcomp.f90", "hunt.f", "locate.f", "polint.f", "polyvl.f"]),
    ("Geodesics", ["geokerr_wrapper.f", "geodesics.f90", "class_geokerr.f90", "camera.f90"]),
    ("Emissivity", ["emis.f90", "polsynchemis.f90", "emis_old.f90"]),
    ("Radiative transfer", ["rad_trans.f90", "radtrans_integrate.f90", "odepack.f90", "odepack_aux.f", "opkda1.f", "opkda2.f", "opkdmain.f"]),
    ("Fluid models", ["fluid.f90", "fluid_model_thindisk.f90", "fluid_model_phatdisk.f90", "fluid_model_numdisk.f90", "fluid_model_hotspot.f90", "fluid_model_hotspot_schnittman.f90", "fluid_model_sphacc.f90", "fluid_model_ffjet.f90", "fluid_model_powerlaw.f90", "fluid_model_sariaf.f90", "fluid_model_toy.f90", "fluid_model_harm.f90", "fluid_model_harm3d.f90", "fluid_model_harmpi.f90", "fluid_model_iharm.f90", "fluid_model_koral.f90", "fluid_model_koral3d.f90", "fluid_model_mb09.f90", "fluid_model_thickdisk.f90"]),
    ("Driver, IO, program", ["grtrans_driver.f90", "pgrtrans.f90", "grtrans.f90", "grtrans_program.f90", "read_inputs.f90", "fits.f90", "fits_prealwin_07102014.f90", "interpolate_aux.f"]),
    ("Fortran tests", ["test_kerr.f90", "test_harm.f90", "test_ffjet.f90", "test_hotspot.f90"]),
    ("Python drivers and tests", ["grtrans_batch.py", "run_grtrans_test_problems_public.py", "unit_tests_public.py", "unit_tests_integration.py", "simple_radtrans_integrate_tests.py", "analytic_pol_rad_trans.py", "namelist.py", "geokerr_interface.py", "read_geodebug_file.py", "pgriter.py", "pgrface.py", "ppslave.py", "ray_integrate.py", "synch_integrate.py", "emitter_observer.py", "quick_radtrans_integrate_test.py", "radtrans_integrate_delo_harm_test.py", "commit_grtrans.py"]),
]


def count_lines(path: Path) -> int:
    with open(path, "rb") as f:
        return sum(1 for _ in f)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", required=True)
    ap.add_argument("--out", default="docs/SOURCE_INVENTORY.md")
    args = ap.parse_args()
    up = Path(args.upstream)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=up,
                            capture_output=True, text=True).stdout.strip()

    lines = []
    lines.append("# GRTRANS source inventory\n")
    lines.append(f"Generated by `scripts/inventory.py` from upstream revision\n`{commit}`. "
                 "Line counts are computed from the checkout at that revision.\n")
    seen = set()
    total = 0
    unknown = []
    for cat, files in CATEGORIES:
        lines.append(f"\n## {cat}\n")
        lines.append("| File | Lang | LOC | Responsibility | Methods / assumptions | Rust target | Deps |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for fn in files:
            p = up / fn
            assert p.exists(), fn
            loc = count_lines(p)
            total += loc
            seen.add(fn)
            lang = "F90" if fn.endswith(".f90") else ("F77" if fn.endswith(".f") else "Py")
            resp, meth, target, deps = META.get(fn, ("?", "?", "?", "?"))
            lines.append(f"| `{fn}` | {lang} | {loc} | {resp} | {meth} | {target} | {deps} |")
    # any source files not categorized
    for p in sorted(list(up.glob("*.f90")) + list(up.glob("*.f")) + list(up.glob("*.py"))):
        if p.name not in seen and p.name not in ("smoke_thindisk.py",):
            unknown.append(p.name)
    if unknown:
        lines.append("\n## Uncategorized files\n")
        for fn in unknown:
            loc = count_lines(up / fn)
            total += loc
            lines.append(f"- `{fn}` ({loc} lines)")
    lines.append(f"\n**Total source lines (categorized Files): {total}**\n")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out} ({total} lines counted)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
