#!/usr/bin/env python3
"""Generate reference fixtures with the upstream Fortran GRTRANS binary.

This script reproduces the runs in run_grtrans_test_problems_public.py of the
upstream repository (jadexter/grtrans @ c76cb11) using the compiled Fortran
binary through the original grtrans_batch.py driver, and stores the resulting
arrays as .npy files plus a JSON manifest with checksums.

Usage:
    GRTRANS_UPSTREAM=/path/to/upstream \
    python3 scripts/gen_reference.py --outdir reference/fixtures [--only NAME]

It must be executed with a Python environment that can import grtrans_batch
(numpy, astropy, pgrtrans f2py module and ./grtrans binary in the upstream
directory).  See docs/REFERENCE_RESULTS.md for the build recipe.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_pickle(path: Path):
    with open(path, "rb") as f:
        return pickle.load(f, encoding="latin1")


def upstream_commit(upstream: Path) -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=upstream, capture_output=True, text=True
    ).stdout.strip()


def main() -> int:
    up = Path(os.environ.get("GRTRANS_UPSTREAM", os.getcwd())).resolve()
    os.chdir(up)
    sys.path.insert(0, str(up))

    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="reference/fixtures")
    ap.add_argument("--only", default=None, help="comma-separated subset")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    outroot = (up / args.outdir).resolve()
    outroot.mkdir(parents=True, exist_ok=True)
    only = set(args.only.split(",")) if args.only else None

    want = lambda name: only is None or name in only

    import grtrans_batch as gr  # noqa: E402  (requires compiled pgrtrans)

    def run_case(name, kwargs, save_keys=("ivals", "ab", "nu", "spec", "lp", "cp"),
                 ref_pickle=None, ref_transform=None, extra_refs=None):
        if not want(name):
            return None
        d = outroot / name
        if (d / "manifest.json").exists() and not args.overwrite:
            print(f"[{name}] exists, skipping (use --overwrite)")
            return json.loads((d / "manifest.json").read_text())
        d.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        x = gr.grtrans()
        x.write_grtrans_inputs("inputs.in", **kwargs)
        x.run_grtrans()
        x.read_grtrans_output()
        wall = time.time() - t0
        manifest = {
            "name": name,
            "upstream_commit": upstream_commit(up),
            "upstream_worktree_dirty": bool(
                subprocess.run(
                    ["git", "status", "--porcelain"], cwd=up, capture_output=True, text=True
                ).stdout.strip()
            ),
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "platform": platform.platform(),
            "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
            "wall_seconds": round(wall, 2),
            "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "inputs": {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                       for k, v in kwargs.items()},
            "files": {},
        }
        # copy exact input decks used
        for fn in ("inputs.in", "files.in"):
            src = up / fn
            if src.exists():
                (d / fn).write_bytes(src.read_bytes())
        # save arrays (npy for Python users, raw LE float64 for Rust tests)
        for key in save_keys:
            arr = getattr(x, key, None)
            if arr is None:
                continue
            arr = np.asarray(arr, dtype="<f8")
            if arr.size == 0:
                continue
            fn = f"{key}.npy"
            np.save(d / fn, arr)
            raw = d / f"{key}.f64.bin"
            raw.write_bytes(arr.tobytes(order="C"))
            manifest["files"][fn] = {
                "sha256": sha256(d / fn),
                "raw": f"{key}.f64.bin",
                "raw_sha256": sha256(raw),
                "raw_dtype": "little-endian float64",
                "shape": list(arr.shape),
                "dtype": str(arr.dtype),
            }
        # comparison against shipped upstream pickle
        if ref_pickle is not None:
            ref = np.asarray(load_pickle(up / ref_pickle))
            if ref_transform is None:
                cur = np.asarray(x.ivals, dtype="<f8")
            else:
                cur = ref_transform(x)
            cur = np.asarray(cur, dtype="<f8")
            if cur.shape != ref.shape:
                print(f"[{name}] WARNING shape mismatch {cur.shape} vs {ref.shape}")
            diff = np.sum(np.abs(cur - ref)) / np.sum(np.abs(ref))
            manifest["upstream_pickle"] = ref_pickle
            manifest["rel_err_vs_upstream_pickle"] = float(diff)
            print(f"[{name}] rel err vs {ref_pickle}: {diff:.3e}  ({wall:.1f}s)")
        (d / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        return manifest

    results = {}

    # --- SPHACC 1D intensity profile and spectrum -------------------------
    results["sphacc"] = run_case(
        "sphacc",
        dict(fname="SPHACC", nfreq=25, nmu=1, fmin=1e8, fmax=1e15,
             ename="SYNCHTHAV", nvals=1, spin=0.0, mbh=1.0, standard=1,
             nn=[10000, 1, 100], gridvals=[0.0, 400.0, 0.0, 0.0],
             uout=0.0025, oname="sphacc_abs.out"),
        ref_pickle="test_grtrans_sphacc_spectrum.p",
        ref_transform=lambda x: np.asarray(x.spec),
    )

    # --- FFJET ------------------------------------------------------------
    results["ffjet"] = run_case(
        "ffjet",
        dict(fname="FFJET", jdfile="m87bl09rfp10xi5a998fluidvars.bin", nfreq=1,
             nmu=1, fmin=3.45e11, fmax=3.45e11, ename="POLSYNCHPL", nvals=4,
             spin=0.998, standard=1, nn=[100, 100, 400], uout=0.01,
             mbh=3.4e9, mumin=0.906, mumax=0.906,
             gridvals=[-40, 20, -20, 40], ntscl=2.0, nrscl=70.0),
        ref_pickle="test_grtrans_ffjet.p",
    )

    results["ffjet_delo"] = run_case(
        "ffjet_delo",
        dict(fname="FFJET", jdfile="m87bl09rfp10xi5a998fluidvars.bin", nfreq=1,
             nmu=1, fmin=3.45e11, fmax=3.45e11, ename="POLSYNCHPL", nvals=4,
             spin=0.998, standard=1, nn=[100, 100, 1600], uout=0.01,
             mbh=3.4e9, mumin=0.906, mumax=0.906,
             gridvals=[-40, 20, -20, 40], iname="delo", ntscl=2.0, nrscl=70.0),
        save_keys=("ivals", "ab", "nu", "spec", "lp", "cp"),
    )

    results["ffjet_formal"] = run_case(
        "ffjet_formal",
        dict(fname="FFJET", jdfile="m87bl09rfp10xi5a998fluidvars.bin", nfreq=1,
             nmu=1, fmin=3.45e11, fmax=3.45e11, ename="POLSYNCHPL", nvals=4,
             spin=0.998, standard=1, nn=[100, 100, 1600], uout=0.01,
             mbh=3.4e9, mumin=0.906, mumax=0.906,
             gridvals=[-40, 20, -20, 40], iname="formal", ntscl=2.0, nrscl=70.0),
        save_keys=("ivals", "ab", "nu", "spec", "lp", "cp"),
    )

    results["ffjet_unpol"] = run_case(
        "ffjet_unpol",
        dict(fname="FFJET", jdfile="m87bl09rfp10xi5a998fluidvars.bin", nfreq=1,
             nmu=1, fmin=3.45e11, fmax=3.45e11, ename="SYNCHPL", nvals=1,
             spin=0.998, standard=1, nn=[100, 100, 400], uout=0.01,
             mbh=3.4e9, mumin=0.906, mumax=0.906,
             gridvals=[-40, 20, -20, 40], ntscl=2.0, nrscl=70.0),
    )

    # --- THINDISK ---------------------------------------------------------
    results["thindisk"] = run_case(
        "thindisk",
        dict(fname="THINDISK", nfreq=25, nmu=1, fmin=2.41e16, fmax=6.31e18,
             ename="BBPOL", nvals=4, spin=0.9, standard=2, nn=[100, 100, 1],
             uout=0.01, mbh=10, mumin=0.26, mumax=0.26,
             gridvals=[-21, 21, -21, 21]),
        ref_pickle="test_grtrans_thindisk.p",
    )

    # --- HARM -------------------------------------------------------------
    results["harm"] = run_case(
        "harm",
        dict(fname="HARM", nfreq=1, nmu=1, fmin=2.3e11, fmax=2.3e11,
             ename="POLSYNCHTH", nvals=1, spin=0.9375, standard=1,
             nn=[150, 150, 400], uout=0.04, mbh=4e6, mdotmin=1.57e15,
             mdotmax=1.57e15, nmdot=1, mumin=0.6428, mumax=0.6428,
             gridvals=[-13, 13, -13, 13], hhfile="dump040", hdfile="dump",
             hindf=40, hnt=1, muval=1.0 / 4.0, gmin=1.0),
        ref_pickle="test_grtrans_harm.p",
    )

    # --- POWERLAW toroidal field -----------------------------------------
    results["powerlaw"] = run_case(
        "powerlaw",
        dict(fname="POWERLAW", nfreq=1, nmu=1, fmin=3.45e11, fmax=3.45e11,
             ename="POLSYNCHTH", nvals=4, spin=0.0, standard=1,
             nn=[200, 200, 1600], uout=0.00005, mbh=4e6, mumin=0.5, mumax=0.5,
             nrotype=1, gridvals=[1200.0, 4000.0, 0.0, 2.0 * np.pi],
             iname="lsoda", srin=3200.0, srout=3300.0, ntscl=5e11,
             sthin=-0.02, sthout=0.02, rcut=4000.0, snscl=1e5, phi0=-0.5,
             sphiin=0.0, gmin=1.0),
        ref_pickle="test_toroidalfield.p",
    )

    (outroot / "generation_manifest.json").write_text(
        json.dumps({k: v for k, v in results.items() if v is not None}, indent=2) + "\n"
    )
    print("done:", ", ".join(k for k, v in results.items() if v))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
