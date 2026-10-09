"""Python API smoke test for grtrans-rs.

Run after building/installing the wheel (maturin develop or pip install of a
built wheel):

    python tests/python/test_python_api.py

Compares a THINDISK image against the regenerated upstream reference.
"""
import sys
from pathlib import Path

import numpy as np

import grtrans_python as gt


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    print("grtrans_python version:", gt.version())

    # single geodesic ray
    x, k, lam = gt.geodesic_ray(
        alpha=-10.4, beta=-9.75, mu0=0.6428, spin=0.9375, standard=1, uout=0.04, nup=25
    )
    assert x.shape == (27, 4) and k.shape == (27, 4) and lam.shape == (27,)
    assert np.all(np.isfinite(x)) and np.all(np.isfinite(k))

    # THINDISK image vs upstream reference
    ab, ivals, nu = gt.run_image(
        fname="THINDISK", nfreq=25, nmu=1, fmin=2.41e16, fmax=6.31e18,
        ename="BBPOL", nvals=4, spin=0.9, standard=2, nn=[100, 100, 1],
        uout=0.01, mbh=10, mumin=0.26, mumax=0.26, gridvals=[-21, 21, -21, 21],
    )
    assert ab.shape == (10000, 2)
    assert ivals.shape == (10000, 4, 25)
    assert nu.shape == (25,)

    ref = np.fromfile(
        root / "reference/fixtures/thindisk/ivals.f64.bin", dtype="<f8"
    ).reshape(10000, 4, 25)
    err = np.abs(ivals - ref).sum() / np.abs(ref).sum()
    print(f"THINDISK relative error via Python API: {err:.3e}")
    if err >= 1e-2:
        print("FAILED: error above upstream tolerance")
        return 1
    print("OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
