"""Run matching Fortran and Rust THINDISK decks and plot residuals.

Usage: python scripts/compare_thindisk_inclinations.py FORTRAN_BINARY RUST_BINARY OUTPUT
Requires numpy, astropy and matplotlib. Use a new output directory.
The upstream ch24_vals.txt must be beside the Fortran binary.
"""
import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess

from astropy.io import fits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np


FIXTURE = Path(__file__).resolve().parents[1] / "reference/fixtures/thindisk"
COSINES = (0.98, 0.8, 0.6, 0.4, 0.26, 0.05)


def run_solver(binary, folder, inputs, table):
    folder.mkdir(parents=True)
    (folder / "inputs.in").write_text(inputs)
    (folder / "files.in").write_text("&files\n ifile='inputs.in',\n ofile='disk.fits',\n/\n")
    shutil.copyfile(table, folder / "ch24_vals.txt")
    with (folder / "run.log").open("w") as log:
        subprocess.run([binary], cwd=folder, check=True, stdout=log,
                       stderr=subprocess.STDOUT, env={**os.environ, "OMP_NUM_THREADS": "1"})
    with fits.open(folder / "disk.fits") as hdus:
        if len(hdus) != 26:
            raise ValueError(f"{folder}: expected the 25 THINDISK frequency images")
        pixels = hdus[0].data.reshape(10000, 2).copy()
        frequencies = np.array([h.header["NU"] for h in hdus[1:]])
        for h in hdus[1:]:
            if (h.header["NX"], h.header["NY"]) != (100, 100):
                raise ValueError(f"{folder}: unexpected camera dimensions")
        cube = np.stack([h.data.reshape(10000, 4) for h in hdus[1:]]).astype(float)
    if not np.isfinite(cube).all():
        raise ValueError(f"{folder}: nonfinite Stokes data")
    return pixels, frequencies, cube


def compare(mu, template, fortran, rust, out):
    inputs = template
    for field in ("mumin", "mumax"):
        original = f"{field}=0.26,"
        if inputs.count(original) != 1:
            raise ValueError(f"Unexpected THINDISK input for {field}")
        inputs = inputs.replace(original, f"{field}={mu},")
    folder = out / f"mu-{mu}"
    table = fortran.parent / "ch24_vals.txt"
    pixels, nu, ref = run_solver(fortran, folder / "fortran", inputs, table)
    coords, freq, trial = run_solver(rust, folder / "rust", inputs, table)
    np.testing.assert_array_equal(pixels, coords)
    np.testing.assert_allclose(freq, nu, rtol=1e-6, atol=0)
    diff = trial - ref
    peak = np.abs(ref[0, :, 0]).max()
    if not peak > 0:
        raise ValueError(f"{folder}: empty reference intensity image")
    # Sort by beta (rows), then alpha (columns), using the FITS camera coordinates.
    order = np.lexsort((pixels[:, 0], pixels[:, 1]))
    residual = (diff[0, order, 0] / peak).reshape(100, 100)
    row = dict(mu=mu, inclination_deg=float(np.degrees(np.arccos(mu))),
               all_stokes_relative_L1=np.abs(diff).sum() / np.abs(ref).sum(),
               map_frequency_Hz=nu[0], max_abs_I_residual_over_peak=np.abs(residual).max())
    return row, residual


def plot(rows, maps, out):
    limit = max(np.max(np.abs(a)) for a in maps)
    # Common sign-preserving scale; keep a finite display range for exact agreement.
    limit = max(limit, np.finfo(float).eps)
    norm = SymLogNorm(linthresh=limit * 1e-3, vmin=-limit, vmax=limit)
    fig, axes = plt.subplots(2, 3, figsize=(9, 6.5), layout="constrained")
    for ax, row, data in zip(axes.flat, rows, maps):
        im = ax.imshow(data, origin="lower", cmap="RdBu_r", norm=norm, interpolation="nearest")
        ax.set_title(f"{row['inclination_deg']:.1f}°")
        ax.label_outer()
    fig.supxlabel("Camera column (increasing α)")
    fig.supylabel("Camera row (increasing β)")
    fig.suptitle(f"Thin disk · Stokes I residuals at {rows[0]['map_frequency_Hz']:.3g} Hz")
    fig.colorbar(im, ax=axes, label="(Rust − Fortran) / peak Fortran I")
    fig.savefig(out / "residual-maps.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    ax.plot([r["inclination_deg"] for r in rows], [r["all_stokes_relative_L1"] for r in rows], "o-")
    ax.set(xlabel="Inclination from disk axis (degrees)", ylabel="Relative L1 error",
           title="Thin disk · all four Stokes components and 25 frequencies")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    fig.savefig(out / "inclination-errors.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("fortran", "rust", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    template = (FIXTURE / "inputs.in").read_text()
    rows, maps = [], []
    for mu in COSINES:
        row, residual = compare(mu, template, args.fortran.resolve(), args.rust.resolve(), args.output)
        rows.append(row)
        maps.append(residual)
    with (args.output / "inclinations.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    plot(rows, maps, args.output)


if __name__ == "__main__":
    main()
