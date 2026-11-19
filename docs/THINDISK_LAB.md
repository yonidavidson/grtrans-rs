# A first polarized image with GRTRANS in Rust

How does viewing angle change the apparent shape of a thin accretion disk? Start by reproducing a reference image. Then change only the viewing angle and compare the same frequency in both runs.

## Reproduce the reference

On Linux or macOS, start with Git and a stable Rust toolchain. The first build needs internet access to download dependencies. This exercise uses a public revision and the input data bundled with it.

```sh
git clone https://github.com/yonidavidson/grtrans-rs.git
cd grtrans-rs
git checkout --detach 87e726ef9c1c475d65d5bcc5a53721081f6da101
cargo test --locked --release -p grtrans --test fits_output -- --nocapture
```

Look for two passing tests and the line giving the THINDISK relative error. The test sums absolute differences between Rust and Fortran over 10,000 pixels, four Stokes components, and 25 frequencies, then divides by the sum of the absolute Fortran values. A Linux run gave 1.14 × 10⁻⁷. The FITS test requires a value below 10⁻⁴. This comparison checks numerical agreement for this particular disk model.

## Change the viewing angle

Run the following from the repository root. Each directory gets its own input deck and FITS file. Keeping them together lets someone reproduce either result.

```sh
cargo build --locked --release -p grtrans
bin="$PWD/target/release/grtrans"
lab=$(mktemp -d)
for view in reference face_onward; do
  mkdir "$lab/$view"
  cp reference/fixtures/thindisk/inputs.in "$lab/$view/inputs.in"
  printf "&files\n ifile='inputs.in',\n ofile='disk.fits',\n/\n" > "$lab/$view/files.in"
done
sed 's/mumin=0.26/mumin=0.8/; s/mumax=0.26/mumax=0.8/' \
  "$lab/reference/inputs.in" > "$lab/face_onward/inputs.in"
for view in reference face_onward; do
  (cd "$lab/$view" && "$bin")
done
printf 'Results in %s\n' "$lab"
```

The angle is measured from the disk axis. Its cosine changes from 0.26 to 0.8, corresponding to about 75° and 37°. Keep spin, frequency, and color scale fixed when comparing the images.

The primary FITS array stores the camera coordinates. Each subsequent extension holds one frequency, recorded in the `NU` header. With Astropy installed, the following reads the total intensity from the first frequency and reports its peak for each view.

```sh
python - "$lab" <<'PY'
from pathlib import Path
import sys
from astropy.io import fits

for path in sorted(Path(sys.argv[1]).glob("*/disk.fits")):
    with fits.open(path) as hdus:
        image = hdus[1]
        stokes = image.data.reshape(-1, 4)
        print(path.parent.name, image.header["NU"], stokes[:, 0].max())
PY
```

The columns are I, Q, U, and V. For a lab report, record the commit, both inclinations, the chosen frequency, and the reference error. Explain why changing the inclination should change the image even when the transfer code is working correctly.

The original Fortran code and numerical methods are by Jason Dexter and collaborators. Credit the original work when using this exercise. [Dexter (2016)](https://doi.org/10.1093/mnras/stw1526) · [Dexter & Agol (2009)](https://doi.org/10.1088/0004-637X/696/2/1616) · [Original GRTRANS repository](https://github.com/jadexter/grtrans)
