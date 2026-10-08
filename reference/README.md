# Reference fixtures

This directory contains reference outputs generated with the **upstream
Fortran GRTRANS** (`jadexter/grtrans` @ `c76cb11`) and used to validate the
Rust port.

```
reference/
├── upstream/                  # everything needed to rebuild upstream
│   ├── patch_namelists.py     # gfortran>=10 compatibility patch
│   ├── gfortran15.patch       # diff produced by the patch
│   └── Makefile.top           # build configuration used
└── fixtures/
    ├── generation_manifest.json
    ├── thindisk/ ...
    ├── ffjet/ ...
    ├── ffjet_delo/ ...
    ├── ffjet_formal/ ...
    ├── ffjet_unpol/ ...
    ├── harm/ ...
    ├── sphacc/ ...
    └── powerlaw/
```

Each case directory contains the input decks (`inputs.in`, `files.in`), the
arrays as `.npy` and raw little-endian float64 `.f64.bin`, and a
`manifest.json` with provenance (upstream commit, compiler and library
versions, input parameters), SHA-256 checksums, and the measured relative
error against the upstream-shipped regression pickle.

Regeneration is documented in `docs/REFERENCE_RESULTS.md` and implemented in
`scripts/gen_reference.py` (requires a built upstream binary and the f2py
module; see `scripts/build_upstream.sh`).

Integrity is checked with:

```bash
python3 scripts/verify_fixtures.py
```

Do not edit fixture files by hand; regenerate them.
