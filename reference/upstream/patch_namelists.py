#!/usr/bin/env python3
"""Move namelist statements after variable declarations for gfortran>=10.

Newer gfortran requires all namelist members to be declared before the
namelist statement.  The upstream code declares namelists immediately
after `implicit none`, before the type declarations.  This script moves
each namelist statement to just before the module-level `contains`
statement, which is semantically identical in Fortran.
"""
import re
import sys

FILES = [
    "fluid_model_hotspot_schnittman.f90",
    "fluid_model_harmpi.f90",
    "fluid_model_mb09.f90",
    "fluid_model_harm3d.f90",
    "fluid_model_harm.f90",
    "fluid_model_iharm.f90",
    "fluid_model_numdisk.f90",
    "fluid_model_hotspot.f90",
    "fluid_model_thindisk.f90",
    "fluid_model_phatdisk.f90",
    "fluid_model_koral.f90",
    "fluid_model_koral3d.f90",
    "fluid_model_ffjet.f90",
    "fluid_model_thickdisk.f90",
    "grtrans_program.f90",
]

NL = re.compile(r"^\s*namelist\s+/", re.IGNORECASE)
CONTAINS = re.compile(r"^\s*contains\s*$", re.IGNORECASE)

for fname in FILES:
    with open(fname) as f:
        lines = f.readlines()
    nlines = [l for l in lines if NL.match(l)]
    assert len(nlines) == 1, (fname, len(nlines))
    lines = [l for l in lines if not NL.match(l)]
    # insertion point: first module-level "contains"
    idx = None
    for i, l in enumerate(lines):
        if CONTAINS.match(l):
            idx = i
            break
    if idx is None:
        # program without contains: insert after the last declaration
        # (a line starting with a type specifier) before the first
        # executable statement
        lastdecl = 0
        for i, l in enumerate(lines):
            s = l.strip().lower()
            if s.startswith(("character", "integer", "real", "logical")):
                lastdecl = i
        idx = lastdecl + 1
    lines.insert(idx, nlines[0])
    with open(fname, "w") as f:
        f.writelines(lines)
    print(f"patched {fname} (inserted at line {idx+1})")
