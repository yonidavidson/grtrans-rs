#!/usr/bin/env python3
"""Extract the SPHACC profile tables from fluid_model_sphacc.f90.

Handles Fortran literals with a leading decimal point (e.g. `.13733`),
which a naive number regex would mis-parse as 13733.
"""
from __future__ import annotations

import re
import struct
import sys
from pathlib import Path

SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "~/dev/private/grtrans-rs-work/upstream/fluid_model_sphacc.f90"
).expanduser()
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(
    "crates/grtrans-physics/src/sphacc_tables.rs"
)

NUM = r"(?:\d+\.?\d*|\.\d+)(?:[dD][+-]?\d+)?"

text = SRC.read_text()
flat = " ".join(text.split())
out = {}
for name in ("tvals", "rvals", "vvals"):
    m = re.search(rf"{name}\s*=\s*\(/(.*?)/\)", flat)
    nums = re.findall(NUM, m.group(1))
    out[name] = [float(n.lower().replace("d", "e")) for n in nums]

uvals = [
    struct.unpack("f", struct.pack("f", (i - 1) * 0.001 + 1.0 / 399.0))[0]
    for i in range(1, 499)
]

lines = [
    "//! Generated from fluid_model_sphacc.f90 (jadexter/grtrans @ c76cb11).",
    "//! Do not edit; regenerate with scripts/extract_sphacc_tables.py.",
    "",
]
for name in ("rvals", "tvals", "vvals"):
    vals = out[name]
    lines.append(f"pub static {name.upper()}: [f64; {len(vals)}] = [")
    for i in range(0, len(vals), 4):
        lines.append("    " + ", ".join(repr(v) for v in vals[i : i + 4]) + ",")
    lines.append("];")
lines.append("")
lines.append(f"pub static UVALS: [f64; {len(uvals)}] = [")
for i in range(0, len(uvals), 4):
    lines.append("    " + ", ".join(repr(v) for v in uvals[i : i + 4]) + ",")
lines.append("];")

OUT.write_text("\n".join(lines) + "\n")
print(f"wrote {OUT}")
