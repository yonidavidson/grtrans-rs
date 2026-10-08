#!/usr/bin/env python3
"""Verify SHA-256 checksums of all reference fixtures.

Exits non-zero on any mismatch or missing file, so CI cannot silently skip
reference tests.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "reference" / "fixtures"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_case(d: Path, failures: list[str]) -> None:
    mfile = d / "manifest.json"
    if not mfile.exists():
        failures.append(f"{d.name}: missing manifest.json")
        return
    m = json.loads(mfile.read_text())
    for fn, meta in m["files"].items():
        p = d / fn
        if not p.exists():
            failures.append(f"{d.name}/{fn}: missing")
            continue
        if sha256(p) != meta["sha256"]:
            failures.append(f"{d.name}/{fn}: sha256 mismatch")
        raw = meta.get("raw")
        if raw:
            rp = d / raw
            if not rp.exists():
                failures.append(f"{d.name}/{raw}: missing")
            elif sha256(rp) != meta["raw_sha256"]:
                failures.append(f"{d.name}/{raw}: sha256 mismatch")
    for extra in ("inputs.in", "files.in"):
        if not (d / extra).exists():
            failures.append(f"{d.name}/{extra}: missing")


def main() -> int:
    failures: list[str] = []
    cases = [d for d in sorted(FIXTURES.iterdir()) if d.is_dir()]
    if not cases:
        print(f"no fixture directories found under {FIXTURES}", file=sys.stderr)
        return 1
    for d in cases:
        verify_case(d, failures)
    if failures:
        print("fixture verification FAILED:")
        for f in failures:
            print("  -", f)
        return 1
    n = len(cases)
    print(f"fixture verification OK ({n} cases)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
