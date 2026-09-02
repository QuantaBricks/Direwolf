#!/usr/bin/env python3
# Copyright (c) 2026 QuantaBricks
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Run every functional-accuracy test in this directory and diff
against references.csv. See README.md for what this set covers and
why a few entries carry a documented, unresolved residual instead of a
tight tolerance.

Usage (from anywhere - paths below are relative to this script's own
location, not the caller's cwd):
    python3 run_benchmark.py [--engine PATH_TO_ENGINE_BINARY]

Exit code is 0 only if every row is within its own tolerance.
"""
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENGINE_DEFAULT = HERE.parent.parent / "Direwolf"


def run_one(engine_path, inp_path, out_path):
    result = subprocess.run(
        [str(engine_path), str(inp_path), str(out_path)],
        capture_output=True, text=True, timeout=600,
    )
    if result.returncode != 0:
        return None, f"Direwolf exited {result.returncode}: {result.stderr.strip()[:200]}"
    text = out_path.read_text()
    m = re.search(r"Total Energy \(Hartree\)\s*=\s*(-?\d+\.\d+)", text)
    if not m:
        return None, "no Total Energy (Hartree) line in output"
    return float(m.group(1)), None


def main():
    engine_path = ENGINE_DEFAULT
    if len(sys.argv) > 1 and sys.argv[1] == "--engine":
        engine_path = Path(sys.argv[2])
    if not engine_path.exists():
        print(f"Direwolf binary not found at {engine_path} - build it first (make, from the repo root).")
        return 1

    with open(HERE / "references.csv") as f:
        rows = list(csv.DictReader(f))

    all_ok = True
    width = max(len(r["functional"]) for r in rows)
    print(f"{'Functional':<{width}}  {'Direwolf (Ha)':>18}  {'Reference (Ha)':>18}  {'diff':>10}  {'tol':>8}  status")
    with tempfile.TemporaryDirectory() as tmpdir:
        for row in rows:
            inp = HERE / row["input"]
            out = Path(tmpdir) / (row["input"] + ".out")
            energy, err = run_one(engine_path, inp, out)
            tol = float(row["tolerance_hartree"])
            if energy is None:
                print(f"{row['functional']:<{width}}  {'ERROR':>18}  {'':>18}  {'':>10}  {tol:>8.1e}  FAIL ({err})")
                all_ok = False
                continue
            ref = float(row["reference_hartree"])
            diff = energy - ref
            ok = (tol == 0.0) or (abs(diff) <= tol)
            status = "ok" if ok else "OUT OF TOLERANCE"
            if not ok:
                all_ok = False
            note = f"  [{row['notes']}]" if row["notes"] else ""
            print(f"{row['functional']:<{width}}  {energy:>18.9f}  {ref:>18.9f}  {diff:>10.2e}  {tol:>8.1e}  {status}{note}")

    print()
    print("All within tolerance." if all_ok else "Some rows are out of tolerance - see above (a KNOWN-residual note doesn't excuse a NEW regression beyond its documented tolerance).")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
