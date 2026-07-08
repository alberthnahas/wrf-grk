#!/usr/bin/env python3
"""Patch WRF projection global attributes onto chem-emission auxinput files.

Required because the EDGAR/FINN/SOMFFN/VPRM preprocessing scripts only write
`:TITLE = "WRF-Chem EMISSIONS"` as a global attribute, but WRF-Chem 4.1.5's
input reader expects the full set of WRF projection/grid attributes
(`DX`, `DY`, `MAP_PROJ`, dimensions, etc.) — without them WRF reports
"Error trying to read metadata" and fails with `Status = -4` on the first
auxinput read (see Bug 11 in docs/simulation.md).

Source of attributes: simulations/IDN_BB_2015/run/wrfinput_d01 (after real.exe).
Targets (auxinput files):
    input_daily/wrfchemi_d01_*       (auxinput5 — anthro, daily 24-frame)
    input/wrffirechemi_d01_*         (auxinput7 — fire, daily 24-frame)
    input/wrfoce_d01_*               (auxinput6 — ocean, monthly)
    input/vprm_input_d01_*           (auxinput15 — VPRM, 8-day MODIS)
Each file's own `:TITLE` is preserved.

Idempotent: re-running just rewrites the same attributes.
"""
import os
import sys
import glob
import netCDF4 as nc

WRF_GRK = "/home/igrk/WRF-GRK"
SIM_DIR = os.path.join(WRF_GRK, "simulations/IDN_BB_2015")
SRC = os.path.join(SIM_DIR, "run", "wrfinput_d01")

if not os.path.exists(SRC):
    print(f"ERROR: {SRC} not found. Run real.exe first.", file=sys.stderr)
    sys.exit(1)

with nc.Dataset(SRC, "r") as s:
    src_attrs = {a: s.getncattr(a) for a in s.ncattrs()}

PATTERNS = [
    os.path.join(SIM_DIR, "input_daily", "wrfchemi_d01_*"),
    os.path.join(SIM_DIR, "input", "wrffirechemi_d01_*"),
    os.path.join(SIM_DIR, "input", "wrfoce_d01_*"),
    os.path.join(SIM_DIR, "input", "vprm_input_d01_*"),
]

files = []
for p in PATTERNS:
    files += sorted(glob.glob(p))
files = [f for f in files if not os.path.isdir(f)]

print(f"Patching {len(files)} files with {len(src_attrs)} globals from {os.path.basename(SRC)}")

patched = 0
for f in files:
    real = os.path.realpath(f)
    try:
        with nc.Dataset(real, "a") as d:
            own_title = d.getncattr("TITLE") if "TITLE" in d.ncattrs() else None
            for a, v in src_attrs.items():
                if a == "TITLE":
                    continue
                d.setncattr(a, v)
            if own_title is not None:
                d.setncattr("TITLE", own_title)
        patched += 1
    except Exception as e:
        print(f"  FAIL {f}: {e}", file=sys.stderr)

print(f"Done. {patched}/{len(files)} files patched.")
