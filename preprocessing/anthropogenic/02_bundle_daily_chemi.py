#!/usr/bin/env python3
"""Bundle hourly single-frame wrfchemi files into daily 24-frame files.

Required because WRF v4.1.5 with `io_style_emissions = 2` and
`frames_per_auxinput5 = 1` aborts at startup with `wrf_get_next_time
Status = -4` on hourly single-frame anthro emission files (see Bug 11
in docs/simulation.md).

Reads:  rawdata of `01_process_edgar.py` →
        simulations/IDN_BB_2015/input/wrfchemi_d01_<date>_<HH>:00:00
Writes: simulations/IDN_BB_2015/input_daily/wrfchemi_d01_<date>_00:00:00
        (one file per day, 24 hourly frames, NETCDF3_CLASSIC)
        and creates symlinks in simulations/IDN_BB_2015/run/.

Idempotent: skips days whose output already exists with 24 frames.
"""
import os
import sys
import glob
from collections import defaultdict
from datetime import datetime
import netCDF4 as nc

WRF_GRK = "/home/igrk/WRF-GRK"
INPUT = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/input")
OUT = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/input_daily")
RUN = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/run")

os.makedirs(OUT, exist_ok=True)

src_files = sorted(glob.glob(os.path.join(INPUT, "wrfchemi_d01_*")))
print(f"Found {len(src_files)} hourly source files")

groups = defaultdict(list)
for f in src_files:
    base = os.path.basename(f).replace("wrfchemi_d01_", "")
    t = datetime.strptime(base, "%Y-%m-%d_%H:%M:%S")
    groups[t.date()].append((t, f))

print(f"{len(groups)} days")
written = 0

for day, items in sorted(groups.items()):
    items.sort()
    if len(items) != 24:
        print(f"  WARN day {day}: {len(items)} frames (expected 24)", file=sys.stderr)
    out_name = f"wrfchemi_d01_{day.strftime('%Y-%m-%d')}_00:00:00"
    out_path = os.path.join(OUT, out_name)

    # Idempotent skip if output already has 24 frames
    if os.path.exists(out_path):
        try:
            with nc.Dataset(out_path) as d:
                if len(d.dimensions["Time"]) == 24:
                    continue
        except Exception:
            pass
        os.remove(out_path)

    template = items[0][1]
    with nc.Dataset(template, "r") as s, nc.Dataset(out_path, "w", format="NETCDF3_CLASSIC") as d:
        for a in s.ncattrs():
            d.setncattr(a, s.getncattr(a))
        for name, dim in s.dimensions.items():
            if name == "Time":
                d.createDimension("Time", len(items))
            else:
                d.createDimension(name, len(dim))
        for name, var in s.variables.items():
            new = d.createVariable(name, var.dtype, var.dimensions)
            for a in var.ncattrs():
                new.setncattr(a, var.getncattr(a))
        for i, (_, fpath) in enumerate(items):
            with nc.Dataset(fpath, "r") as fr:
                for name, var in fr.variables.items():
                    if "Time" in var.dimensions:
                        d.variables[name][i, ...] = var[0, ...]
                    elif i == 0:
                        d.variables[name][...] = var[...]

    link = os.path.join(RUN, out_name)
    if os.path.islink(link) or os.path.exists(link):
        os.remove(link)
    os.symlink(out_path, link)
    written += 1

print(f"\nDone. {written} new daily files written to {OUT}")
print(f"Symlinks created/updated in {RUN}")
print()
print("Next step: run preprocessing/utils/01_patch_emiss_globals.py to add")
print("WRF projection metadata so wrf.exe accepts the files.")
