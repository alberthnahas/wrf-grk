#!/usr/bin/env python3
"""Create start-date-aligned copies of monthly/8-day auxinput files.

WRF-Chem with `io_style_emissions = 2` does strict time matching: it opens
auxinput files by name and then verifies the internal `Times` variable
matches the model clock. Symlinking a non-matching name to a real file
fails because the internal `Times` is wrong (see Bug 12, Bug 15 in
docs/simulation.md).

This script makes real copies (not symlinks) of the nearest source file
and rewrites the internal `Times` variable to the simulation start. WRF
projection globals are also patched in.

Usage: edit START_DATE below to match your simulation start, then run.
The simulation start MUST be a real MODIS 8-day composite date for the
VPRM read schedule to land on real files (Bug 15) — i.e. one of:
    2015-07-20, 07-28, 08-05, 08-13, 08-21, 08-29, 09-06, 09-14, 09-22,
    09-30, 10-08, 10-16, 10-24, 11-01, 11-09, 11-17, 11-25.

If your simulation start is itself a MODIS composite, the VPRM file
already exists with the correct internal Times — no copy needed.
The wrfoce file is monthly and almost never matches the start date;
this script always creates a start-date-aligned copy of the nearest
month's file.

Idempotent: overwrites existing target if needed.
"""
import os
import sys
import shutil
from datetime import datetime
from glob import glob
import netCDF4 as nc

WRF_GRK = "/home/igrk/WRF-GRK"
SIM_DIR = os.path.join(WRF_GRK, "simulations/IDN_BB_2015")
INPUT = os.path.join(SIM_DIR, "input")
RUN = os.path.join(SIM_DIR, "run")
WRFINPUT = os.path.join(RUN, "wrfinput_d01")

START_DATE = datetime(2015, 7, 28, 0)
END_DATE   = datetime(2015, 12, 1, 0)
# auxinput6_interval_m from namelist (oce read interval)
OCE_INTERVAL_DAYS = 30

START_STR = START_DATE.strftime("%Y-%m-%d_%H:%M:%S")

# WRF projection globals (for patching the new files)
if not os.path.exists(WRFINPUT):
    print(f"ERROR: {WRFINPUT} not found. Run real.exe first.", file=sys.stderr)
    sys.exit(1)
with nc.Dataset(WRFINPUT) as s:
    src_attrs = {a: s.getncattr(a) for a in s.ncattrs() if a != "TITLE"}


def nearest_source_file(pattern_glob, start):
    """Pick the source file with date <= start, else the earliest available."""
    files = sorted(glob(pattern_glob))
    if not files:
        return None
    candidates = []
    for f in files:
        base = os.path.basename(f)
        # extract YYYY-MM-DD between the prefix and "_00:00:00"
        try:
            date_str = base.split("_d01_")[1].split("_00:00:00")[0]
            d = datetime.strptime(date_str, "%Y-%m-%d")
        except Exception:
            continue
        candidates.append((d, f))
    candidates.sort()
    le = [c for c in candidates if c[0] <= start]
    return (le[-1][1] if le else candidates[0][1])


def make_aligned_copy(src, dst_name, dst_dir, target_time_str=None):
    """Copy src → dst_dir/dst_name and rewrite the internal Times to target_time_str.

    target_time_str defaults to START_STR for backward compatibility, but callers
    should pass the actual per-file date — otherwise files at later read times
    will all carry Times = START, making WRF reject them with `wrf_get_next_time
    Status = -4` (see Bug 19/21).
    """
    dst = os.path.join(dst_dir, dst_name)
    if os.path.islink(dst) or os.path.exists(dst):
        os.remove(dst)
    shutil.copy2(src, dst)
    target = target_time_str if target_time_str is not None else START_STR
    with nc.Dataset(dst, "a") as d:
        # Rewrite Times[0] to the requested per-file target
        t = d.variables["Times"]
        for i, c in enumerate(target):
            t[0, i] = c
        # Patch projection globals
        own_title = d.getncattr("TITLE") if "TITLE" in d.ncattrs() else None
        for a, v in src_attrs.items():
            d.setncattr(a, v)
        if own_title is not None:
            d.setncattr("TITLE", own_title)
    return dst


def main():
    print(f"Simulation start: {START_STR}")
    print()

    # wrfoce — create copies at every auxinput6 read date (start + N × interval).
    # WRF computes the filename from current sim time, so a file must exist at
    # every read-time stamp, not just at start. With monthly source data and a
    # 30-day read interval, mid-month read-times like 2015-08-27 won't match any
    # source file unless we create one (see Bug 19).
    from datetime import timedelta
    interval = timedelta(days=OCE_INTERVAL_DAYS)
    t = START_DATE
    while t <= END_DATE:
        ts = t.strftime("%Y-%m-%d_%H:%M:%S")
        src_oce = nearest_source_file(os.path.join(INPUT, "wrfoce_d01_*"), t)
        if src_oce:
            dst = make_aligned_copy(src_oce, f"wrfoce_d01_{ts}", RUN, target_time_str=ts)
            print(f"  wrfoce[{ts}]: copied {os.path.basename(src_oce)} → {os.path.basename(dst)} (Times rewritten to {ts})")
        else:
            print(f"  WARN[{ts}]: no wrfoce source file found", file=sys.stderr)
        t += interval

    # vprm — only create copy if start is NOT a real MODIS composite
    vprm_src_exact = os.path.join(INPUT, f"vprm_input_d01_{START_STR}")
    if os.path.exists(vprm_src_exact):
        # Real composite — just ensure run-dir symlink
        link = os.path.join(RUN, f"vprm_input_d01_{START_STR}")
        if os.path.islink(link) or os.path.exists(link):
            os.remove(link)
        os.symlink(vprm_src_exact, link)
        print(f"  vprm: real MODIS composite exists ({os.path.basename(vprm_src_exact)}); symlinked into run/")
    else:
        src_vprm = nearest_source_file(os.path.join(INPUT, "vprm_input_d01_*"), START_DATE)
        if src_vprm:
            dst = make_aligned_copy(src_vprm, f"vprm_input_d01_{START_STR}", RUN, target_time_str=START_STR)
            print(f"  vprm: copied {os.path.basename(src_vprm)} → {os.path.basename(dst)} (Times rewritten)")
            print(f"  WARN: start {START_DATE.date()} is not a MODIS 8-day composite — restart segments may fail (see Bug 15)", file=sys.stderr)


if __name__ == "__main__":
    main()
