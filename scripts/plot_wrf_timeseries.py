#!/usr/bin/env python3
"""Create WRF-Chem essential-variable timeseries from wrfout files.

The script reads wrfout directly and writes:
  - one CSV with hourly min/mean/max/p95 statistics
  - one PNG per variable showing hourly mean and max

For 3D fields, the default statistic is computed at the lowest model level.
This keeps the product useful for monitoring GHG and meteorology without
accidentally doing expensive full-column diagnostics.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


DEFAULT_WRFDIR = Path("/home/igrk/WRF-GRK/simulations/IDN_BB_2015/run")
DEFAULT_OUTDIR = Path("/home/igrk/WRF-GRK/simulations/IDN_BB_2015/plots/wrf_diagnostics/timeseries")

GHG = [
    "CO2_ANT", "CO2_BIO", "CO2_OCE", "CO2_BCK", "CO2_BBU", "CO2_TST",
    "CH4_ANT", "CH4_BIO", "CH4_BCK", "CH4_BBU", "CH4_TST",
    "CO_ANT", "CO_BCK", "CO_BBU", "CO_TST",
]

METEO = [
    "T2", "PSFC", "U10", "V10", "PBLH", "HFX", "LH", "SWDOWN", "GLW", "OLR",
    "RAINC", "RAINNC", "SST", "TSK", "QVAPOR", "T", "P", "PB", "U", "V", "W",
    "CLDFRA", "QCLOUD", "QRAIN", "QICE", "QSNOW",
]

DERIVED = ["CO2_TOTAL", "CH4_TOTAL", "CO_TOTAL", "WIND10"]

TOTAL_COMPONENTS = {
    "CO2_TOTAL": ["CO2_ANT", "CO2_BIO", "CO2_OCE", "CO2_BCK", "CO2_BBU", "CO2_TST"],
    "CH4_TOTAL": ["CH4_ANT", "CH4_BIO", "CH4_BCK", "CH4_BBU", "CH4_TST"],
    "CO_TOTAL": ["CO_ANT", "CO_BCK", "CO_BBU", "CO_TST"],
}

LABELS = {
    "CO2_TOTAL": "Total surface CO2 (ppm)",
    "CH4_TOTAL": "Total surface CH4 (ppb)",
    "CO_TOTAL": "Total surface CO (ppb)",
    "WIND10": "10 m wind speed (m s-1)",
}


def decode_times(ds: nc.Dataset) -> list[str]:
    return ["".join(row.astype(str)) for row in ds.variables["Times"][:]]


def wrfout_files(wrfdir: Path, complete_only: bool) -> list[str]:
    files = sorted(glob.glob(str(wrfdir / "wrfout_d01_*")))
    if not complete_only:
        return files
    keep = []
    for f in files:
        try:
            with nc.Dataset(f) as ds:
                t = decode_times(ds)
            if len(t) == 24 and t[0].endswith("00:00:00") and t[-1].endswith("23:00:00"):
                keep.append(f)
            elif len(t) == 1 and t[0].endswith("00:00:00"):
                keep.append(f)
        except Exception:
            pass
    return keep


def surface_or_2d(var):
    arr = var[:]
    dims = var.dimensions
    if "bottom_top" in dims:
        axis = dims.index("bottom_top")
        arr = np.take(arr, 0, axis=axis)
    elif "bottom_top_stag" in dims:
        axis = dims.index("bottom_top_stag")
        arr = np.take(arr, 0, axis=axis)
    return np.asarray(arr)


def read_field(ds: nc.Dataset, name: str):
    if name in TOTAL_COMPONENTS:
        parts = []
        for v in TOTAL_COMPONENTS[name]:
            if v not in ds.variables:
                return None
            parts.append(surface_or_2d(ds.variables[v]).astype("f8"))
        return sum(parts)
    if name == "WIND10":
        if "U10" not in ds.variables or "V10" not in ds.variables:
            return None
        return np.hypot(ds.variables["U10"][:], ds.variables["V10"][:])
    if name not in ds.variables:
        return None
    return surface_or_2d(ds.variables[name])


def stats_for(arr):
    flat = arr.reshape((arr.shape[0], -1))
    return {
        "min": np.nanmin(flat, axis=1),
        "mean": np.nanmean(flat, axis=1),
        "p95": np.nanpercentile(flat, 95, axis=1),
        "max": np.nanmax(flat, axis=1),
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--wrfdir", type=Path, default=DEFAULT_WRFDIR)
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument("--group", choices=["all", "ghg", "meteo"], default="all")
    p.add_argument("--variables", nargs="*", help="Explicit variable list; overrides --group.")
    p.add_argument("--include-partial", action="store_true")
    args = p.parse_args()

    variables = args.variables
    if not variables:
        variables = (DERIVED + GHG + METEO) if args.group == "all" else (DERIVED + GHG if args.group == "ghg" else METEO)

    files = wrfout_files(args.wrfdir, complete_only=not args.include_partial)
    if not files:
        raise SystemExit(f"No wrfout files found in {args.wrfdir}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    csv_path = args.outdir / "wrf_hourly_essential_timeseries.csv"
    rows = []
    series = {v: {"time": [], "mean": [], "max": []} for v in variables}

    for i, f in enumerate(files, 1):
        print(f"[{i}/{len(files)}] {os.path.basename(f)}", flush=True)
        with nc.Dataset(f) as ds:
            times = decode_times(ds)
            for v in variables:
                arr = read_field(ds, v)
                if arr is None:
                    continue
                s = stats_for(arr)
                for ti, time in enumerate(times):
                    rows.append({
                        "time": time,
                        "file": os.path.basename(f),
                        "variable": v,
                        "min": float(s["min"][ti]),
                        "mean": float(s["mean"][ti]),
                        "p95": float(s["p95"][ti]),
                        "max": float(s["max"][ti]),
                    })
                    series[v]["time"].append(time)
                    series[v]["mean"].append(float(s["mean"][ti]))
                    series[v]["max"].append(float(s["max"][ti]))

    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["time", "file", "variable", "min", "mean", "p95", "max"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {csv_path}")

    for v, data in series.items():
        if not data["time"]:
            continue
        x = np.arange(len(data["time"]))
        fig, ax = plt.subplots(figsize=(12, 4.5))
        ax.plot(x, data["mean"], color="#1f77b4", lw=1.5, label="domain mean")
        ax.plot(x, data["max"], color="#d62728", lw=0.9, alpha=0.75, label="domain max")
        tick_step = max(1, len(x) // 10)
        ax.set_xticks(x[::tick_step])
        ax.set_xticklabels([data["time"][j][:10] for j in x[::tick_step]], rotation=35, ha="right")
        ax.set_title(LABELS.get(v, v), fontweight="bold")
        ax.set_ylabel(v)
        ax.grid(True, ls="--", lw=0.5, alpha=0.4)
        ax.legend(loc="best")
        fig.tight_layout()
        fig.savefig(args.outdir / f"timeseries_{v}.png", dpi=160)
        plt.close(fig)
    print(f"Wrote PNG plots to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
