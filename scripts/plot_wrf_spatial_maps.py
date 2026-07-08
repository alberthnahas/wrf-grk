#!/usr/bin/env python3
"""Create spatial maps for WRF-Chem essential meteorology and GHG variables.

The layout follows the local BMKG hotspot-map style: wide Indonesia extent,
province/land outlines, gridlines, an inset information panel, and clear source
metadata. Colormaps are selected per variable family rather than reusing the
hotspot colors.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import cartopy.crs as ccrs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import netCDF4 as nc
import numpy as np
from shapely.geometry import shape


DEFAULT_WRFDIR = Path("/home/igrk/WRF-GRK/simulations/IDN_BB_2015/run")
DEFAULT_OUTDIR = Path("/home/igrk/WRF-GRK/simulations/IDN_BB_2015/plots/wrf_diagnostics/maps")
GEOJSON = Path("/home/igrk/WRF-GRK/geojson")

GHG = [
    "CO2_TOTAL", "CO2_BCK", "CO2_ANT", "CO2_BIO", "CO2_OCE", "CO2_BBU", "CO2_TST",
    "CH4_TOTAL", "CH4_BCK", "CH4_ANT", "CH4_BIO", "CH4_BBU", "CH4_TST",
    "CO_TOTAL", "CO_BCK", "CO_ANT", "CO_BBU", "CO_TST",
]

METEO = [
    "T2", "PSFC", "WIND10", "U10", "V10", "PBLH", "HFX", "LH", "SWDOWN", "GLW", "OLR",
    "RAINC", "RAINNC", "SST", "TSK", "QVAPOR", "CLDFRA",
]

TOTAL_COMPONENTS = {
    "CO2_TOTAL": ["CO2_ANT", "CO2_BIO", "CO2_OCE", "CO2_BCK", "CO2_BBU", "CO2_TST"],
    "CH4_TOTAL": ["CH4_ANT", "CH4_BIO", "CH4_BCK", "CH4_BBU", "CH4_TST"],
    "CO_TOTAL": ["CO_ANT", "CO_BCK", "CO_BBU", "CO_TST"],
}

LABELS = {
    "CO2_TOTAL": ("Total CO2", "ppm"),
    "CH4_TOTAL": ("Total CH4", "ppb"),
    "CO_TOTAL": ("Total CO", "ppb"),
    "WIND10": ("10 m wind speed", "m s-1"),
}


def decode_times(ds: nc.Dataset) -> list[str]:
    return ["".join(row.astype(str)) for row in ds.variables["Times"][:]]


def latest_complete_file(wrfdir: Path) -> str:
    files = sorted(glob.glob(str(wrfdir / "wrfout_d01_*")))
    complete = []
    for f in files:
        try:
            with nc.Dataset(f) as ds:
                t = decode_times(ds)
            if len(t) == 24 and t[-1].endswith("23:00:00"):
                complete.append(f)
        except Exception:
            pass
    if not complete:
        raise FileNotFoundError(f"No complete wrfout files in {wrfdir}")
    return complete[-1]


def find_time_file(wrfdir: Path, target: str) -> tuple[str, int, str]:
    if target == "latest-complete":
        f = latest_complete_file(wrfdir)
        with nc.Dataset(f) as ds:
            t = decode_times(ds)
        return f, len(t) - 1, t[-1]
    for f in sorted(glob.glob(str(wrfdir / "wrfout_d01_*"))):
        with nc.Dataset(f) as ds:
            times = decode_times(ds)
        if target in times:
            return f, times.index(target), target
    raise ValueError(f"Could not find time {target}")


def surface_or_2d(var, tidx):
    dims = var.dimensions
    arr = var[tidx]
    if "bottom_top" in dims:
        axis = dims.index("bottom_top") - 1
        arr = np.take(arr, 0, axis=axis)
    elif "bottom_top_stag" in dims:
        axis = dims.index("bottom_top_stag") - 1
        arr = np.take(arr, 0, axis=axis)
    return np.asarray(arr)


def field(ds, name, tidx):
    if name in TOTAL_COMPONENTS:
        parts = []
        for v in TOTAL_COMPONENTS[name]:
            if v not in ds.variables:
                return None
            parts.append(surface_or_2d(ds.variables[v], tidx).astype("f8"))
        return sum(parts)
    if name == "WIND10":
        return np.hypot(ds.variables["U10"][tidx], ds.variables["V10"][tidx])
    if name not in ds.variables:
        return None
    return surface_or_2d(ds.variables[name], tidx)


def label_units(ds, name):
    if name in LABELS:
        return LABELS[name]
    if name in ds.variables:
        desc = getattr(ds.variables[name], "description", name)
        units = getattr(ds.variables[name], "units", "")
        return desc if desc else name, units
    return name, ""


def style_for(name, arr):
    lname = name.lower()
    if "bbu" in lname or "ant" in lname or "bio" in lname or "oce" in lname:
        return "magma", None
    if lname.startswith("co2") or lname.startswith("ch4") or lname.startswith("co_") or lname in ("co_total",):
        return "viridis", None
    if name in ("T2", "TSK", "SST"):
        return "coolwarm", None
    if name in ("HFX", "LH", "U10", "V10"):
        vmax = np.nanpercentile(np.abs(arr), 98)
        return "RdBu_r", TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)
    if name in ("RAINC", "RAINNC", "QVAPOR", "PBLH", "WIND10"):
        return "YlGnBu", None
    if name in ("PSFC", "SWDOWN", "GLW", "OLR"):
        return "plasma", None
    return "viridis", None


def geojson_geoms(path: Path):
    if not path.exists():
        return []
    with path.open() as fh:
        data = json.load(fh)
    if data.get("type") == "FeatureCollection":
        return [shape(f["geometry"]) for f in data.get("features", []) if f.get("geometry")]
    if data.get("type") == "GeometryCollection":
        return [shape(g) for g in data.get("geometries", [])]
    if "type" in data:
        return [shape(data)]
    return []


def add_overlays(ax, extent):
    proj = ccrs.PlateCarree()
    world = geojson_geoms(GEOJSON / "world_without_idn.json")
    prov = geojson_geoms(GEOJSON / "indonesia_38prov.geojson")
    peat = geojson_geoms(GEOJSON / "Indonesia_peat_lands.json")
    if world:
        ax.add_geometries(world, proj, facecolor="#eeeeee", edgecolor="#111111", linewidth=0.45, zorder=2)
    if prov:
        ax.add_geometries(prov, proj, facecolor="#d9d9d9", edgecolor="#111111", linewidth=0.65, zorder=3)
    if peat:
        ax.add_geometries(peat, proj, facecolor="none", edgecolor="#7a3b1a", linewidth=0.45, alpha=0.8, zorder=5)
    ax.set_extent(extent, crs=proj)
    gl = ax.gridlines(draw_labels=True, xlocs=np.arange(95, 146, 5), ylocs=np.arange(-15, 16, 5),
                      linewidth=0.45, color="gray", alpha=0.65, linestyle="--")
    gl.top_labels = False
    gl.right_labels = False


def robust_limits(arr, name):
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return None, None
    if name in ("CO2_BCK", "CH4_BCK", "CO_BCK", "CO2_TOTAL", "CH4_TOTAL", "CO_TOTAL", "T2", "SST", "TSK", "PSFC"):
        return np.nanpercentile(finite, 2), np.nanpercentile(finite, 98)
    return np.nanpercentile(finite, 1), np.nanpercentile(finite, 99)


def plot_one(ds, name, arr, lats, lons, timestamp, src_name, outdir):
    title, units = label_units(ds, name)
    cmap, norm = style_for(name, arr)
    vmin, vmax = robust_limits(arr, name)
    fig = plt.figure(figsize=(14, 9))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    add_overlays(ax, [95, 143, -15, 15])
    mesh = ax.pcolormesh(lons, lats, arr, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm,
                         vmin=None if norm else vmin, vmax=None if norm else vmax, shading="auto", zorder=4)
    cb = fig.colorbar(mesh, ax=ax, orientation="vertical", shrink=0.78, pad=0.015)
    cb.set_label(f"{name} ({units})" if units else name)

    finite = arr[np.isfinite(arr)]
    info = (
        f"{name}\n"
        f"{timestamp}\n"
        f"min {np.nanmin(finite):.3g}\n"
        f"mean {np.nanmean(finite):.3g}\n"
        f"max {np.nanmax(finite):.3g}"
    )
    ax.text(0.98, 0.98, info, transform=ax.transAxes, va="top", ha="right",
            fontsize=10, family="monospace",
            bbox=dict(facecolor="#fde5b5", edgecolor="black", linewidth=1.0, alpha=0.92))
    ax.text(0.03, 0.06,
            f"WRF-GHG IDN_BB_2015\nSource: {src_name}\nGeoJSON: provinces + peatlands",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=11, fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75))
    ax.set_title(title, fontsize=15, fontweight="bold", pad=14)
    fig.tight_layout()
    safe_time = timestamp.replace(":", "").replace("_", "-")
    fig.savefig(outdir / f"map_{name}_{safe_time}.png", dpi=180)
    plt.close(fig)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--wrfdir", type=Path, default=DEFAULT_WRFDIR)
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument("--time", default="latest-complete", help="'latest-complete' or exact WRF Times string")
    p.add_argument("--group", choices=["all", "ghg", "meteo"], default="all")
    p.add_argument("--variables", nargs="*", help="Explicit variables; overrides --group.")
    p.add_argument("--list-vars", action="store_true")
    args = p.parse_args()

    variables = args.variables
    if not variables:
        variables = GHG + METEO if args.group == "all" else (GHG if args.group == "ghg" else METEO)

    f, tidx, timestamp = find_time_file(args.wrfdir, args.time)
    with nc.Dataset(f) as ds:
        available = []
        for v in variables:
            if v in TOTAL_COMPONENTS or v == "WIND10" or v in ds.variables:
                available.append(v)
        if args.list_vars:
            print("Source:", f)
            print("Time:", timestamp)
            print("\n".join(available))
            return 0

        args.outdir.mkdir(parents=True, exist_ok=True)
        lats = ds.variables["XLAT"][tidx]
        lons = ds.variables["XLONG"][tidx]
        for i, v in enumerate(available, 1):
            arr = field(ds, v, tidx)
            if arr is None:
                continue
            print(f"[{i}/{len(available)}] {v}", flush=True)
            plot_one(ds, v, arr, lats, lons, timestamp, os.path.basename(f), args.outdir)
    print(f"Wrote maps to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
