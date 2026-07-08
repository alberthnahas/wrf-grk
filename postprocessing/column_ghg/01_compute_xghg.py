#!/usr/bin/env python3
"""
Post-processing: dry-air column-averaged GHG mole fractions from WRF-Chem chem_opt=17.

Outputs:
  XCO2  [ppm]  – pressure-weighted dry-air column average (total + 6 components)
  XCH4  [ppb]  – pressure-weighted dry-air column average (total + 5 components)
  XCO   [ppb]  – pressure-weighted dry-air column average (total + 4 components)

Comparable to: GOSAT, OCO-2, TROPOMI (without applying averaging kernels)

Usage:
    python 01_compute_xghg.py [--start YYYY-MM-DD] [--end YYYY-MM-DD] [--plot]

Output:
    simulations/IDN_BB_2015/output/xghg/xghg_<YYYY-MM-DD_HH>.nc  (one file per timestep)
    OR one combined file: xghg_IDN_BB_2015.nc  (all timesteps, --combine flag)
"""

import argparse
import os
import sys
import glob
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr
import netCDF4 as nc

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]
WRFOUT_DIR = ROOT / "simulations/IDN_BB_2015/run"
OUTPUT_DIR = ROOT / "simulations/IDN_BB_2015/output/xghg"
PLOT_DIR   = ROOT / "simulations/IDN_BB_2015/output/plots/xghg"

# Tracer names in wrfout (all ppmv dry-air mole fractions)
CO2_TRACERS = ["co2_ant", "co2_bio", "co2_oce", "co2_bck", "co2_bbu", "co2_tst"]
CH4_TRACERS = ["ch4_ant", "ch4_bio", "ch4_bck", "ch4_bbu", "ch4_tst"]
CO_TRACERS  = ["co_ant",  "co_bck",  "co_bbu",  "co_tst"]

# "background" tracer names that represent the boundary-condition field (not a source)
CO2_BG = "co2_bck"
CH4_BG = "ch4_bck"
CO_BG  = "co_bck"

# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def pressure_layer_thickness(ds):
    """
    Compute dry-air pressure layer thickness [Pa] for each WRF model level.

    Uses the exact terrain-following coordinate:
        dp[k] = (znw[k] - znw[k+1]) * (PSFC - P_TOP)

    Then the dry-air correction:
        dp_dry[k] = dp[k] / (1 + qvapor[k])

    Returns
    -------
    dp_dry : np.ndarray, shape (Time, nz, ny, nx)
    """
    znw  = ds["ZNW"].isel(Time=0).values          # (nz+1,) staggered eta
    psfc = ds["PSFC"].values                       # (Time, ny, nx)
    p_top = float(ds["P_TOP"].isel(Time=0).values) # scalar [Pa]
    qvap = ds["QVAPOR"].values                     # (Time, nz, ny, nx) [kg/kg]

    # Layer pressure thickness (total air)
    deta = znw[:-1] - znw[1:]                      # (nz,)
    deta_4d = deta[np.newaxis, :, np.newaxis, np.newaxis]
    psfc_4d = psfc[:, np.newaxis, :, :]
    dp = deta_4d * (psfc_4d - p_top)              # (Time, nz, ny, nx) [Pa]

    # Dry-air pressure thickness
    dp_dry = dp / (1.0 + qvap)                    # [Pa]
    return dp_dry


def column_average(tracer_3d, dp_dry):
    """
    Pressure-weighted dry-air column average.

    Parameters
    ----------
    tracer_3d : np.ndarray, shape (Time, nz, ny, nx)
    dp_dry    : np.ndarray, shape (Time, nz, ny, nx)

    Returns
    -------
    xgas : np.ndarray, shape (Time, ny, nx)
    """
    return (tracer_3d * dp_dry).sum(axis=1) / dp_dry.sum(axis=1)


def load_tracer(ds, name):
    """Load tracer, return zeros if not present (e.g. co2_tst when not used)."""
    if name in ds:
        return ds[name].values
    else:
        sample = ds[list(ds.data_vars)[0]].values
        print(f"  [WARN] tracer '{name}' not found, substituting zeros", flush=True)
        return np.zeros_like(sample)


def process_wrfout(wrfout_path, output_dir, plot=False):
    """
    Process one wrfout file (may contain multiple timesteps).
    Writes one output NC file per timestep.
    """
    print(f"Processing: {wrfout_path}", flush=True)
    ds = xr.open_dataset(wrfout_path, engine="netcdf4", decode_times=False)

    # --- coordinates ---
    xlat  = ds["XLAT"].isel(Time=0).values   # (ny, nx)
    xlong = ds["XLONG"].isel(Time=0).values  # (ny, nx)

    # Parse WRF times
    wrf_times_raw = ds["Times"].values        # (Time, 19) bytes
    wrf_times = []
    for t in wrf_times_raw:
        ts = bytes(t).decode("utf-8").replace("_", " ")
        wrf_times.append(datetime.strptime(ts, "%Y-%m-%d %H:%M:%S"))

    # --- pressure layer thicknesses ---
    dp_dry = pressure_layer_thickness(ds)     # (Time, nz, ny, nx)

    for ti, dt in enumerate(wrf_times):
        dp_i = dp_dry[[ti]]                   # (1, nz, ny, nx)

        # --- CO2 column [ppm] ---
        co2_parts = {}
        co2_total = np.zeros((1, ds.dims["bottom_top"], ds.dims["south_north"], ds.dims["west_east"]))
        for name in CO2_TRACERS:
            arr = load_tracer(ds, name)[[ti]]
            co2_parts[name] = column_average(arr, dp_i)[0]       # (ny, nx)
            co2_total += arr
        xco2_total = column_average(co2_total, dp_i)[0]          # (ny, nx) [ppm]

        # --- CH4 column [ppb] ---
        ch4_parts = {}
        ch4_total = np.zeros_like(co2_total)
        for name in CH4_TRACERS:
            arr = load_tracer(ds, name)[[ti]]
            ch4_parts[name] = column_average(arr, dp_i)[0] * 1e3  # ppmv → ppb
            ch4_total += arr
        xch4_total = column_average(ch4_total, dp_i)[0] * 1e3    # ppb

        # --- CO column [ppb] ---
        co_parts = {}
        co_total = np.zeros_like(co2_total)
        for name in CO_TRACERS:
            arr = load_tracer(ds, name)[[ti]]
            co_parts[name] = column_average(arr, dp_i)[0] * 1e3   # ppmv → ppb
            co_total += arr
        xco_total = column_average(co_total, dp_i)[0] * 1e3       # ppb

        # --- Write NetCDF ---
        out_name = output_dir / f"xghg_{dt.strftime('%Y-%m-%d_%H')}.nc"
        write_output(out_name, dt, xlat, xlong,
                     xco2_total, co2_parts,
                     xch4_total, ch4_parts,
                     xco_total,  co_parts)

        if plot:
            make_plots(out_name, dt)

    ds.close()


def write_output(out_path, dt, xlat, xlong,
                 xco2, co2_parts,
                 xch4, ch4_parts,
                 xco,  co_parts):
    """Write one timestep to NetCDF."""
    ny, nx = xlat.shape
    with nc.Dataset(out_path, "w", format="NETCDF4") as f:
        f.title = "WRF-Chem chem_opt=17 column-averaged GHG (XCO2, XCH4, XCO)"
        f.simulation = "IDN_BB_2015"
        f.timestamp = dt.strftime("%Y-%m-%d %H:%M:%S")
        f.created   = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
        f.note      = ("Dry-air pressure-weighted column average. "
                       "XCO2 [ppm], XCH4 [ppb], XCO [ppb]. "
                       "Comparable to GOSAT/OCO-2/TROPOMI (no averaging kernels applied).")

        f.createDimension("south_north", ny)
        f.createDimension("west_east", nx)

        def add_coord(name, data, units, desc):
            v = f.createVariable(name, "f4", ("south_north", "west_east"), zlib=True)
            v.units       = units
            v.description = desc
            v[:] = data

        add_coord("XLAT",  xlat,  "degrees_north", "Latitude")
        add_coord("XLONG", xlong, "degrees_east",  "Longitude")

        def add_2d(name, data, units, desc, longname=""):
            v = f.createVariable(name, "f4", ("south_north", "west_east"),
                                 zlib=True, complevel=4, fill_value=-9999.0)
            v.units       = units
            v.description = desc
            if longname:
                v.long_name = longname
            v[:] = data

        # XCO2
        add_2d("XCO2",     xco2, "ppm", "Dry-air column-averaged CO2 (total all tracers)",
               "XCO2 = co2_ant+co2_bio+co2_oce+co2_bck+co2_bbu+co2_tst column average")
        for name, arr in co2_parts.items():
            add_2d(f"XCO2_{name.upper()}", arr, "ppm",
                   f"Column average of {name} component")

        # XCH4
        add_2d("XCH4",     xch4, "ppb", "Dry-air column-averaged CH4 (total all tracers)",
               "XCH4 = ch4_ant+ch4_bio+ch4_bck+ch4_bbu+ch4_tst column average")
        for name, arr in ch4_parts.items():
            add_2d(f"XCH4_{name.upper()}", arr, "ppb",
                   f"Column average of {name} component")

        # XCO
        add_2d("XCO",      xco,  "ppb", "Dry-air column-averaged CO (total all tracers)",
               "XCO = co_ant+co_bck+co_bbu+co_tst column average")
        for name, arr in co_parts.items():
            add_2d(f"XCO_{name.upper()}", arr, "ppb",
                   f"Column average of {name} component")

    print(f"  Wrote: {out_path}", flush=True)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_plots(nc_path, dt):
    """Generate PNG maps for XCO2, XCH4, XCO."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
    except ImportError:
        print("  [WARN] matplotlib/cartopy not available — skipping plots", flush=True)
        return

    ds = xr.open_dataset(nc_path)
    lat  = ds["XLAT"].values
    lon  = ds["XLONG"].values
    proj = ccrs.PlateCarree()

    plot_vars = [
        ("XCO2",  "XCO₂",  "ppm",  "RdYlBu_r", None,  None),
        ("XCH4",  "XCH₄",  "ppb",  "RdYlBu_r", None,  None),
        ("XCO",   "XCO",   "ppb",  "hot_r",    None,  None),
    ]

    plot_out = PLOT_DIR / dt.strftime("%Y%m")
    plot_out.mkdir(parents=True, exist_ok=True)

    for varname, title, units, cmap, vmin, vmax in plot_vars:
        if varname not in ds:
            continue
        data = ds[varname].values
        fig, ax = plt.subplots(figsize=(12, 6),
                               subplot_kw={"projection": proj})
        ax.set_extent([89, 146, -16, 16], crs=proj)
        ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=0)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.4, linestyle=":")

        im = ax.pcolormesh(lon, lat, data,
                           transform=proj, cmap=cmap,
                           vmin=vmin, vmax=vmax, shading="auto")
        cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
        cbar.set_label(units)

        gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="gray", alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlocator = mticker.FixedLocator(range(90, 150, 10))
        gl.ylocator = mticker.FixedLocator(range(-15, 20, 5))

        ax.set_title(f"{title}  –  {dt.strftime('%Y-%m-%d %H:%M UTC')}", fontsize=12)
        outfile = plot_out / f"{varname.lower()}_{dt.strftime('%Y%m%d_%H')}.png"
        fig.savefig(outfile, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Plot: {outfile}", flush=True)

    # --- Fire contribution plot: anomaly from background ---
    try:
        bbu_co2  = ds["XCO2_CO2_BBU"].values  if "XCO2_CO2_BBU"  in ds else None
        bbu_ch4  = ds["XCH4_CH4_BBU"].values  if "XCH4_CH4_BBU"  in ds else None
        bbu_co   = ds["XCO_CO_BBU"].values    if "XCO_CO_BBU"    in ds else None

        if bbu_co is not None:
            fig, axes = plt.subplots(1, 3, figsize=(18, 5),
                                     subplot_kw={"projection": proj})
            items = [
                (axes[0], bbu_co2,  "Fire ΔXCo₂", "ppm",  "YlOrRd"),
                (axes[1], bbu_ch4,  "Fire ΔXCH₄", "ppb",  "YlOrRd"),
                (axes[2], bbu_co,   "Fire ΔXCO",   "ppb",  "YlOrRd"),
            ]
            for ax, data, ttl, units, cmap in items:
                if data is None:
                    continue
                ax.set_extent([89, 146, -16, 16], crs=proj)
                ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=0)
                ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
                im = ax.pcolormesh(lon, lat, data, transform=proj,
                                   cmap=cmap, vmin=0, shading="auto")
                cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
                cbar.set_label(units)
                ax.set_title(f"{ttl}  {dt.strftime('%Y-%m-%d %H:%M')}", fontsize=10)
            fig.suptitle("Biomass burning GHG column enhancement", fontsize=12)
            outfile = plot_out / f"fire_xghg_{dt.strftime('%Y%m%d_%H')}.png"
            fig.savefig(outfile, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"  Plot: {outfile}", flush=True)
    except Exception as e:
        print(f"  [WARN] fire plot failed: {e}", flush=True)

    ds.close()


# ---------------------------------------------------------------------------
# Combine output files into a single NetCDF with time dimension
# ---------------------------------------------------------------------------

def combine_output(output_dir, combined_path):
    """Merge all hourly xghg_*.nc files into one file with a time axis."""
    files = sorted(glob.glob(str(output_dir / "xghg_????-??-??_??.nc")))
    if not files:
        print("No files to combine.", flush=True)
        return

    print(f"Combining {len(files)} files → {combined_path}", flush=True)
    datasets = []
    times = []
    for f in files:
        ts_str = Path(f).stem.replace("xghg_", "").replace("_", " ") + ":00:00"
        times.append(datetime.strptime(ts_str, "%Y-%m-%d %H:%M:%S"))
        ds = xr.open_dataset(f)
        # Add expand_dims for time
        ds = ds.expand_dims("time")
        datasets.append(ds)

    combined = xr.concat(datasets, dim="time")
    combined["time"] = [np.datetime64(t) for t in times]
    combined.to_netcdf(combined_path, format="NETCDF4",
                       encoding={v: {"zlib": True, "complevel": 4}
                                 for v in combined.data_vars
                                 if v not in ("XLAT", "XLONG")})
    print(f"  Done: {combined_path}", flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compute XCO2, XCH4, XCO from WRF-Chem output")
    parser.add_argument("--start",   default="2015-07-25", help="Start date YYYY-MM-DD (default: 2015-07-25)")
    parser.add_argument("--end",     default="2015-12-01", help="End date YYYY-MM-DD (default: 2015-12-01)")
    parser.add_argument("--wrfdir",  default=str(WRFOUT_DIR), help="Directory containing wrfout files")
    parser.add_argument("--outdir",  default=str(OUTPUT_DIR), help="Output directory for xghg NC files")
    parser.add_argument("--plot",    action="store_true", help="Generate PNG plots")
    parser.add_argument("--combine", action="store_true", help="Combine hourly files into one NC at end")
    args = parser.parse_args()

    wrfdir  = Path(args.wrfdir)
    outdir  = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.plot:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)

    # Find wrfout files
    all_files = sorted(glob.glob(str(wrfdir / "wrfout_d01_*")))
    if not all_files:
        print(f"ERROR: No wrfout files found in {wrfdir}", file=sys.stderr)
        sys.exit(1)

    # Filter by date range
    t_start = datetime.strptime(args.start, "%Y-%m-%d")
    t_end   = datetime.strptime(args.end,   "%Y-%m-%d")
    selected = []
    for f in all_files:
        # filename: wrfout_d01_YYYY-MM-DD_HH:MM:SS
        try:
            ts_str = Path(f).name[10:].replace("_", " ")
            ts = datetime.strptime(ts_str, "%Y-%m-%d %H:%M:%S")
            if t_start <= ts <= t_end:
                selected.append(f)
        except ValueError:
            selected.append(f)  # can't parse — include anyway

    print(f"Found {len(selected)} wrfout files in [{args.start}, {args.end}]", flush=True)

    # Process each file
    for f in selected:
        try:
            process_wrfout(Path(f), outdir, plot=args.plot)
        except Exception as e:
            print(f"  [ERROR] {f}: {e}", file=sys.stderr)
            import traceback; traceback.print_exc()

    # Optionally combine
    if args.combine:
        combined_path = outdir / "xghg_IDN_BB_2015.nc"
        combine_output(outdir, combined_path)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
