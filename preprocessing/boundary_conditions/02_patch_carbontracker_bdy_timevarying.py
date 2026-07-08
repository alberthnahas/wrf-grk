#!/usr/bin/env python3
"""Patch time-varying CarbonTracker CO2/CH4 lateral BCs into wrfbdy_d01.

WRF boundary files carry both side values and tendencies:
  <VAR>_BXS/BXE/BYS/BYE   boundary values at each wrfbdy time
  <VAR>_BTXS/BTXE/BTYS/BTYE tendencies to the next wrfbdy time

The older 01_process_carbontracker.py wrote a single constant-in-time CT field
and set all BT arrays to zero. This script uses each wrfbdy Times entry, patches
all four sides, and computes BT from adjacent boundary times.

CO is intentionally left untouched here because this workspace does not contain
a CarbonTracker CO mole-fraction dataset; CO_BCK remains the documented 80 ppb
fallback from 01_process_carbontracker.py.
"""
import os
from datetime import datetime, timedelta

import netCDF4 as nc
import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator

WRF_GRK = "/home/igrk/WRF-GRK"
CT_DIR = os.path.join(WRF_GRK, "rawdata/carbontracker")
RUN_DIR = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/run")
WRFINPUT = os.path.join(RUN_DIR, "wrfinput_d01")
WRFBDY = os.path.join(RUN_DIR, "wrfbdy_d01")


def find_ct_files(species, year, month):
    prefix = "CT2025" if species == "CO2" else "CTCH4_2024"
    pattern = f"{prefix}.molefrac_glb3x2_{year}-{month:02d}"
    return sorted(
        os.path.join(CT_DIR, f)
        for f in os.listdir(CT_DIR)
        if f.startswith(pattern)
    )


def pick_ct_file(species, target_dt):
    files = find_ct_files(species, target_dt.year, target_dt.month)
    if target_dt.day <= 3:
        prev = target_dt.replace(day=1) - timedelta(days=1)
        files = find_ct_files(species, prev.year, prev.month) + files
    if target_dt.day >= 28:
        nxt = (target_dt.replace(day=28) + timedelta(days=4)).replace(day=1)
        files = files + find_ct_files(species, nxt.year, nxt.month)
    if not files:
        raise FileNotFoundError(f"No CarbonTracker files for {species} near {target_dt}")

    dated = []
    for f in files:
        date_s = os.path.basename(f).rsplit(".", 1)[0][-10:]
        try:
            dated.append((abs((datetime.strptime(date_s, "%Y-%m-%d") - target_dt).total_seconds()), f))
        except ValueError:
            pass
    candidates = [f for _, f in sorted(dated)] if dated else files
    for f in candidates:
        try:
            with nc.Dataset(f):
                pass
            return f
        except Exception as exc:
            print(f"  WARN: skipping unreadable {os.path.basename(f)} ({exc})")
    raise OSError(f"No readable CarbonTracker files for {species} near {target_dt}")


def load_ct_field(species, target_dt):
    path = pick_ct_file(species, target_dt)
    ds = xr.open_dataset(path, cache=False)
    try:
        needle = "co2" if species == "CO2" else "ch4"
        names = [v for v in ds.data_vars if needle in v.lower() and "total" in v.lower()]
        if not names:
            names = [v for v in ds.data_vars if needle in v.lower()]
        if not names:
            raise KeyError(f"No {species} variable found in {path}")
        var = ds[names[0]]

        lat = ds["latitude"].values if "latitude" in ds else ds["lat"].values
        lon = ds["longitude"].values if "longitude" in ds else ds["lon"].values

        data = var.values
        if "time" in var.dims:
            tidx = 0
            if "time" in ds.coords and ds["time"].size > 1:
                times = ds["time"].values
                target64 = np.datetime64(target_dt)
                tidx = int(np.argmin(np.abs(times - target64)))
            data = data[tidx]
        data = np.asarray(data).squeeze()

        if "pressure" in ds:
            p_edge = ds["pressure"].values
            p_edge_t = p_edge[0] if p_edge.ndim == 4 else p_edge
            p_mid = 0.5 * (p_edge_t[:-1] + p_edge_t[1:])
            lev = p_mid.mean(axis=(1, 2)) / 100.0
        elif "at" in ds and "bt" in ds and "surf_pressure" in ds:
            at = ds["at"].values
            bt = ds["bt"].values
            psfc = ds["surf_pressure"].values[0]
            p_edge = at[:, None, None] + bt[:, None, None] * psfc[None, :, :]
            p_mid = 0.5 * (p_edge[:-1] + p_edge[1:])
            lev = p_mid.mean(axis=(1, 2)) / 100.0
        elif "level" in ds:
            lev = ds["level"].values.astype(float)
        else:
            lev = None

        if lat[0] > lat[-1]:
            lat = lat[::-1]
            if data.ndim == 3:
                data = data[:, ::-1, :]
            elif data.ndim == 2:
                data = data[::-1, :]

        return data, lat, lon, lev
    finally:
        ds.close()


def boundary_points(wrf_lat, wrf_lon, wrf_pres, bdy_width):
    pieces = {
        "BXS": (wrf_lat[:, :bdy_width], wrf_lon[:, :bdy_width], wrf_pres[:, :, :bdy_width]),
        "BXE": (wrf_lat[:, -bdy_width:], wrf_lon[:, -bdy_width:], wrf_pres[:, :, -bdy_width:]),
        "BYS": (wrf_lat[:bdy_width, :], wrf_lon[:bdy_width, :], wrf_pres[:, :bdy_width, :]),
        "BYE": (wrf_lat[-bdy_width:, :], wrf_lon[-bdy_width:, :], wrf_pres[:, -bdy_width:, :]),
    }
    return pieces


def regrid_to_boundary(data, lat, lon, lev, target_lat, target_lon, target_pres):
    nz = target_pres.shape[0]
    shp2 = target_lat.shape
    pts = np.column_stack([target_lat.ravel(), target_lon.ravel()])

    if data.ndim == 2:
        interp = RegularGridInterpolator((lat, lon), data.astype(np.float32), bounds_error=False, fill_value=None)
        surface = interp(pts).astype(np.float32)
        return np.repeat(surface[None, :], nz, axis=0).reshape((nz,) + shp2)

    on_boundary = np.empty((data.shape[0], pts.shape[0]), dtype=np.float32)
    for k in range(data.shape[0]):
        interp = RegularGridInterpolator((lat, lon), data[k].astype(np.float32), bounds_error=False, fill_value=None)
        on_boundary[k] = interp(pts).astype(np.float32)

    if lev is None:
        return np.repeat(on_boundary[:1], nz, axis=0).reshape((nz,) + shp2)

    lev_use = np.asarray(lev, dtype=np.float64)
    vals_use = on_boundary
    if float(lev_use[0]) < float(lev_use[-1]):
        lev_use = lev_use[::-1]
        vals_use = vals_use[::-1]

    out = np.empty((nz, pts.shape[0]), dtype=np.float32)
    target_p = target_pres.reshape(nz, -1)
    xp = lev_use[::-1]
    for col in range(pts.shape[0]):
        out[:, col] = np.interp(target_p[:, col], xp, vals_use[::-1, col])
    return out.reshape((nz,) + shp2)


def to_bdy_order(suffix, arr):
    # arr is (bottom_top, south_north, west_east) for the boundary strip.
    # WRF boundary variables use (bdy_width, bottom_top, horizontal).
    if suffix in ("BXS", "BXE"):
        return np.transpose(arr, (2, 0, 1))
    return np.transpose(arr, (1, 0, 2))


with nc.Dataset(WRFINPUT) as ds:
    wrf_lat = ds.variables["XLAT"][0, :, :]
    wrf_lon = ds.variables["XLONG"][0, :, :]
    wrf_pres = (ds.variables["P"][0, :, :, :] + ds.variables["PB"][0, :, :, :]) / 100.0

with nc.Dataset(WRFBDY, "r+") as ds:
    times = [bytes(row).decode() for row in ds.variables["Times"][:]]
    dts = [datetime.strptime(t, "%Y-%m-%d_%H:%M:%S") for t in times]
    bdy_width = len(ds.dimensions["bdy_width"])
    pieces = boundary_points(wrf_lat, wrf_lon, wrf_pres, bdy_width)

    for species in ("CO2", "CH4"):
        print(f"Patching time-varying {species}_BCK for {len(dts)} wrfbdy times")
        previous = None
        previous_idx = None
        for idx, dt in enumerate(dts):
            data, lat, lon, lev = load_ct_field(species, dt)
            if species == "CO2" and np.nanmax(data) < 0.01:
                data = data * 1e6
            if species == "CH4" and np.nanmax(data) < 0.01:
                data = data * 1e9

            current = {}
            for suffix, (tlat, tlon, tpres) in pieces.items():
                vals = to_bdy_order(suffix, regrid_to_boundary(data, lat, lon, lev, tlat, tlon, tpres))
                ds.variables[f"{species}_BCK_{suffix}"][idx] = vals
                current[suffix] = vals

            if previous is not None:
                dt_seconds = (dt - dts[previous_idx]).total_seconds()
                for suffix, bt_suffix in (("BXS", "BTXS"), ("BXE", "BTXE"), ("BYS", "BTYS"), ("BYE", "BTYE")):
                    ds.variables[f"{species}_BCK_{bt_suffix}"][previous_idx] = (
                        current[suffix] - previous[suffix]
                    ) / dt_seconds
            previous = current
            previous_idx = idx

            if idx % 50 == 0 or idx == len(dts) - 1:
                print(f"  {species}: patched {idx + 1}/{len(dts)} ({times[idx]})")

        # No future boundary time is available for the last tendency.
        for bt_suffix in ("BTXS", "BTXE", "BTYS", "BTYE"):
            ds.variables[f"{species}_BCK_{bt_suffix}"][-1] = 0.0

print("Done.")
