#!/usr/bin/env python3
"""Extract compact WRF-Chem essentials from segment wrfout files.

Default target:
  simulations/IDN_BB_2015/run/wrfessentials_seg1_idnbb2015.nc

The extraction keeps:
  - time/grid/eta coordinates needed for column work
  - essential meteorology and surface fields
  - all active chem_opt=17 GHG tracer fields
  - VPRM/ocean diagnostic flux fields if present

It intentionally refuses to run unless sufficient free space is available.
The threshold is conservative because the output is a single compressed NetCDF
created alongside large wrfout files.
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np


DEFAULT_WRFDIR = Path("/home/igrk/WRF-GRK/simulations/IDN_BB_2015/run")
DEFAULT_OUT = DEFAULT_WRFDIR / "wrfessentials_seg1_idnbb2015.nc"

COORD_VARS = [
    "Times", "XTIME", "XLAT", "XLONG", "ZNU", "ZNW",
]

ESSENTIAL_MET = [
    # 3D dynamics / thermodynamics
    "U", "V", "W", "T", "P", "PB", "PH", "PHB", "QVAPOR",
    "QCLOUD", "QRAIN", "QICE", "QSNOW", "CLDFRA",
    # 2D surface / column context
    "PSFC", "T2", "Q2", "U10", "V10", "TSK", "SST", "XLAND",
    "HFX", "LH", "PBLH", "UST", "SWDOWN", "GLW", "OLR",
    "RAINC", "RAINNC", "SNOW", "SNOWH",
    "VEGFRA", "LAI", "ALBEDO", "EMISS", "TMN", "COSZEN",
    "TSLB", "SMOIS", "SH2O", "GRDFLX",
]

GHG_TRACERS = [
    "CO2_ANT", "CO2_BIO", "CO2_OCE", "CO2_BCK", "CO2_BBU", "CO2_TST",
    "CH4_ANT", "CH4_BIO", "CH4_BCK", "CH4_BBU", "CH4_TST",
    "CO_ANT", "CO_BCK", "CO_BBU", "CO_TST",
]

GHG_FLUX_DIAGNOSTICS = [
    "EBIO_GEE", "EBIO_RES", "EBIO_CO2OCE", "EBIO_CH4WET", "EBIO_CH4SOIL",
    "RAD_VPRM", "LAMBDA_VPRM", "ALPHA_VPRM", "RESP_VPRM",
    # Some WRF builds preserve lowercase names in wrfout.
    "ebio_gee", "ebio_res", "ebio_co2oce", "ebio_ch4wet", "ebio_ch4soil",
    "rad_vprm", "lambda_vprm", "alpha_vprm", "resp_vprm",
]


def decode_times(ds: nc.Dataset) -> list[str]:
    return ["".join(row.astype(str)) for row in ds.variables["Times"][:]]


def complete_or_end_segment_file(path: str) -> bool:
    """Keep full 24-frame files and the 1-frame segment-end wrfout."""
    with nc.Dataset(path) as ds:
        times = decode_times(ds)
    if len(times) == 24 and times[0].endswith("00:00:00") and times[-1].endswith("23:00:00"):
        return True
    if len(times) == 1 and times[0] == "2015-08-29_00:00:00":
        return True
    return False


def selected_files(wrfdir: Path, pattern: str) -> list[str]:
    files = sorted(glob.glob(str(wrfdir / pattern)))
    return [f for f in files if complete_or_end_segment_file(f)]


def copy_attrs(src, dst) -> None:
    for attr in src.ncattrs():
        try:
            dst.setncattr(attr, src.getncattr(attr))
        except Exception:
            pass


def create_dimensions(out: nc.Dataset, first: nc.Dataset, total_time: int) -> None:
    for name, dim in first.dimensions.items():
        size = total_time if name == "Time" else len(dim)
        out.createDimension(name, size)


def compression_chunks(src_var):
    """Use time-slab chunks for practical hourly access and bounded memory."""
    if not src_var.dimensions or src_var.dimensions[0] != "Time":
        return None
    chunks = []
    for i, dim in enumerate(src_var.dimensions):
        if i == 0:
            chunks.append(1)
        else:
            chunks.append(len(src_var.get_dims()[i]))
    return tuple(chunks)


def create_variable(out: nc.Dataset, src_var, name: str, complevel: int, least_significant_digit: int | None):
    kwargs = {}
    if src_var.dtype.kind in ("f", "i", "u") and name != "Times":
        kwargs.update(zlib=True, complevel=complevel, shuffle=True)
        chunks = compression_chunks(src_var)
        if chunks is not None:
            kwargs["chunksizes"] = chunks
        if src_var.dtype.kind == "f" and least_significant_digit is not None:
            kwargs["least_significant_digit"] = least_significant_digit
    fill = getattr(src_var, "_FillValue", None)
    if fill is not None:
        kwargs["fill_value"] = fill
    dst = out.createVariable(name, src_var.dtype, src_var.dimensions, **kwargs)
    copy_attrs(src_var, dst)
    return dst


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wrfdir", type=Path, default=DEFAULT_WRFDIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--pattern", default="wrfout_d01_2015-*_00:00:00")
    parser.add_argument("--min-free-gb", type=float, default=320.0,
                        help="Refuse to run if output filesystem has less free space.")
    parser.add_argument("--complevel", type=int, default=5,
                        help="NetCDF4/zlib compression level, 1-9. Higher is smaller but slower.")
    parser.add_argument("--least-significant-digit", type=int, default=4,
                        help="Decimal digits retained for float quantization. Use -1 to disable.")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.complevel <= 9:
        raise ValueError("--complevel must be in [1, 9]")
    lsd = None if args.least_significant_digit < 0 else args.least_significant_digit

    files = selected_files(args.wrfdir, args.pattern)
    if not files:
        print(f"ERROR: no complete wrfout files found in {args.wrfdir}", file=sys.stderr)
        return 1

    free_gb = shutil.disk_usage(args.output.parent).free / 1024**3
    if free_gb < args.min_free_gb:
        print(
            f"ERROR: free space is {free_gb:.1f} GiB, below --min-free-gb {args.min_free_gb:.1f} GiB. "
            "Transfer/delete segment wrfout first, or rerun with a lower threshold if you accept the risk.",
            file=sys.stderr,
        )
        return 2

    total_time = 0
    time_ranges = []
    for f in files:
        with nc.Dataset(f) as ds:
            times = decode_times(ds)
        total_time += len(times)
        time_ranges.append((os.path.basename(f), len(times), times[0], times[-1]))

    with nc.Dataset(files[0]) as first:
        keep = []
        for name in COORD_VARS + ESSENTIAL_MET + GHG_TRACERS + GHG_FLUX_DIAGNOSTICS:
            if name in first.variables and name not in keep:
                keep.append(name)

        missing = [name for name in GHG_TRACERS if name not in first.variables]
        print(f"Input files: {len(files)}")
        print(f"Total frames: {total_time}")
        print(f"Output: {args.output}")
        print(f"Variables kept: {len(keep)}")
        print(f"Compression: zlib=True, shuffle=True, complevel={args.complevel}, least_significant_digit={lsd}")
        if missing:
            print(f"WARNING: missing expected GHG tracers: {', '.join(missing)}")
        if args.dry_run:
            for row in time_ranges[:3] + ([("...", "...", "...", "...")] if len(time_ranges) > 6 else []) + time_ranges[-3:]:
                print(row)
            print("Kept variables:", ", ".join(keep))
            return 0

        tmp = args.output.with_suffix(args.output.suffix + ".tmp")
        if tmp.exists():
            tmp.unlink()
        if args.output.exists():
            raise FileExistsError(f"Refusing to overwrite existing output: {args.output}")

        with nc.Dataset(tmp, "w", format="NETCDF4") as out:
            create_dimensions(out, first, total_time)
            copy_attrs(first, out)
            out.setncattr("TITLE", "WRF-Chem IDN_BB_2015 seg1 essentials for met/GHG/column analysis")
            out.setncattr("source_files", "\n".join(os.path.basename(f) for f in files))
            out.setncattr("extracted_variables", ", ".join(keep))
            out.setncattr(
                "column_requirements_note",
                "Dry-air total-column GHG can be computed from GHG tracers, QVAPOR, PSFC, ZNW, XLAT/XLONG, and P_TOP global attr.",
            )
            out_vars = {
                name: create_variable(out, first.variables[name], name, args.complevel, lsd)
                for name in keep
            }

            offset = 0
            for idx, f in enumerate(files, 1):
                with nc.Dataset(f) as src:
                    ntime = len(src.dimensions["Time"])
                    print(f"[{idx}/{len(files)}] {os.path.basename(f)} frames={ntime}", flush=True)
                    sl = slice(offset, offset + ntime)
                    for name in keep:
                        src_var = src.variables[name]
                        dst_var = out_vars[name]
                        if src_var.dimensions and src_var.dimensions[0] == "Time":
                            dst_var[sl, ...] = src_var[:]
                        elif offset == 0:
                            dst_var[...] = src_var[:]
                    offset += ntime

        tmp.rename(args.output)
        print(f"Done: {args.output}")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
