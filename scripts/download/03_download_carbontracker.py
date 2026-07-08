#!/usr/bin/env python3
"""
Download CarbonTracker CO2 and CH4 mole fractions using aria2c.
CO2: CT2025 (daily, glb3x2, flat directory)
CH4: CT-CH4-2025 (daily, glb3x2, year/month subdirs)
Period: 2015-06-01 to 2015-12-01
Uses aria2c: 16 connections per file, 8 files in parallel.
"""
import os
import datetime
import subprocess
import tempfile

OUT_DIR = "/home/igrk/WRF-GRK/rawdata/carbontracker"
os.makedirs(OUT_DIR, exist_ok=True)

NOAA_BASE = "https://gml.noaa.gov/aftp/products/carbontracker"
CO2_BASE  = f"{NOAA_BASE}/co2/CT2025/molefractions/co2_total"
CO2_PREFIX = "CT2025"
CH4_BASE  = f"{NOAA_BASE}/ch4/CT-CH4-2025/molefractions"
CH4_PREFIX = "CTCH4_2024"

START = datetime.date(2015, 6, 1)
END   = datetime.date(2015, 12, 1)

# Build list of (url, filename) pairs — skip already-complete files
entries = []

day = START
while day <= END:
    datestr = day.strftime("%Y-%m-%d")
    fname = f"{CO2_PREFIX}.molefrac_glb3x2_{datestr}.nc"
    dest  = os.path.join(OUT_DIR, fname)
    if not (os.path.exists(dest) and os.path.getsize(dest) > 50e6):  # CO2 files ~93MB
        entries.append((f"{CO2_BASE}/{fname}", fname))
    day += datetime.timedelta(days=1)

day = START
while day <= END:
    datestr = day.strftime("%Y-%m-%d")
    fname = f"{CH4_PREFIX}.molefrac_glb3x2_{datestr}.nc"
    dest  = os.path.join(OUT_DIR, fname)
    url   = f"{CH4_BASE}/{day.year}/{day.strftime('%m')}/{fname}"
    if not (os.path.exists(dest) and os.path.getsize(dest) > 20e6):  # CH4 files ~35MB
        entries.append((url, fname))
    day += datetime.timedelta(days=1)

print(f"Files to download: {len(entries)}")
already = sum(1 for f in os.listdir(OUT_DIR) if f.endswith(".nc"))
print(f"Already present:   {already}")

if not entries:
    print("All files already downloaded.")
else:
    # Write aria2c input file: url\n  out=filename\n
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tf:
        for url, fname in entries:
            tf.write(f"{url}\n  out={fname}\n")
        input_file = tf.name

    cmd = [
        "aria2c",
        f"--input-file={input_file}",
        f"--dir={OUT_DIR}",
        "--max-concurrent-downloads=4",   # 4 files at once
        "--split=2",                       # 2 segments per file
        "--max-connection-per-server=2",   # 2 connections per server
        "--min-split-size=20M",
        "--retry-wait=10",
        "--max-tries=10",
        "--file-allocation=none",
        "--check-certificate=false",
        "--auto-file-renaming=false",      # never create .1.nc duplicates
        "--console-log-level=notice",
        "--summary-interval=30",
    ]

    print("Starting aria2c download...")
    result = subprocess.run(cmd)
    os.unlink(input_file)

    if result.returncode != 0:
        print(f"aria2c exited with code {result.returncode}")
    else:
        print("aria2c download complete.")

nc_files = [f for f in os.listdir(OUT_DIR) if f.endswith(".nc")]
print(f"\nTotal NetCDF files in {OUT_DIR}: {len(nc_files)}")
