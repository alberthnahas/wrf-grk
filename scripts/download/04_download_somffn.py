#!/usr/bin/env python3
"""
Download SOMFFN ocean CO2 flux climatology
Landschutzer et al. - ocean pCO2 and air-sea flux
"""
import os
import requests
from tqdm import tqdm

OUT_DIR = "/home/igrk/WRF-GRK/rawdata/ocean"
os.makedirs(OUT_DIR, exist_ok=True)

# SOMFFN v2023 from NOAA OCADS (accession 0160558)
URLS = [
    "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0160558/MPI_SOM-FFN_v2023/MPI_SOM-FFN_v2023_NCEI_OCADS.nc",
]

# Fallback: v2022
FALLBACK_URLS = [
    "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0160558/MPI_SOM-FFN_v2022/MPI_SOM-FFN_v2022_NCEI_OCADS.nc",
]

def download_file(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 1e5:
        print(f"[SKIP] {os.path.basename(dest)} already exists")
        return True
    print(f"Downloading {os.path.basename(dest)}...")
    try:
        r = requests.get(url, stream=True, timeout=300)
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(dest, "wb") as f, tqdm(total=total, unit="B", unit_scale=True) as bar:
            for chunk in r.iter_content(65536):
                f.write(chunk)
                bar.update(len(chunk))
        print(f"  -> {dest}")
        return True
    except Exception as e:
        print(f"  ERROR: {e}")
        if os.path.exists(dest):
            os.remove(dest)
        return False

success = False
for url in URLS:
    fname = os.path.basename(url)
    if download_file(url, os.path.join(OUT_DIR, fname)):
        success = True
        break

if not success:
    print("Primary URL failed, trying fallback...")
    for url in FALLBACK_URLS:
        fname = os.path.basename(url)
        if download_file(url, os.path.join(OUT_DIR, fname)):
            success = True
            break

if not success:
    print("ERROR: Could not download SOMFFN data.")
    print("Please download manually from:")
    print("  https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/")
else:
    print("\nSOMFFN download complete.")
    for f in os.listdir(OUT_DIR):
        print(f"  {f}")
