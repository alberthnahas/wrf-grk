#!/usr/bin/env python3
"""
Download EDGAR v8.0 (GHG) + HTAP v3 (CO) gridded emissions for 2015.
Files come as ZIP archives; NetCDF files are extracted to rawdata/edgar/.

CO2: EDGAR v8.0 FT2022 GHG — per-sector annual files (~115KB each zip, ~20 sectors)
CH4: EDGAR v8.0 FT2022 GHG — per-sector annual files (~115KB each zip, ~22 sectors)
     AWB (agricultural waste burning) excluded from CH4 — covered by FINN fire emissions.
CO:  EDGAR HTAP v3 — monthly sector gridmaps 0.1x0.1 (~365MB zip)

TOTALS files (CO2/CH4) are also kept as validation/fallback.
"""
import os
import zipfile
import requests
from tqdm import tqdm

OUT_DIR = "/home/igrk/WRF-GRK/rawdata/edgar"
os.makedirs(OUT_DIR, exist_ok=True)

GHG_BASE = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v80_FT2022_GHG"
HTAP_BASE = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/htap_v3/gridmaps_01x01/emissions"

# --- CO2 sectors (EDGAR v8.0, all except biogenic)
CO2_SECTORS = [
    "ENE", "IND", "CHE", "NMM", "IRO", "NFE", "NEU",
    "TRO", "RCO", "REF_TRF", "PRO_FFF", "PRU_SOL",
    "SWD_INC", "AGS",
    "TNR_Ship", "TNR_Aviation_LTO", "TNR_Aviation_CDS",
    "TNR_Aviation_CRS", "TNR_Other",
    # TNR_Aviation_SPS: not available in EDGAR v8.0 (404)
]

# --- CH4 sectors (EDGAR v8.0)
# AWB (agricultural waste burning) deliberately excluded — covered by FINN fire emissions
# TNR_Aviation_SPS: not available in EDGAR v8.0 (404)
CH4_SECTORS = [
    "ENE", "PRO_FFF", "PRO_COAL", "PRO_GAS", "PRO_OIL",
    "IND", "CHE", "IRO", "REF_TRF", "RCO", "TRO",
    "ENF", "MNM", "AGS",
    "SWD_LDF", "SWD_INC", "WWT",
    "TNR_Ship", "TNR_Aviation_LTO", "TNR_Aviation_CDS",
    "TNR_Aviation_CRS", "TNR_Other",
]

# Build file list
FILES = [
    # TOTALS (kept for validation/fallback)
    (f"{GHG_BASE}/CO2/TOTALS/emi_nc/v8.0_FT2022_GHG_CO2_2015_TOTALS_emi_nc.zip",
     "v8.0_FT2022_GHG_CO2_2015_TOTALS_emi_nc.zip", 10_000_000),
    (f"{GHG_BASE}/CH4/TOTALS/emi_nc/v8.0_FT2022_GHG_CH4_2015_TOTALS_emi_nc.zip",
     "v8.0_FT2022_GHG_CH4_2015_TOTALS_emi_nc.zip", 10_000_000),
    # HTAP v3 CO (monthly, sector-resolved)
    (f"{HTAP_BASE}/CO/edgar_HTAPv3_2015_CO.zip",
     "edgar_HTAPv3_2015_CO.zip", 200_000_000),
]
# Add per-sector CO2 files (~115KB each)
for sec in CO2_SECTORS:
    fname = f"v8.0_FT2022_GHG_CO2_2015_{sec}_emi_nc.zip"
    url = f"{GHG_BASE}/CO2/{sec}/emi_nc/{fname}"
    FILES.append((url, fname, 50_000))
# Add per-sector CH4 files (~115KB each)
for sec in CH4_SECTORS:
    fname = f"v8.0_FT2022_GHG_CH4_2015_{sec}_emi_nc.zip"
    url = f"{GHG_BASE}/CH4/{sec}/emi_nc/{fname}"
    FILES.append((url, fname, 50_000))


def download_file(url, dest, min_size=10_000_000):
    if os.path.isfile(dest) and os.path.getsize(dest) > min_size:
        print(f"[SKIP] {os.path.basename(dest)} already downloaded")
        return
    print(f"Downloading {os.path.basename(dest)}...")
    r = requests.get(url, stream=True, timeout=300)
    if r.status_code == 404:
        print(f"  [SKIP] 404 Not Found — sector not available in EDGAR v8.0")
        return
    r.raise_for_status()
    total = int(r.headers.get("content-length", 0))
    with open(dest, "wb") as f, tqdm(total=total, unit="B", unit_scale=True) as bar:
        for chunk in r.iter_content(65536):
            f.write(chunk)
            bar.update(len(chunk))
    print(f"  saved -> {dest}")


def extract_nc(zip_path, out_dir):
    print(f"Extracting {os.path.basename(zip_path)}...")
    with zipfile.ZipFile(zip_path, "r") as zf:
        nc_files = [n for n in zf.namelist() if n.endswith(".nc")]
        for nc in nc_files:
            dest = os.path.join(out_dir, os.path.basename(nc))
            if os.path.isfile(dest) and os.path.getsize(dest) > 1e4:
                print(f"  [SKIP] {os.path.basename(nc)} already extracted")
                continue
            print(f"  extracting {os.path.basename(nc)}")
            with zf.open(nc) as src, open(dest, "wb") as dst:
                dst.write(src.read())


for url, fname, min_size in FILES:
    zip_path = os.path.join(OUT_DIR, fname)
    download_file(url, zip_path, min_size)
    extract_nc(zip_path, OUT_DIR)

print("\nEDGAR download complete.")
print("NetCDF files in:", OUT_DIR)
for f in sorted(f for f in os.listdir(OUT_DIR) if f.endswith(".nc")):
    size = os.path.getsize(os.path.join(OUT_DIR, f)) / 1e6
    print(f"  {f}  ({size:.1f} MB)")
