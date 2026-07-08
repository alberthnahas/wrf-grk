#!/usr/bin/env python3
"""
Download MODIS MOD09A1 v061 (8-day surface reflectance) and MCD12Q1 v061 (land cover) for VPRM
Period: 2015 (full year for EVI/LSWI min/max)
Coverage: Indonesia bounding box 95E-141E, 11S-6N

Uses NASA CMR API + LP DAAC Earthdata Cloud
Token read from rawdata/modis/nasa_credentials.txt (line: token: <value>)
"""
import os
import sys
import requests
from tqdm import tqdm

OUT_DIR = "/home/igrk/WRF-GRK/rawdata/modis"
CREDS_FILE = os.path.join(OUT_DIR, "nasa_credentials.txt")
os.makedirs(OUT_DIR, exist_ok=True)

# Read token from credentials file
TOKEN = os.environ.get("EARTHDATA_TOKEN", "")
if not TOKEN and os.path.exists(CREDS_FILE):
    with open(CREDS_FILE) as f:
        for line in f:
            if line.startswith("token:"):
                TOKEN = line.split(":", 1)[1].strip()
                break

if not TOKEN:
    print(f"ERROR: No NASA Earthdata token. Add 'token: ...' to {CREDS_FILE}")
    sys.exit(1)

print("Using NASA Earthdata token.")
HEADERS = {"Authorization": f"Bearer {TOKEN}"}
CMR_BASE = "https://cmr.earthdata.nasa.gov/search/granules.json"

# Indonesia bounding box: west, south, east, north
BBOX = "95,-11,141,6"
YEAR = 2015


def cmr_search_all(short_name, version, temporal, bbox, page_size=200):
    """Return all HDF download URLs for product within bbox and temporal range."""
    params = {
        "short_name": short_name,
        "version": version,
        "temporal[]": temporal,
        "bounding_box": bbox,
        "page_size": page_size,
    }
    urls = []
    page = 1
    while True:
        params["page_num"] = page
        r = requests.get(CMR_BASE, params=params, headers=HEADERS, timeout=30)
        r.raise_for_status()
        entries = r.json().get("feed", {}).get("entry", [])
        if not entries:
            break
        for entry in entries:
            for link in entry.get("links", []):
                href = link.get("href", "")
                if href.endswith(".hdf") and "lpdaac.earthdatacloud" in href:
                    urls.append(href)
                    break
        if len(entries) < page_size:
            break
        page += 1
    return urls


def download_file(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 1e5:
        print(f"[SKIP] {os.path.basename(dest)}")
        return True
    try:
        r = requests.get(url, headers=HEADERS, stream=True, timeout=600,
                         allow_redirects=True)
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(dest, "wb") as f, tqdm(total=total, unit="B", unit_scale=True,
                                          desc=os.path.basename(dest)[:45], leave=False) as bar:
            for chunk in r.iter_content(65536):
                f.write(chunk)
                bar.update(len(chunk))
        return True
    except Exception as e:
        print(f"  ERROR {os.path.basename(dest)}: {e}")
        if os.path.exists(dest):
            os.remove(dest)
        return False


# ---- MOD09A1.061: 8-day composites, full year 2015 ----
print(f"\n=== MOD09A1.061 (8-day surface reflectance, {YEAR}) ===")
urls = cmr_search_all("MOD09A1", "061", f"{YEAR}-01-01,{YEAR}-12-31", BBOX)
print(f"Found {len(urls)} granules")
for url in urls:
    fname = os.path.basename(url)
    download_file(url, os.path.join(OUT_DIR, fname))

# ---- MCD12Q1.061: Annual land cover 2015 ----
print(f"\n=== MCD12Q1.061 (annual land cover {YEAR}) ===")
urls_lc = cmr_search_all("MCD12Q1", "061", f"{YEAR}-01-01,{YEAR}-12-31", BBOX, page_size=50)
print(f"Found {len(urls_lc)} granules")
for url in urls_lc:
    fname = os.path.basename(url)
    download_file(url, os.path.join(OUT_DIR, fname))

print("\nMODIS download complete.")
hdf_files = [f for f in os.listdir(OUT_DIR) if f.endswith(".hdf")]
print(f"  {len(hdf_files)} HDF files in {OUT_DIR}")
