#!/usr/bin/env python3
"""
Download ERA5 pressure-level and single-level data for WRF-GRK
Period: 2015-07-25 to 2015-11-30 (met spin-up starts Jul 25)
Domain: Indonesia + buffer [16N, 88E, -16S, 147E]
"""
import cdsapi
import os

# CDS credentials (do not modify ~/.cdsapirc)
CDS_URL = "https://cds.climate.copernicus.eu/api"
CDS_KEY = "cdc594b9-8754-417f-83f8-97259710653d"

OUT_DIR = "/home/igrk/WRF-GRK/rawdata/era5"
os.makedirs(OUT_DIR, exist_ok=True)

AREA = [16, 88, -16, 147]  # N, W, S, E

PRESSURE_VARS = [
    "geopotential", "temperature", "u_component_of_wind",
    "v_component_of_wind", "vertical_velocity", "specific_humidity",
    "fraction_of_cloud_cover", "specific_cloud_ice_water_content",
    "specific_cloud_liquid_water_content",
]

PRESSURE_LEVELS = [
    "1", "2", "3", "5", "7", "10", "20", "30", "50", "70",
    "100", "125", "150", "175", "200", "225", "250", "300",
    "350", "400", "450", "500", "550", "600", "650", "700",
    "750", "775", "800", "825", "850", "875", "900", "925",
    "950", "975", "1000",
]

SINGLE_VARS = [
    "10m_u_component_of_wind", "10m_v_component_of_wind",
    "2m_dewpoint_temperature", "2m_temperature",
    "land_sea_mask", "mean_sea_level_pressure",
    "sea_ice_cover", "sea_surface_temperature",
    "skin_temperature", "snow_depth",
    "soil_temperature_level_1", "soil_temperature_level_2",
    "soil_temperature_level_3", "soil_temperature_level_4",
    "surface_pressure", "volumetric_soil_water_layer_1",
    "volumetric_soil_water_layer_2", "volumetric_soil_water_layer_3",
    "volumetric_soil_water_layer_4",
]

MONTHS = ["07", "08", "09", "10", "11"]
YEAR = "2015"

c = cdsapi.Client(url=CDS_URL, key=CDS_KEY)

for month in MONTHS:
    # ---- Pressure levels ----
    out_pl = os.path.join(OUT_DIR, f"era5_pl_{YEAR}{month}.grib")
    if os.path.exists(out_pl) and os.path.getsize(out_pl) > 1e6:
        print(f"[SKIP] {out_pl} already exists")
    else:
        print(f"Downloading ERA5 pressure levels: {YEAR}-{month}...")
        c.retrieve(
            "reanalysis-era5-pressure-levels",
            {
                "product_type": "reanalysis",
                "variable": PRESSURE_VARS,
                "pressure_level": PRESSURE_LEVELS,
                "year": YEAR,
                "month": month,
                "day": [f"{d:02d}" for d in range(1, 32)],
                "time": [f"{h:02d}:00" for h in range(0, 24, 6)],
                "area": AREA,
                "format": "grib",
            },
            out_pl,
        )
        print(f"  -> {out_pl}")

    # ---- Single levels ----
    out_sl = os.path.join(OUT_DIR, f"era5_sl_{YEAR}{month}.grib")
    if os.path.exists(out_sl) and os.path.getsize(out_sl) > 1e6:
        print(f"[SKIP] {out_sl} already exists")
    else:
        print(f"Downloading ERA5 single levels: {YEAR}-{month}...")
        c.retrieve(
            "reanalysis-era5-single-levels",
            {
                "product_type": "reanalysis",
                "variable": SINGLE_VARS,
                "year": YEAR,
                "month": month,
                "day": [f"{d:02d}" for d in range(1, 32)],
                "time": [f"{h:02d}:00" for h in range(0, 24, 6)],
                "area": AREA,
                "format": "grib",
            },
            out_sl,
        )
        print(f"  -> {out_sl}")

print("\nERA5 download complete.")
