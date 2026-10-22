# Data

## province_rainfall_monthly.csv

Monthly rainfall (mm/month) for the 38 provinces of Indonesia, January 1985 to April 2025 (484 months). The first column is the month (`YYYY-MM`); the other columns are provinces, keyed by their ISO 3166-2:ID code (for example `ID-AC` for Aceh). There are no missing values.

Each value is the mean of all CHIRPS v2.0 pixels (0.05 degree) that touch the province polygon (rasterised with `all_touched=True`), taken from the global monthly CHIRPS v2.0 NetCDF file.

- CHIRPS v2.0: Climate Hazards Center, UC Santa Barbara, https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_monthly/netcdf/
  Funk, C., Peterson, P., Landsfeld, M., et al. (2015). The climate hazards infrared precipitation with stations: a new environmental record for monitoring extremes. *Scientific Data* 2, 150066. https://doi.org/10.1038/sdata.2015.66
- Province polygons: "Provinsi" layer (38 provinces, 2023) of https://github.com/Alf-Anas/batas-administrasi-indonesia. According to that repository, the boundaries come from the GIS portal of the Directorate-General of Population and Civil Registration (Dukcapil), Ministry of Home Affairs of Indonesia. The polygons are not included here.

In the paper, each series is converted to standardised anomalies by calendar month (subtract the mean and divide by the sample standard deviation of that calendar month); see `monthly_anomalies` in `scripts/pattern_tools.py`.

## provinces.csv

One row per province, in the same order as the columns of the rainfall file:

| column | meaning |
|---|---|
| `region_id` | ISO 3166-2:ID code |
| `kode_prov` | two-digit province code (Kepmendagri No. 300.2.2-2138 of 2025) |
| `name` | English province name |
| `lat`, `lon` | polygon centroid in degrees (WGS 84) |

Centroids were computed in the Equal Earth equal-area projection (EPSG:8857) and converted back to latitude and longitude. The geographic weight matrix $W_{\mathrm{geo}}$ of the paper uses the four nearest provinces by great-circle (haversine) distance between these centroids, with weight 1/4 each.

## Terms of use

CHIRPS data are in the public domain. The rainfall series in this folder are released under the Creative Commons Attribution 4.0 International license (CC BY 4.0); please cite Funk et al. (2015) as well. The centroid coordinates are derived from the polygon layer named above, which does not state a license; they are provided only so that the results can be reproduced, and users should credit that source.
