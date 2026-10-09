# SDHI-precipitation

Short-Duration High-Intensity (SDHI) Precipitation from a Convection-Permitting Regional Climate Model (CP-RCM)

> **Status (January 2026):** Data and scripts are preliminary and under peer review, provided for reviewers. Use with caution; they may change.

**Please cite:** Gensini, V. A., A. M. Haberlie, and W. S. Ashley, 2023: Convection-permitting simulations of historical and possible future climate over the contiguous United States. _Climate Dynamics_, 60, 109–126.

---

## Setup

1. Download the supporting data (~20 GB) with [`00_download_data.ipynb`](https://github.com/skye-leake/SDHI-Precip/blob/main/00_download_data.ipynb)
2. Extract: `tar -xzf SDHI_Precip_supporting_data.tar.gz`
3. `cd SDHI_Precip_supporting_data`
4. Launch `jupyter notebook` (or `jupyter lab`) **from this directory**; all notebooks use paths relative to it
5. If running elsewhere, set `DATA_ROOT` in each notebook's config block

**Epochs:** HIST (1990–2005), MID4.5/MID8.5 (2040–2055), EOC4.5/EOC8.5 (2085–2100); 15 water years each.

---

## Pipeline

Notebooks marked _demo_ run on a laptop-scale subset (50×50 cells, 2 years) to verify methods only; results are not scientifically valid.

| #   | Notebook                       | Purpose                                                                                                         | Output                                                                        |
| --- | ------------------------------ | --------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- |
| 00  | Data download                  | Fetch supporting data archive                                                                                   | `SDHI_Precip_supporting_data.tar.gz`                                          |
| 01  | Percentile extraction _(demo)_ | Bin 15-min precipitation into percentile bands; count occurrences and sum accumulation per cell, per water year | `DIST_*`, `Accum_*` NetCDF                                                    |
| 02  | Percentile extraction _(full)_ | As 01, full domain, all epochs                                                                                  | NetCDF by scenario/year/season                                                |
| 03  | IDF curves _(demo)_            | Annual maxima → Gumbel fit → bootstrap CIs (7×7 metro domains)                                                  | IDF tables, diagnostics                                                       |
| 04  | IDF curves _(full)_            | As 03, 21×21 cell (~80 km) metro domains, all epochs                                                            | IDF tables, bootstrap CI NetCDF                                               |
| 30  | Single-panel map               | Precipitation maps with NCA region outlines                                                                     | Map figures                                                                   |
| 31  | Multipanel data prep           | Departures/percent change vs. HIST; Mann-Whitney U + FDR                                                        | `multipanel_percentile/`, `multipanel_{RCP}/`                                 |
| 32  | Multipanel figure              | Assemble publication figure from 31                                                                             | PNG/PDF                                                                       |
| 35  | Departure small multiples      | Grid of annual and seasonal departure-from-HIST maps                                                            | Figure 3                                                                      |
| 36  | Annual heatmap                 | Regional % change in occurrence/accumulation by band (NCA regions)                                              | Figure 2 (2×2 heatmap), QC bar charts                                         |
| 37  | Seasonal heatmap               | As 36, by DJF/MAM/JJA/SON                                                                                       | Seasonal heatmaps, bar charts                                                 |
| 40  | Bootstrap CI tables            | Median IDF intensity ± CI half-width across 441 cells                                                           | Tables S1–S3 (15-min, 1-hr, 24-hr; 8 cities × 5 scenarios × 6 return periods) |

**Dependencies:** 31 → 32; 02 → 30, 31, 35, 36, 37; 04 → 40.
