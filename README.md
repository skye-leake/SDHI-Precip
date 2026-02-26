# Analysis Pipeline: README

A high-level summary of each notebook in the pipeline, in order of execution.

---

## Preprocessing / Data Extraction

### `01`: Percentile Interval Extraction *(Demo)*
- **Purpose:** Bins 15-minute WRF precipitation into percentile-based intensity bands (intervals); counts occurrences and sums accumulations per band, per grid cell, per water year
- **Inputs:** Demo WRF NetCDF files (50×50 spatial subset), pre-calculated historical percentiles
- **Outputs:** Per-band NetCDF files: occurrence counts (`DIST_*`) and precipitation sums (`Accum_*`)
- **Note:** 2-year, laptop-scale version for method verification only; not scientifically valid at full scale

### `02`: Percentile Interval Extraction *(Full Domain)*
- **Purpose:** Same methodology as `01`, applied to the full WRF domain across all five epochs (HIST, MID4.5, MID8.5, EOC4.5, EOC8.5)
- **Inputs:** Full-domain WRF NetCDF files, historical percentiles dataset
- **Outputs:** Per-band NetCDF files: occurrence counts and precipitation sums, by scenario/year/season

---

## IDF Curve Analysis

### `03`: IDF Curve Analysis *(Demo)*
- **Purpose:** Derives Intensity-Duration-Frequency (IDF) curves for target metro areas from annual precipitation maxima; fits Gumbel distributions and quantifies uncertainty via bootstrapping
- **Inputs:** Demo WRF precipitation data (50×50 subset, 2 years); 7×7 grid cell metro domains
- **Outputs:** IDF tables and diagnostic plots
- **Note:** Demo sample is artificially augmented via duplication; outputs demonstrate methodology only: not scientifically valid

### `04`: IDF Curve Analysis *(Full Domain)*
- **Purpose:** Same methodology as `03` applied epoch by epoch across all metro areas
- **Inputs:** Full WRF precipitation data; 21×21 grid cell (~80×80 km) metro domains
- **Outputs:** IDF intensity tables by city/scenario/return period/duration; bootstrap confidence interval NetCDF files; diagnostic plots

---

## Figure Generation

### `30`: Single-Panel Map Figure
- **Purpose:** Generates single-panel maps of precipitation fields with NCA region outlines and cartographic formatting
- **Inputs:** Zarr/NetCDF precipitation data
- **Outputs:** Single-panel map figures

### `31`: Multipanel Figure Data Preparation
- **Purpose:** Pre-computes all data needed for multipanel figures: calculates departures/percentages from historical baseline and runs Mann-Whitney U + FDR statistical tests; saves results as intermediate NetCDF files
- **Inputs:** Zarr data for historical and future scenarios
- **Outputs:** Intermediate NetCDF files (`multipanel_percentile/`, `multipanel_{RCP}/`)

### `32`: Multipanel Figure Generation
- **Purpose:** Loads pre-computed statistics from `31` and assembles the final publication-ready multipanel precipitation figure
- **Inputs:** Intermediate NetCDF files from `31`
- **Outputs:** High-resolution multipanel figure (PNG/PDF)

---

## Percentile Band Heatmaps

### `36`: Annual Precipitation Percentile Heatmap
- **Purpose:** Calculates regional percent changes in occurrence and accumulation by percentile band for annual aggregations across all NCA regions; generates summary bar charts and the primary heatmap figure
- **Inputs:** Zarr stores (annual): occurrence counts and precipitation sums, all scenarios/epochs
- **Outputs:** Bar chart figures (intermediate QC); Figure 2: 2×2 heatmap grid (RCP4.5/8.5 × Occurrences/Accumulation, MID and EOC)

### `37`: Seasonal Precipitation Percentile Heatmap
- **Purpose:** Same as `36`, stratified by season (DJF, MAM, JJA, SON)
- **Inputs:** Zarr stores (seasonal subsets): same structure as `36`
- **Outputs:** Seasonal bar chart figures; 2×2 heatmap grids per season

---

## Supplementary Tables

### `40`: Bootstrap Confidence Interval Table Generation
- **Purpose:** Extracts and formats bootstrap uncertainty estimates from IDF results into publication-ready LaTeX tables; also computes relative uncertainty metrics by duration and scenario
- **Inputs:** Bootstrap result NetCDF files from `04` (1000-iteration CI bounds, 441 cells per metro)
- **Outputs:** Tables S1–S3: median IDF intensity ± CI half-width for 15-min, 1-hr, and 24-hr durations across 8 cities, 5 scenarios, and 6 return periods; diagnostic uncertainty figures