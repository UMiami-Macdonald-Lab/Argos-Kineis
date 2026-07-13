# Argos–Kineis Satellite Telemetry Analysis

This repository contains the R code used for the analyses presented in:

**Exploring Ecological Inferences Across Evolving Satellite Telemetry Platforms: A Comparison of Argos and Kineis Data Streams**

Julia Saltzman, Christopher A. Searcy, Candy K. Real, and Catherine Macdonald

## Overview

The purpose of this study was to evaluate how the transition from the legacy Argos satellite system to the integrated Argos–Kineis constellation affects animal telemetry data and the ecological conclusions derived from those data.

Using satellite telemetry data from sharks tagged in South Florida, the analyses compare Argos- and Kineis-derived locations across several components of movement ecology, including:

* Data yield and transmission frequency
* Temporal spacing and clustering of transmissions
* State-space modeled tracks
* Individual-level home-range estimates
* Species-level activity-space estimates
* Spatial overlap among utilization distributions
* Location quality class distributions
* Positional uncertainty and continuous error estimates
* Behavioral-state classifications

The repository provides a reproducible R workflow for generating the statistical analyses, model diagnostics, figures, and supplemental tables associated with the manuscript.

## Repository Contents

### `MASTER_Argos_Kineis_All_Analyses_and_Supplementals.R`

This is the primary analysis script. It combines the individual analysis workflows into a single ordered pipeline.

The script includes:

1. Importing and cleaning raw Argos and Kineis telemetry records
2. Removing invalid locations, pre-deployment records, and locations occurring on land
3. Removing duplicate timestamps within each shark and satellite system
4. Fitting daily state-space models using `aniMotum`
5. Estimating individual-level 50% and 90% kernel utilization distributions
6. Estimating pooled species-level activity spaces
7. Comparing KDE area among satellite and processing workflows
8. Calculating Jaccard similarity among KDE polygons
9. Comparing transmission intervals and transmission rates
10. Evaluating temporal clustering using the coefficient of variation of transmission intervals
11. Comparing Argos and Kineis location quality class distributions
12. Comparing Kineis continuous error radii with conventional Argos class-based error assumptions
13. Summarizing behavioral-state classifications and matched-day agreement
14. Exporting statistical tables, model diagnostics, figures, spatial files, and supplemental outputs

## Analytical Workflows

Four telemetry workflows are compared throughout the analyses:

* **Argos–Raw:** cleaned locations generated through legacy Argos processing
* **Argos–SSM:** Argos locations regularized with a state-space model
* **Kineis–Raw:** cleaned locations generated through Kineis processing
* **Kineis–SSM:** Kineis locations regularized with a state-space model

State-space modeled tracks are predicted at 24-hour intervals using a correlated random walk process model implemented in `aniMotum`.

Kernel utilization distributions are estimated at the:

* **50% level**, representing core-use areas
* **90% level**, representing broader activity or home-range areas

## Required Input Files

The master script expects the raw Argos and Kineis telemetry files and the associated spatial boundary files to be available locally.

The primary expected inputs are:

```text
Argos-Sharks.xls
Doppler-Sharks.xlsx
Detailed_Florida_State_Boundary.shp
Detailed_Florida_State_Boundary.shx
Detailed_Florida_State_Boundary.dbf
Detailed_Florida_State_Boundary.prj
```

The Florida boundary shapefile is used to identify and remove locations that occur on land.

Input filenames and paths can be modified in the settings section at the beginning of the master script.

## Data Availability

The raw telemetry data used in this study are not included in this repository.

The datasets generated and analyzed during this study are part of ongoing graduate student research and are not publicly available at this time. The data contain fine-scale animal location information from ecologically sensitive nursery and reproductive habitats. Public release of these locations could increase the risk of disturbance or exploitation of critical habitats; therefore, access to the raw telemetry data is restricted.

Data may be made available by the corresponding author upon reasonable request.

Because the raw data are restricted, the full workflow cannot be reproduced directly from this repository without obtaining authorized access to the required telemetry files. The code is provided to document the analytical workflow and support transparency and reproducibility of the methods.

## Software Requirements

The analyses were conducted in R. The master script uses packages including:

```r
dplyr
tidyr
readr
readxl
lubridate
tibble
sf
sp
adehabitatHR
aniMotum
lme4
lmerTest
emmeans
multcomp
multcompView
effectsize
performance
ggplot2
patchwork
purrr
writexl
```

Some analyses may require additional package dependencies that are installed automatically with these packages.

## Running the Analysis

1. Download or clone the repository.
2. Place the authorized input files in the project directory.
3. Open the R project or set the working directory to the repository folder.
4. Confirm that the filenames in the script match the local input filenames.
5. Run:

```r
source("MASTER_Argos_Kineis_All_Analyses_and_Supplementals.R")
```

The script is organized sequentially because later analyses depend on objects and files generated in earlier sections.

## Outputs

The workflow generates several types of outputs, including:

### Cleaned telemetry data

```text
Argos-CLEAN_postdeploy_water.csv
Doppler-CLEAN_postdeploy_water.csv
```

### State-space model outputs

```text
Argos-SSM_daily.csv
Doppler-SSM_daily.csv
SSM_excluded_tracks.csv
```

### Spatial outputs

```text
KDE_polygons_50_90_all_methods.gpkg
KDE_areas_50_90_all_methods.csv
Jaccard_overlap_by_track_level.csv
```

### Statistical outputs

The script exports model summaries, Type III ANOVA tables, estimated marginal means, post-hoc comparisons, effect sizes, intraclass correlation coefficients, and descriptive statistics as `.csv` files.

### Figures

The workflow produces publication-ready `.png` and `.pdf` figures for:

* Transmission spacing and transmission yield
* Individual-level KDE areas
* Species-level activity spaces
* Individual- and species-level spatial overlap
* Location quality class distributions
* Kineis positional error
* Kineis-to-Argos expected-error ratios
* Model diagnostics
* Behavioral-state comparisons

## Reproducibility Notes

Some analytical outputs may vary slightly across software versions because of updates to model-fitting algorithms, package defaults, spatial geometry handling, and random procedures.

For reproducibility, users should record:

```r
sessionInfo()
```

The master script also exports session information when the complete workflow finishes.

## Citation

When using or adapting code from this repository, please cite the associated manuscript:

Saltzman, J., Searcy, C. A., Real, C. K., and Macdonald, C. *Exploring Ecological Inferences Across Evolving Satellite Telemetry Platforms: A Comparison of Argos and Kineis Data Streams.*

The complete publication citation and DOI will be added following publication.

## Contact

Questions regarding the code or requests for access to the underlying data should be directed to the corresponding author. Data requests will be considered based on the proposed use, research objectives, relevant permits, and the need to protect ecologically sensitive animal locations.
