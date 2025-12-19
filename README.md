# Panel Counterfactual Estimators for Policy Shocks under Latent Trends

**Causal ML Final Project**

Submitted by: Giorgio Coppola  

## Writeup

📄 **Read the writeup [in PDF](writeup/writeup.pdf) or [in HTML](writeup/writeup.html)**

## Overview

Simulation-based comparison of panel causal estimators under latent factor trends (violations of parallel trends). 

**Research question:** How do TWFE, Matrix Completion, Synthetic DiD, and TROP compare in bias, coverage, and power when outcomes are sparse and units follow latent interactive trends?

**Application context:** Italy–Libya MoU's effect on migrant mortality in the Central Mediterranean (May 2017 operational onset). The DGP is calibrated to match real IOM data patterns, but the goal is methodological.

**Setup:** 56 grid cells × 35 months, binary outcome (any deaths in cell-month), treatment = proximity to Libya × post-May 2017.

## Reproduce
```r
source("R/run_analysis.R")
```
```bash
quarto render writeup/writeup.qmd
```

## Repository structure
```
├── R/
│   ├── 01_dgp.R           # Data generating process
│   ├── 02_calibration.R   # Calibration from IOM data
│   ├── 03_estimators.R    # TWFE, MC, SDID, TROP
│   ├── 04_simulation.R    # Power and scenario analysis
│   └── run_analysis.R     # Main script
├── data/
│   └── calibration_targets.rds
├── figures/               # Output plots
├── output/                # Summary tables
└── writeup/
    └── writeup.qmd        # Writeup in various formats
```

## Data

Analysis uses pre-computed calibration targets (`data/calibration_targets.rds`) derived from [IOM Missing Migrants Project](https://missingmigrants.iom.int/) data. These contain summary statistics only (event rates, seasonality indices)—no raw incident data. Fully reproducible without access to the original dataset.