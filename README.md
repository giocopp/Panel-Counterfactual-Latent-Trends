# Panel Counterfactual Estimators for Policy Shocks under Latent Trends

**Causal ML Final Project**

Submitted by: Giorgio Coppola  

## Writeup

📝 **Read the writeup [in PDF](writeup/writeup.pdf) or [in HTML](writeup/writeup.html)**

## Overview

Simulation comparing panel causal estimators under latent factor trends (violations of parallel trends). 

**Research question:** How do Two-Ways Fixed Effect (TWFE), Matrix Completion (MC), Synthetic Difference in Difference (SynthDiD), and Triply Robust Panel Estimator (TROP) compare in bias, coverage, and power when outcomes are sparse/binary and units follow latent interactive trends correlated with treatment assignment?

**Application context:** Italy–Libya MoU's effect on migrant mortality in the Central Mediterranean (May 2017 operational onset). The DGP is calibrated to match real IOM data patterns, but the goal is methodological.

**Setup:** 56 grid cells × 35 months, binary outcome (any deaths in cell-month), treatment = proximity to Libya × post-May 2017.

## Reproduce Results  

To reproduce the analysis, run the following commands in R:
```r
source("R/00_setup.R")
source("R/run_analysis.R")
```

To render the writeup (from the project root):
```bash
quarto render writeup/writeup.qmd --to pdf
quarto render writeup/writeup.qmd --to html
```

## Repository structure

```
├── R/
│   ├── 00_setup.R         # Prepare environment
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

## DGP Features

The data generating process includes:

- **Interactive fixed effects**: Unit loadings × time factors
- **Correlated loadings**: Near-Libya units have systematically different loadings, violating parallel trends
- **Sparse binary outcome**: ~6% event rate (any deaths in cell-month)
- **Seasonality**: Calibrated from IOM data

Key parameters:
- `factor_strength`: Controls magnitude of latent trends (0 = none, 0.5 = moderate, 0.9 = strong)
- `loading_correlation`: Controls correlation between loadings and treatment (0 = uncorrelated, 0.5 = moderate)

## Key Findings

1. **MC is best overall** — unbiased, highest power, properly sized
2. **TWFE is surprisingly robust** — only small bias despite correlated loadings
3. **TROP trades power for robustness** — good coverage, lower power (*note: custom, simplified implementation*)
4. **SynthDiD struggles with binary outcomes** — persistent negative bias (~2pp); future work could include testing with continuous outcomes

## Data

Analysis uses a pre-computed calibrated synthetic dataset (`data/calibration_targets.rds`) derived from [IOM Missing Migrants Project](https://missingmigrants.iom.int/) data. These contain summary statistics only (event rates, seasonality indices).

## Dependencies

**Required:**
- tidyverse, fixest, Matrix, patchwork, scales

**Optional (for SynthDiD and parallel processing):**
- synthdid, furrr, progressr

Install all dependencies:
```r
source("R/00_setup.R")
```
