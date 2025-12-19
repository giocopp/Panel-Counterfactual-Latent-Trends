# Panel Counterfactual Estimators for Policy Shocks under Latent Trends
**Causal ML Final Project (simulation tutorial + reproducible code)**

This project is a **simulation-based power analysis** comparing modern panel causal estimators under **latent interactive trends** (violations of parallel trends). The motivating application is the **Italy–Libya MoU** and the operational shift around **May 2017**, but the goal is methodological: *when do different panel counterfactual estimators work well?*

## Research question (course project version)
> In a spatial unit×month panel, how do TWFE, Matrix Completion, Synthetic DiD, and TROP compare in bias, coverage, RMSE, and power when the data exhibit latent factor trends and sparse outcomes?

## Setup 
- **Units:** 1°×1° grid cells in the Central Mediterranean (56 cells by default)
- **Time:** months (Apr 2015–Feb 2018 by default)
- **Outcome (primary):** `Y_any = 1{deaths_observed > 0}` (deadly-event indicator per cell-month)
- **Treatment exposure:** `near_libya = 1{distance_to_tripoli ≤ 200 km}`
- **Post period:** `post_treat = 1{t ≥ May 2017}`
- **Treatment indicator:** `D_it = near_libya_i × post_treat_t`
- **Estimand:** ATT on `Y_any` for treated units in the post period

## Methods compared
- **TWFE** (baseline DiD with unit and time fixed effects)
- **Matrix Completion (MC)** (low-rank imputation of the untreated counterfactual)
- **Synthetic DiD (SDID)** (weights + DiD)
- **TROP** (hybrid weighting + low-rank + TWFE; designed to be robust under certain misspecifications)

## Calibration using real data
The `R/run_analysis.R` will:
1) aggregate incidents to the same **grid×month** panel (Central Mediterranean only),  
2) compute the **pre-period target event rate** for `Y_any`,  
3) **calibrate `baseline_mortality`** so the simulated pre-period `mean(Y_any)` matches the IOM target, and  
4) replace the default sinusoidal seasonality with **month-of-year multipliers** learned from IOM.

## Repository structure

```
causal_ml_project/
├── R/
│   ├── 00_setup.R         
│   ├── 01_dgp.R           
│   ├── 02_estimators.R    
│   ├── 03_simulation.R    
│   ├── 04_calibration.R   
│   └── run_analysis.R     
├── data/
├── figures/
│   ├── dgp_overview.pdf
│   ├── power_curves.pdf
│   └── main_results.pdf
├── output/
│   ├── power_summary.csv
│   └── table1_main_results.csv
├── writeup/
│   └── writeup.qmd
└── README.md
```

## Run instructions
From the project root:

```r
source("R/00_setup.R")
source("R/run_analysis.R")
```

Then render the writeup:

```bash
quarto render writeup.qmd
```

Outputs are written to `figures/` and `output/`.
