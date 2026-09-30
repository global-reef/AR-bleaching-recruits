# Coral bleaching in juvenile recruits and mature reef corals

Analysis repository for a collaborative study between Global Reef and the Aow Thai Marine Ecology Centre (ATMEC), comparing bleaching responses of juvenile corals recruiting to artificial reefs (AR) and mature corals on natural reefs (NR) during the 2024 mass bleaching event in Thailand.

## Overview

This study examines whether juvenile coral recruits on artificial reef structures differ in bleaching response from mature corals on adjacent natural reefs.

Coral condition was quantified from CPCE image points at:

- Koh Tao: Mango Bay and Tanote Bay, sampled across Start, Middle, Late, and Post bleaching periods
- Rayong: sampled once during the Middle bleaching period

Artificial reef corals represent juvenile recruits, while natural reef corals represent mature colonies.

Five coral-condition states are retained:

- H: Healthy
- PBL: Partially bleached
- FBL: Fully bleached
- PRK: Partially recently killed
- FRK: Fully recently killed

These states are treated as nominal rather than ordinal.

## Research questions

1. At the temporally matched Middle bleaching period, does coral condition differ between juvenile AR and mature NR corals, and does this difference vary between Koh Tao and Rayong?

2. Across the bleaching event at Koh Tao, do AR and NR coral-condition trajectories differ through time, and does this pattern vary between sites?

3. Within coral genera, does five-state condition composition differ between AR and NR at the Middle bleaching period?

4. Within coral genera at Koh Tao, do temporal condition trajectories differ between AR and NR?

## Analysis framework

The analysis uses two complementary approaches.

### Q1-Q2: population-level coral condition

Overall coral-condition differences are analysed using frequentist multinomial mixed models implemented with `mclogit::mblogit`.

CPCE points are treated as subsamples within quadrats, with `quadrat_id` included as a random intercept to account for clustering.

- Q1: `condition ~ reef_type * location`
- Q2: `condition ~ reef_type * bleaching_period * site`

Estimated condition probabilities and AR-NR contrasts are obtained using `emmeans`.

### Q3-Q4: genus-specific condition composition

Genus-specific analyses use quadrat-level proportions of the five condition states.

Bray-Curtis PERMANOVA is used to test:

- Q3: AR-NR differences within genus at the Middle bleaching period
- Q4: reef type × bleaching period interactions within genus at Koh Tao

Permutation tests use 99,999 permutations with site-restricted permutations where appropriate.

PERMDISP is used to assess multivariate dispersion, and false discovery rate is controlled using Benjamini-Hochberg adjustment.

## Repository structure

```text
00_SETUP.R
00.1_CPCE-reshape.R
01_CLEAN.R
02_EXPLORE.R
03_MODEL.R
04_RESULTS.R
05_PLOTS.R
99_SCRATCH.R

data_raw/
data_clean/
docs/

Analysis_YYYY.MM.DD/
├── eda/
├── fits/
├── plots/
├── stats/
└── tables/
```

### Scripts

`00_SETUP.R`  
Loads packages, defines directories, creates dated analysis folders, and stores shared plotting settings and colour palettes.

`00.1_CPCE-reshape.R`  
Reshapes the original AR and NR CPCE spreadsheets into a common long-format dataset.

`01_CLEAN.R`  
Standardises substrate and health codes, assigns coral taxonomy, defines bleaching periods and sampling identifiers, and creates the analysis-ready five-state coral-condition response.

`02_EXPLORE.R`  
Evaluates sampling balance, replication, missingness, genus overlap, condition-state sparsity, and model estimability.

`03_MODEL.R`  
Fits the four primary analyses and saves model fits, diagnostics, permutation-test results, and genus-level support tables.

`04_RESULTS.R`  
Generates estimated probabilities, post-hoc contrasts, adjusted genus-level results, and manuscript/supplementary tables.

`05_PLOTS.R`  
Produces the four final manuscript figures in TIFF and PNG formats.

`99_SCRATCH.R`  
Temporary exploratory code not used in the reproducible analysis pipeline.

## Output structure

Each analysis run is written to a dated directory:

`Analysis_YYYY.MM.DD/`

with:

- `eda/`: exploratory and model-support summaries
- `fits/`: fitted model and permutation objects
- `stats/`: detailed analytical outputs and diagnostics
- `tables/`: manuscript and supplementary tables
- `plots/`: manuscript figures

This preserves previous analysis states while keeping the current workflow reproducible.

## Reproducibility

Analyses are run sequentially from:

1. `00_SETUP.R`
2. `00.1_CPCE-reshape.R`
3. `01_CLEAN.R`
4. `02_EXPLORE.R`
5. `03_MODEL.R`
6. `04_RESULTS.R`
7. `05_PLOTS.R`

Permutation analyses use:

```r
set.seed(42)
n_perm <- 99999
```

Core packages include:

- tidyverse
- mclogit
- emmeans
- vegan
- ggalluvial
- readxl
- janitor

## Interpretation notes

Artificial reef and natural reef comparisons are also juvenile-mature comparisons and therefore reef type and life stage cannot be separated statistically.

Exact quadrats and colonies were not repeatedly tracked through time. Temporal patterns therefore represent population-level shifts in coral condition rather than individual colony transitions.

Rayong was sampled only during the Middle bleaching period. All Koh Tao-Rayong comparisons are therefore restricted to the corresponding Middle period.


## Collaboration
This project is a collaboration between:
- Global Reef, Koh Tao, Thailand
- Aow Thai Marine Ecology Centre (ATMEC), Rayong, Thailand

## License
This project is private and not licensed for redistribution. For collaboration inquiries, please contact scarlett@global-reef.com.
Affiliation: Global Reef, Koh Tao, Thailand