### 00. SETUP ####

#### Packages ####

suppressPackageStartupMessages({
  library(readxl)
  library(janitor)
  library(glmmTMB)
  library(DHARMa)
  library(emmeans)
  library(mclogit)
  library(tidyverse)
})

#### Analysis date ####

analysis_date <- "2026.08.25" # change every time 

#### Directories ####

data_raw_dir   <- "data_raw"
data_clean_dir <- "data_clean"
r_dir          <- "R"
docs_dir       <- "docs"

dir.create(data_clean_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(docs_dir, showWarnings = FALSE, recursive = TRUE)

#### Output structure ####

output_dir <- paste0("Analysis_", analysis_date)

fits_dir   <- file.path(output_dir, "fits")
plots_dir  <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")
stats_dir  <- file.path(output_dir, "stats")
eda_dir    <- file.path(output_dir, "eda")

purrr::walk(
  c(output_dir, fits_dir, plots_dir, tables_dir, stats_dir, eda_dir),
  ~ dir.create(.x, showWarnings = FALSE, recursive = TRUE)
)

#### Plot theme ####

theme_clean <- theme_minimal(base_family = "Arial") +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )


#### Colour palettes ####


#### Shared helper functions ####