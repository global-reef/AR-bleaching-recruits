### 00. SETUP ####

#### Packages ####

suppressPackageStartupMessages({
  library(readxl)
  library(janitor)
  library(glmmTMB)  
  library(coin)
  library(vegan)
  library(DHARMa)
  library(emmeans)
  library(mclogit)
  library(tidyverse)
})

#### Analysis date ####

analysis_date <- "2026.09.30" # change every time 

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
# colour blind accessible palettes 

## health condition 
# option A
condition_palette2 <- c(
  H = "#2A2B59",
  PBL = "#3B4F7A",
  FBL = "#5E6FA3",
  PRK = "#8A8BC3",
  FRK = "#C5B4D9"
)
# option B 
condition_palette <- c(
  H = "#E8D8B5",
  PBL = "#B9C7B5",
  FBL = "#7FA9A8",
  PRK = "#477F91",
  FRK = "#1F4E6B"
)

#### Shared helper functions ####