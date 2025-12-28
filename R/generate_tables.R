#!/usr/bin/env Rscript
# R/generate_tables.R
#
# Reproduces Table 1 for main, and Tables S1, and S2 for Supplement.
#
# Outputs:
# - `results/table_comparison_baseline_mixed.tex`
# - `results/table_comparison_baseline_exp.tex`
# - `results/table_comparison_baseline_rational.tex`
#

suppressPackageStartupMessages({
  library(dplyr)
  library(reticulate)
  library(knitr)
  library(kableExtra)
})

project_root <- getwd() # assumes running from root

# set venv and load config.py
use_virtualenv(file.path(project_root, ".venv"), required = TRUE)
py_run_string(sprintf("import sys; sys.path.append('%s')", file.path(project_root, "src")))
config <- import("config")
RESULT <- config$RESULT

# load median NSSR distances
nssr_path <- file.path(project_root, "results", "model_comparison", 
                       paste0("median_nssr_distances_", RESULT, ".csv"))
df_nssr <- read.csv(nssr_path, check.names = FALSE)


generate_comparison_tex <- function(variant_display_name, variant_code, out_path) {
  
  # load model selection summary results for this variant 
  bf_filename <- paste0("baseline_vs_", variant_code, "_summary_", RESULT, ".csv")
  bf_path <- file.path(project_root, "results", "model_selection", bf_filename)  
    
  df_bf <- read.csv(bf_path, check.names = FALSE)

  # merge nssr and bf by location
  df_merged <- df_nssr %>%
    inner_join(df_bf, by = "Location")

  # find bayes factor column
  bf_col_name <- names(df_merged)[grep("^BF", names(df_merged))]

  # output columns:
  # state_abbr, NSSR_baseline, NSSR_variant, BF_eq, Evidence
  df_formatted <- df_merged %>%
    select(
      State    = Abbr, 
      NSSE_B   = Baseline,
      NSSE_V   = all_of(variant_display_name),
      BF_val   = all_of(bf_col_name), 
      Evidence = Interpretation
    ) %>%
    arrange(desc(as.numeric(BF_val))) %>%
    rowwise() %>%
    # bold the lowest NSSR in a row
    mutate(
      bold_B = as.numeric(NSSE_B) < as.numeric(NSSE_V),
      
      # baseline
      NSSE_B_str = sprintf("%.4f", as.numeric(NSSE_B)),
      NSSE_B_tex = ifelse(bold_B, paste0("\\textbf{", NSSE_B_str, "}"), NSSE_B_str),
      
      # behavior variant
      NSSE_V_str = sprintf("%.4f", as.numeric(NSSE_V)),
      NSSE_V_tex = ifelse(!bold_B, paste0("\\textbf{", NSSE_V_str, "}"), NSSE_V_str),
      
      # bayes Ffctor, handle infinity
      BF_tex = ifelse(is.infinite(as.numeric(BF_val)), 
                      "$\\infty$", 
                      sprintf("%.2f", as.numeric(BF_val)))
    ) %>%
    ungroup() %>%
    select(State, NSSE_B = NSSE_B_tex, NSSE_V = NSSE_V_tex, BF_eq = BF_tex, Evidence)

  # map variants names to letters for subscripts
  subscript <- switch(variant_code, "mixed" = "M", "exp" = "E", "rational" = "R")
  col_names <- c(
    "State", 
    "Median $NSSE_{B}$", 
    paste0("Median $NSSE_{", subscript, "}$"), 
    "$BF_{eq}$", 
    "Evidence"
  )

  # generate tabular without the table environment
  k_latex <- df_formatted %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      escape = FALSE, 
      col.names = col_names,
      align = "lrrrl",
      table.envir = NULL
    )

  # clean up latex
  k_latex_str <- as.character(k_latex)
  
  # export to tex file
  writeLines(k_latex_str, out_path)
  cat("Generated table:", out_path, "\n")
}

# save tables to results 
res_dir <- file.path(project_root, "results")

# Main Paper Table 1 (Mixed)
generate_comparison_tex("Behavioral (Mixed)", "mixed", 
  file.path(res_dir, "table_comparison_baseline_mixed.tex"))

# Supplement Table S1 (Exponential)
generate_comparison_tex("Behavioral (Exponential)", "exp", 
  file.path(res_dir, "table_comparison_baseline_exp.tex"))

# Supplement Table S2 (Rational)
generate_comparison_tex("Behavioral (Rational)", "rational", 
  file.path(res_dir, "table_comparison_baseline_rational.tex"))

