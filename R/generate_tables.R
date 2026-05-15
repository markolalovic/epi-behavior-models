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
# ordered by decreasing Bayes factor in favor of the behavioral model
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
  # df_merged <- df_nssr %>%
  #   inner_join(df_bf, by = "Location")
  # restrict to current selected locations
  selected_locs <- unlist(config$LOCATIONS)

  df_nssr_selected <- df_nssr %>%
    filter(Abbr %in% selected_locs)

  missing_nssr <- setdiff(selected_locs, df_nssr_selected$Abbr)
  if (length(missing_nssr) > 0) {
    stop(paste("Missing NSSE rows for:", paste(missing_nssr, collapse = ", ")))
  }

  missing_bf <- setdiff(df_nssr_selected$Location, df_bf$Location)
  if (length(missing_bf) > 0) {
    stop(paste("Missing model-selection rows for:", paste(missing_bf, collapse = ", ")))
  }

  # merge nssr and bf by full location name
  df_merged <- df_nssr_selected %>%
    inner_join(df_bf, by = "Location")

  if (nrow(df_merged) != length(selected_locs)) {
    stop("Merged table does not contain exactly the selected locations.")
  }  

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

generate_excluded_locations_tex <- function(out_path) {
  
  selection_path <- file.path(project_root, "results", "state_selection_summary.csv")
  
  if (!file.exists(selection_path)) {
    stop(paste("Missing state-selection summary:", selection_path))
  }
  
  df_sel <- read.csv(selection_path, check.names = FALSE)
  
  selected_locs <- unlist(config$LOCATIONS)
  excluded_locs <- unlist(config$LOCATIONS_BAD)
  
  accidentally_selected <- intersect(selected_locs, excluded_locs)
  if (length(accidentally_selected) > 0) {
    stop(paste(
      "Locations appear in both LOCATIONS and LOCATIONS_BAD:",
      paste(accidentally_selected, collapse = ", ")
    ))
  }
  
  df_excluded <- df_sel %>%
    filter(location %in% excluded_locs)
  
  missing_rows <- setdiff(excluded_locs, df_excluded$location)
  if (length(missing_rows) > 0) {
    stop(paste(
      "Missing state-selection rows for:",
      paste(missing_rows, collapse = ", ")
    ))
  }
  
  df_table <- df_excluded %>%
    mutate(
      # Use the final config exclusion list as authoritative.
      objective_status = "excluded",
      
      # CO was excluded because of a reporting/data anomaly.
      reason = case_when(
        location == "CO" ~ "reporting anomaly",
        TRUE ~ reason
      ),
      
      primary_reason = case_when(
        reason == "insufficient mortality signal" ~ "Low signal",
        reason == "no completed mortality wave" ~ "Incomplete wave",
        reason == "non-coherent mortality wave" ~ "Non-coherent wave",
        reason == "reporting anomaly" ~ "Reporting anomaly",
        TRUE ~ reason
      ),
      
      reason_order = case_when(
        primary_reason == "Low signal" ~ 1,
        primary_reason == "Incomplete wave" ~ 2,
        primary_reason == "Non-coherent wave" ~ 3,
        primary_reason == "Reporting anomaly" ~ 4,
        TRUE ~ 5
      ),
      
      diagnostic_note = case_when(
        location == "CO" ~ "reporting discontinuity",
        
        primary_reason == "Low signal" &
          as.numeric(peak) < 10 ~ "peak $<10$",
        
        primary_reason == "Low signal" &
          as.numeric(total_deaths) < 50 ~ "$\\sum_t y(t)<50$",
        
        primary_reason == "Incomplete wave" &
          as.numeric(days_after_peak) < 21 ~ "peak near window end",
        
        primary_reason == "Incomplete wave" &
          as.numeric(tail_to_peak) > 0.40 ~ "tail/peak $>0.40$",
        
        primary_reason == "Non-coherent wave" &
          as.numeric(high_blocks) > 1 ~ "multiple high-mortality blocks",
        
        primary_reason == "Non-coherent wave" &
          as.numeric(n_outlier_days) > 3 ~ "large local deviations",
        
        TRUE ~ ""
      ),
      
      peak_fmt = sprintf("%.2f", as.numeric(peak)),
      tail_to_peak_fmt = sprintf("%.2f", as.numeric(tail_to_peak)),
      
      # Block count is only relevant for the non-coherent-wave exclusion.
      blocks_fmt = case_when(
        primary_reason == "Non-coherent wave" ~ as.character(high_blocks),
        TRUE ~ "--"
      )
    ) %>%
    arrange(reason_order, location) %>%
    select(
      State = location,
      `Primary reason` = primary_reason,
      Peak = peak_fmt,
      `Tail/peak` = tail_to_peak_fmt,
      Blocks = blocks_fmt,
      Note = diagnostic_note
    )
  
  col_names <- c(
    "State",
    "Primary reason",
    "$\\max_t y(t)$",
    "Tail/peak",
    "Blocks",
    "Diagnostic note"
  )
  
  k_latex <- df_table %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      col.names = col_names,
      align = "llrrll",
      table.envir = NULL
    )
  
  k_latex_str <- as.character(k_latex)
  
  # Remove row spacing inserted by kable/booktabs, if present.
  k_latex_str <- gsub("\\\\addlinespace\n?", "", k_latex_str)
  
  writeLines(k_latex_str, out_path)
  cat("Generated table:", out_path, "\n")
}

# NOTE: new table for sensitivity analysis
generate_sensitivity_mixed_tex <- function(out_path) {
  
  mixed_col <- "Behavioral (Mixed)"
  
  selected_locs <- unlist(config$LOCATIONS)
  excluded_locs <- unlist(config$LOCATIONS_BAD)
  
  # model-selection summaries issues:
  # - RESULT 1 has full location names
  # - RESULT 2 has abbreviations and different column names
  bf1_path <- file.path(
    project_root, "results", "model_selection",
    "baseline_vs_mixed_summary_1.csv"
  )
  bf2_path <- file.path(
    project_root, "results", "model_selection",
    "baseline_vs_mixed_summary_2.csv"
  )
  
  bf1 <- read.csv(bf1_path, check.names = FALSE)
  bf2 <- read.csv(bf2_path, check.names = FALSE)
  
  bf1_col <- names(bf1)[grep("^BF", names(bf1))]
  
  # standardize RESULT 1 format: Location is full name
  bf1_std <- bf1 %>%
    transmute(
      Location = Location,
      BF_eq = as.numeric(.data[[bf1_col]]),
      Evidence = Interpretation
    ) %>%
    inner_join(
      df_nssr %>% select(Location, Abbr),
      by = "Location"
    )
  
  # Standardize RESULT 2 format: location is abbreviation
  bf2_std <- bf2 %>%
    filter(behavior_variant == "behavior_mixed") %>%
    transmute(
      Abbr = location,
      BF_eq = as.numeric(BF_eq),
      Evidence = Interpretation
    ) %>%
    inner_join(
      df_nssr %>% select(Location, Abbr),
      by = "Abbr"
    )
  
  # Combine model-selection results for all 51 locations
  # if a location appears twice, keep the first occurrence
  bf_all <- bind_rows(bf1_std, bf2_std) %>%
    distinct(Abbr, .keep_all = TRUE)
  
  df_table <- df_nssr %>%
    select(
      Location,
      Abbr,
      Baseline,
      all_of(mixed_col)
    ) %>%
    inner_join(
      bf_all %>% select(Abbr, BF_eq, Evidence),
      by = "Abbr"
    ) %>%
    mutate(
      Set = case_when(
        Abbr %in% selected_locs ~ "Included",
        Abbr %in% excluded_locs ~ "Excluded",
        TRUE ~ "Other"
      ),
      Set_order = case_when(
        Set == "Included" ~ 1,
        Set == "Excluded" ~ 2,
        TRUE ~ 3
      ),
      Avg_NSSE = (as.numeric(Baseline) + as.numeric(.data[[mixed_col]])) / 2
    ) %>%
    arrange(Set_order, desc(Avg_NSSE), Abbr) %>%
    rowwise() %>%
    mutate(
      bold_B = as.numeric(Baseline) < as.numeric(.data[[mixed_col]]),
      
      NSSE_B_str = sprintf("%.4f", as.numeric(Baseline)),
      NSSE_M_str = sprintf("%.4f", as.numeric(.data[[mixed_col]])),
      Avg_NSSE_str = sprintf("%.4f", as.numeric(Avg_NSSE)),
      
      NSSE_B_tex = ifelse(
        bold_B,
        paste0("\\textbf{", NSSE_B_str, "}"),
        NSSE_B_str
      ),
      NSSE_M_tex = ifelse(
        !bold_B,
        paste0("\\textbf{", NSSE_M_str, "}"),
        NSSE_M_str
      ),
      BF_tex = ifelse(
        is.infinite(as.numeric(BF_eq)),
        "$\\infty$",
        sprintf("%.2f", as.numeric(BF_eq))
      )
    ) %>%
    ungroup() %>%
    select(
      State = Abbr,
      Set,
      NSSE_B = NSSE_B_tex,
      NSSE_M = NSSE_M_tex,
      Avg_NSSE = Avg_NSSE_str,
      BF_eq = BF_tex,
      Evidence
    )
  
  n_included <- sum(df_table$Set == "Included")
  
  col_names <- c(
    "State",
    "Set",
    "Median $NSSE_B$",
    "Median $NSSE_M$",
    "Avg. NSSE",
    "$BF_{eq}$",
    "Evidence"
  )
  
  k_latex <- df_table %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      col.names = col_names,
      align = "llrrrrl",
      table.envir = NULL
    ) %>%
    row_spec(n_included, extra_latex_after = "\\midrule")
  
  k_latex_str <- as.character(k_latex)
  k_latex_str <- gsub("\\\\addlinespace\n?", "", k_latex_str)
  
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

# Supplement table: excluded locations and state-selection diagnostics
generate_excluded_locations_tex(
  file.path(res_dir, "table_excluded_locations.tex")
)

# NOTE: New Supplement sensitivity table mixed model comparison across all 51 locations
generate_sensitivity_mixed_tex(
  file.path(res_dir, "table_sensitivity_mixed.tex")
)
