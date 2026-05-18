#!/usr/bin/env Rscript
# R/generate_tables.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(reticulate)
  library(knitr)
  library(kableExtra)
})

project_root <- getwd()

use_virtualenv(file.path(project_root, ".venv"), required = TRUE)
py_run_string(sprintf("import sys; sys.path.append('%s')", file.path(project_root, "src")))
config <- import("config")
RESULT <- config$RESULT

nssr_path <- file.path(
  project_root, "results", "model_comparison",
  paste0("median_nssr_distances_", RESULT, ".csv")
)
df_nssr <- read.csv(nssr_path, check.names = FALSE)


generate_comparison_tex <- function(variant_display_name, variant_code, out_path) {

  selected_locs <- unlist(config$LOCATIONS)

  bf_filename <- paste0("baseline_vs_", variant_code, "_summary_", RESULT, ".csv")
  bf_path <- file.path(project_root, "results", "model_selection", bf_filename)
  df_bf <- read.csv(bf_path, check.names = FALSE)

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

  df_merged <- df_nssr_selected %>%
    inner_join(df_bf, by = "Location")

  if (nrow(df_merged) != length(selected_locs)) {
    stop("Merged comparison table does not contain exactly the selected locations.")
  }

  bf_col_name <- names(df_merged)[grep("^BF", names(df_merged))]

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
    mutate(
      bold_B = as.numeric(NSSE_B) < as.numeric(NSSE_V),

      NSSE_B_str = sprintf("%.4f", as.numeric(NSSE_B)),
      NSSE_B_tex = ifelse(bold_B, paste0("\\textbf{", NSSE_B_str, "}"), NSSE_B_str),

      NSSE_V_str = sprintf("%.4f", as.numeric(NSSE_V)),
      NSSE_V_tex = ifelse(!bold_B, paste0("\\textbf{", NSSE_V_str, "}"), NSSE_V_str),

      BF_tex = ifelse(
        is.infinite(as.numeric(BF_val)),
        "$\\infty$",
        sprintf("%.2f", as.numeric(BF_val))
      )
    ) %>%
    ungroup() %>%
    select(State, NSSE_B = NSSE_B_tex, NSSE_V = NSSE_V_tex, BF_eq = BF_tex, Evidence)

  subscript <- switch(variant_code, "mixed" = "M", "exp" = "E", "rational" = "R")

  col_names <- c(
    "State",
    "Median $NSSE_{B}$",
    paste0("Median $NSSE_{", subscript, "}$"),
    "$BF_{eq}$",
    "Evidence"
  )

  k_latex <- df_formatted %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      col.names = col_names,
      align = "lrrrl",
      table.envir = NULL
    )

  writeLines(as.character(k_latex), out_path)
  cat("Generated table:", out_path, "\n")
}


generate_excluded_locations_tex <- function(out_path) {

  state_path <- file.path(project_root, "results", "state_selection_summary.csv")
  df_state <- read.csv(state_path, check.names = FALSE)

  excluded_locs <- unlist(config$LOCATIONS_BAD)

  df_table <- df_state %>%
    filter(location %in% excluded_locs) %>%
    mutate(
      reason_order = case_when(
        reason == "insufficient mortality signal" ~ 1,
        reason == "incomplete mortality wave" ~ 2,
        reason == "reporting discontinuity" ~ 3,
        TRUE ~ 4
      ),
      reason_short = case_when(
        reason == "insufficient mortality signal" ~ "Low signal",
        reason == "incomplete mortality wave" ~ "Incomplete wave",
        reason == "reporting discontinuity" ~ "Reporting anomaly",
        TRUE ~ reason
      ),
      diagnostic_note = case_when(
        reason == "insufficient mortality signal" & peak <= 10 ~
          "peak $\\leq 10$",
        reason == "insufficient mortality signal" & total_deaths < 50 ~
          "total $<50$",
        reason == "incomplete mortality wave" & days_after_peak < 21 ~
          "peak too late",
        reason == "incomplete mortality wave" & tail_to_peak > 0.40 ~
          "tail/peak $>0.4$",
        reason == "reporting discontinuity" ~
          "reporting discontinuity",
        TRUE ~ reason
      )
    ) %>%
    arrange(reason_order, location) %>%
    transmute(
      Location = location,
      `Primary reason` = reason_short,
      `max_t_y` = sprintf("%.2f", as.numeric(peak)),
      `Tail/peak` = sprintf("%.2f", as.numeric(tail_to_peak)),
      `Diagnostic note` = diagnostic_note
    )

  col_names <- c(
    "Location",
    "Primary reason",
    "$\\max_t y(t)$",
    "Tail/peak",
    "Diagnostic note"
  )

  k_latex <- df_table %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      col.names = col_names,
      align = "llrrl",
      table.envir = NULL
    )

  k_latex_str <- as.character(k_latex)
  k_latex_str <- gsub("\\\\addlinespace\n?", "", k_latex_str)

  writeLines(k_latex_str, out_path)
  cat("Generated table:", out_path, "\n")
}


generate_sensitivity_table_tex <- function(
  variant_display_name,
  variant_code,
  out_path
) {

  all_locs <- unlist(config$LOCATIONS_ALL)
  selected_locs <- unlist(config$LOCATIONS)
  excluded_locs <- unlist(config$LOCATIONS_BAD)

  bf1_path <- file.path(
    project_root, "results", "model_selection",
    paste0("baseline_vs_", variant_code, "_summary_1.csv")
  )

  bf2_path <- file.path(
    project_root, "results", "model_selection",
    paste0("baseline_vs_", variant_code, "_summary_2.csv")
  )

  bf1 <- read.csv(bf1_path, check.names = FALSE)
  bf2 <- read.csv(bf2_path, check.names = FALSE)

  bf1_col <- names(bf1)[grep("^BF", names(bf1))]

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

  bf2_std <- bf2 %>%
    filter(behavior_variant == paste0("behavior_", variant_code)) %>%
    transmute(
      Abbr = location,
      BF_eq = as.numeric(BF_eq),
      Evidence = Interpretation
    ) %>%
    inner_join(
      df_nssr %>% select(Location, Abbr),
      by = "Abbr"
    )

  bf_all <- bind_rows(bf1_std, bf2_std) %>%
    distinct(Abbr, .keep_all = TRUE)

  r0_path <- file.path(project_root, "data", "plotting", "r0_boxplot_data.csv")
  final_path <- file.path(project_root, "data", "plotting", "final_size_boxplot_data.csv")

  r0_df <- read.csv(r0_path, check.names = FALSE)
  final_df <- read.csv(final_path, check.names = FALSE)

  r0_wide <- r0_df %>%
    filter(model_name %in% c("Baseline", variant_display_name)) %>%
    select(location, model_name, middle) %>%
    pivot_wider(
      names_from = model_name,
      values_from = middle
    ) %>%
    transmute(
      Abbr = location,
      Delta_R0 = as.numeric(.data[[variant_display_name]]) - as.numeric(Baseline)
    )

  final_wide <- final_df %>%
    filter(model_name %in% c("Baseline", variant_display_name)) %>%
    select(location, model_name, middle) %>%
    pivot_wider(
      names_from = model_name,
      values_from = middle
    ) %>%
    transmute(
      Abbr = location,
      Delta_A = as.numeric(Baseline) - as.numeric(.data[[variant_display_name]])
    )

  df_table <- df_nssr %>%
    filter(Abbr %in% all_locs) %>%
    select(
      Location,
      Abbr,
      Baseline,
      all_of(variant_display_name)
    ) %>%
    inner_join(
      bf_all %>% select(Abbr, BF_eq, Evidence),
      by = "Abbr"
    ) %>%
    inner_join(r0_wide, by = "Abbr") %>%
    inner_join(final_wide, by = "Abbr") %>%
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
      )
    )

  missing_locs <- setdiff(all_locs, df_table$Abbr)
  if (length(missing_locs) > 0) {
    stop(paste(
      "Sensitivity table is missing locations:",
      paste(missing_locs, collapse = ", ")
    ))
  }

  df_table <- df_table %>%
    arrange(Set_order, desc(as.numeric(BF_eq)), Abbr) %>%
    rowwise() %>%
    mutate(
      bold_B = as.numeric(Baseline) < as.numeric(.data[[variant_display_name]]),

      NSSE_B_str = sprintf("%.4f", as.numeric(Baseline)),
      NSSE_V_str = sprintf("%.4f", as.numeric(.data[[variant_display_name]])),

      NSSE_B_tex = ifelse(
        bold_B,
        paste0("\\textbf{", NSSE_B_str, "}"),
        NSSE_B_str
      ),
      NSSE_V_tex = ifelse(
        !bold_B,
        paste0("\\textbf{", NSSE_V_str, "}"),
        NSSE_V_str
      ),

      Delta_R0_tex = sprintf("%.2f", as.numeric(Delta_R0)),
      Delta_A_tex  = sprintf("%.3f", as.numeric(Delta_A)),

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
      NSSE_V = NSSE_V_tex,
      Delta_R0 = Delta_R0_tex,
      Delta_A = Delta_A_tex,
      BF_eq = BF_tex,
      Evidence
    )

  n_included <- sum(df_table$Set == "Included")

  subscript <- switch(
    variant_code,
    "mixed" = "M",
    "exp" = "E",
    "rational" = "R"
  )

  col_names <- c(
    "State",
    "Set",
    "Median $NSSE_B$",
    paste0("Median $NSSE_", subscript, "$"),
    paste0("$R_{0,", subscript, "}-R_{0,B}$"),
    paste0("$A_B-A_", subscript, "$"),
    "$BF_{eq}$",
    "Evidence"
  )

  k_latex <- df_table %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      col.names = col_names,
      align = "llrrrrrl",
      table.envir = NULL
    ) %>%
    row_spec(n_included, extra_latex_after = "\\midrule")

  k_latex_str <- as.character(k_latex)
  k_latex_str <- gsub("\\\\addlinespace\n?", "", k_latex_str)

  writeLines(k_latex_str, out_path)
  cat("Generated table:", out_path, "\n")
}


res_dir <- file.path(project_root, "results")

generate_comparison_tex(
  "Behavioral (Mixed)", "mixed",
  file.path(res_dir, "table_comparison_baseline_mixed.tex")
)

generate_comparison_tex(
  "Behavioral (Exponential)", "exp",
  file.path(res_dir, "table_comparison_baseline_exp.tex")
)

generate_comparison_tex(
  "Behavioral (Rational)", "rational",
  file.path(res_dir, "table_comparison_baseline_rational.tex")
)

generate_excluded_locations_tex(
  file.path(res_dir, "table_excluded_locations.tex")
)

generate_sensitivity_table_tex(
  "Behavioral (Mixed)", "mixed",
  file.path(res_dir, "table_sensitivity_mixed.tex")
)

generate_sensitivity_table_tex(
  "Behavioral (Exponential)", "exp",
  file.path(res_dir, "table_sensitivity_exp.tex")
)

generate_sensitivity_table_tex(
  "Behavioral (Rational)", "rational",
  file.path(res_dir, "table_sensitivity_rational.tex")
)

