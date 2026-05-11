#!/usr/bin/env Rscript
# R/plot_fits_data_bad_locations.R
#
# Produces: 
# - figures/extra/bad_locations_fits_mixed.pdf
# - figures/extra/bad_locations_fits_exp.pdf
# - figures/extra/bad_locations_fits_rational.pdf
#
# Each as a single figure for 32 excluded locations:
#   - 6 x 6 grid
#   - Baseline vs Behavioral
#   - 90% predictive intervals + medians + observed data points
#

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reticulate)
  library(tikzDevice)
})

project_root <- getwd() # assuming we are running it from the root

# load venv and config.py
use_virtualenv(file.path(project_root, ".venv"), required = TRUE)
py_run_string(sprintf("import sys; sys.path.append('%s')", file.path(project_root, "src")))
config <- import("config")

# data path
fits_path <- file.path(project_root, "data", "plotting", "fits_data.csv")

#
# ----- for each behavior variant ---
#
create_fits_grid_bad_locations <- function(behavior_model_name, output_filename) {
    
  # paths and order
  excluded_locs <- unlist(config$LOCATIONS_BAD)
  loc_codes <- sort(excluded_locs) # NOTE: in alphabetical order
  
  stopifnot(length(loc_codes) == length(excluded_locs))
  
  out_dir  <- file.path(project_root, "figures", "extra")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tex_file <- file.path(out_dir, paste0(output_filename, ".tex"))
  pdf_file <- file.path(out_dir, paste0(output_filename, ".pdf"))
  
  # models and colors
  model_order_2 <- c("Baseline", behavior_model_name)
  
  py_models <- config$MODELS
  disp <- vapply(py_models, function(m) m$display_name, character(1))
  cols <- vapply(py_models, function(m) m$color, character(1))
  names(cols) <- disp
  
  model_colors_2 <- cols[model_order_2]
  
  # theme config
  font_eps <- 0
  FONT_SIZES <- list(
    base        = 7 + font_eps,
    axis_title  = 8 + font_eps,
    strip_text  = 7 + font_eps,
    legend_text = 7 + font_eps
  )
  
  theme_publication <- theme_bw(base_size = FONT_SIZES$base) +
    theme(
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      axis.title         = element_text(size = FONT_SIZES$axis_title),
      strip.background   = element_rect(fill = "gray90", colour = "black", linewidth = 0.4),
      strip.text         = element_text(size = FONT_SIZES$strip_text, face = "bold"),
      strip.placement    = "outside",
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(size = FONT_SIZES$legend_text),
      legend.key.width   = unit(1.0, "cm"),
      legend.key.height  = unit(0.6, "cm"),
      plot.margin        = margin(6, 6, 6, 6)
    )
  
  # load data and normalize model names
  fits_df <- read.csv(fits_path) %>%
    mutate(date = as.Date(date))
  
  fits_df$model_name <- gsub(
    "Behavioral \\(Exp\\)",
    "Behavioral (Exponential)",
    fits_df$model_name
  )
  
  fits_df <- fits_df %>%
    filter(
      location %in% loc_codes,
      model_name %in% model_order_2
    )
  
  fits_df$location <- factor(fits_df$location, levels = loc_codes)
  fits_df$model_name <- factor(fits_df$model_name, levels = model_order_2)
  
  if (length(unique(fits_df$location)) != length(loc_codes)) {
    missing_locs <- setdiff(loc_codes, unique(as.character(fits_df$location)))
    stop(
      paste(
        "Mismatch between excluded locations and locations in fits_df. Missing:",
        paste(missing_locs, collapse = ", ")
      )
    )
  }
  
  # legend entries
  legend_breaks <- c(model_order_2, "Observed data")
  color_values <- c(model_colors_2, "Observed data" = "black")
  fill_values  <- c(model_colors_2, "Observed data" = NA)
  linetype_values <- c(
    setNames(rep("solid", length(model_order_2)), model_order_2),
    "Observed data" = "blank"
  )
  
  # plot
  p_fits <- ggplot() +
    geom_ribbon(
      data = fits_df,
      aes(x = date, ymin = lower_90, ymax = upper_90, fill = model_name),
      alpha = 0.40,
      show.legend = TRUE
    ) +
    geom_line(
      data = fits_df,
      aes(x = date, y = median, color = model_name, linetype = model_name),
      linewidth = 0.55,
      show.legend = TRUE
    ) +
    geom_point(
      data = fits_df,
      aes(x = date, y = obs, color = "Observed data"),
      size = 0.25,
      shape = 19,
      show.legend = TRUE
    ) +
    facet_wrap(~ location, ncol = 6, scales = "free_y") +
    scale_color_manual(
      name = NULL,
      values = color_values,
      breaks = legend_breaks,
      drop = FALSE
    ) +
    scale_fill_manual(
      name = NULL,
      values = fill_values,
      breaks = legend_breaks,
      drop = FALSE
    ) +
    scale_linetype_manual(
      name = NULL,
      values = linetype_values,
      breaks = legend_breaks,
      drop = FALSE
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        override.aes = list(
          fill      = c(
            adjustcolor(model_colors_2[1], alpha.f = 0.3),
            adjustcolor(model_colors_2[2], alpha.f = 0.3),
            NA
          ),
          linetype  = c("solid", "solid", "blank"),
          shape     = c(NA, NA, 19),
          linewidth = c(0.55, 0.55, 0)
        )
      ),
      fill = "none",
      linetype = "none"
    ) +
    labs(x = "Date", y = "Daily deaths (7-day average)") +
    theme_publication
  
  # export via tikz
  tikz(
    tex_file,
    width = 12.5,
    height = 10.5,
    standAlone = TRUE,
    engine = "luatex",
    packages = c(
      "\\usepackage{amsmath}",
      "\\usepackage{tikz}",
      "\\usepackage[active,tightpage,psfixbb]{preview}",
      "\\PreviewEnvironment{pgfpicture}",
      "\\setlength\\PreviewBorder{0pt}",
      "\\renewcommand{\\familydefault}{\\sfdefault}"
    )
  )
  print(p_fits)
  dev.off()
  
  # compile and clean up
  system(paste(
    "lualatex -halt-on-error -interaction=nonstopmode -output-directory",
    shQuote(out_dir),
    shQuote(tex_file)
  ))
  
  aux_base <- sub("\\.tex$", "", basename(tex_file))
  invisible(lapply(
    file.path(out_dir, paste0(aux_base, c(".aux", ".log"))),
    function(f) if (file.exists(f)) file.remove(f)
  ))
  
  cat("Saved:", pdf_file, "\n")
}

# Mixed
create_fits_grid_bad_locations("Behavioral (Mixed)", "bad_locations_fits_mixed")

# Exponential
create_fits_grid_bad_locations("Behavioral (Exponential)", "bad_locations_fits_exp")

# Rational
create_fits_grid_bad_locations("Behavioral (Rational)", "bad_locations_fits_rational")

