#!/usr/bin/env Rscript
# R/plot_fits_data.R

# Produces: 
# - figures/main/F1_final_fits_mixed.pdf
# - figures/supplement/S4_final_fits_exp.pdf
# - figures/supplement/S5_final_fits_rational.pdf
#
# Each as a single figure for 20 selected locations:
#   - 4 rows x 5 columns
#   - Baseline vs Behavioral (Mixed)
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
create_fits_grid <- function(behavior_model_name, output_filename, order_filename, out_subdir) {
  
  # paths and order
  selected_locs <- unlist(config$LOCATIONS)

  order_path <- file.path(project_root, "data", "plotting", order_filename)
  order_df  <- read.csv(order_path, stringsAsFactors = FALSE)

  loc_codes <- order_df$Location
  loc_codes <- loc_codes[loc_codes %in% selected_locs]

  missing_from_order <- setdiff(selected_locs, loc_codes)
  if (length(missing_from_order) > 0) {
    loc_codes <- c(loc_codes, missing_from_order)
  }

  stopifnot(length(loc_codes) == length(selected_locs))

  out_dir  <- file.path(project_root, "figures", out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tex_file <- file.path(out_dir, paste0(output_filename, ".tex"))
  pdf_file <- file.path(out_dir, paste0(output_filename, ".pdf"))

  # models and colors
  model_order_2 <- c("Baseline", behavior_model_name)
  
  py_models <- config$MODELS
  disp <- vapply(py_models, function(m) m$display_name, character(1))
  cols <- vapply(py_models, function(m) m$color,        character(1))
  names(cols) <- disp
  
  model_colors_2 <- cols[model_order_2]

  # theme config
  font_eps <- 1
  FONT_SIZES <- list(
    base       = 8 + font_eps,
    axis_title = 9 + font_eps,
    strip_text = 8 + font_eps,
    legend_text= 8 + font_eps
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
      legend.key.width   = unit(1.2, "cm"),
      legend.key.height  = unit(0.7, "cm"),
      plot.margin        = margin(6, 6, 6, 6)
    )

  # load data and normalize model names
  fits_df <- read.csv(fits_path) %>% mutate(date = as.Date(date))
  
  # # map display name, Exp -> Exponential
  # target_suffix <- gsub("Behavioral ", "", behavior_model_name)
  # fits_df$model_name <- gsub("\\(Exp\\)", paste0("(", target_suffix, ")"), fits_df$model_name)
  # normalize model names, if old plotting data used Exp
  fits_df$model_name <- gsub(
    "Behavioral \\(Exp\\)",
    "Behavioral (Exponential)",
    fits_df$model_name
  )

  # fits_df <- fits_df %>% filter(model_name %in% model_order_2)
  # fits_df$location <- factor(fits_df$location, levels = loc_codes)
  # new filtering
  fits_df <- fits_df %>%
    filter(
      location %in% loc_codes,
      model_name %in% model_order_2
    )

  fits_df$location <- factor(fits_df$location, levels = loc_codes)
  fits_df$model_name <- factor(fits_df$model_name, levels = model_order_2)  
  if (length(unique(fits_df$location)) != length(selected_locs)) {
    stop("Mismatch between selected locations and locations in fits_df.")
  }

  # legend entries
  legend_breaks <- c(model_order_2, "Observed data")
  color_values <- c(model_colors_2, "Observed data" = "black")
  fill_values  <- c(model_colors_2, "Observed data" = NA)
  linetype_values <- c(setNames(rep("solid", length(model_order_2)), model_order_2), "Observed data" = "blank")

  # plot it
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
      linewidth = 0.6,
      show.legend = TRUE
    ) +
    geom_point(
      data = fits_df,
      aes(x = date, y = obs, color = "Observed data"),
      size = 0.3,
      shape = 19,
      show.legend = TRUE
    ) +
    # facet_wrap(~ location, ncol = 5, scales = "free_y") +
    facet_wrap(~ location, ncol = 6, scales = "free_y") +
    scale_color_manual(name = NULL, values = color_values, breaks = legend_breaks, drop = FALSE) +
    scale_fill_manual(name = NULL, values = fill_values, breaks = legend_breaks, drop = FALSE) +
    scale_linetype_manual(name = NULL, values = linetype_values, breaks = legend_breaks, drop = FALSE) +
    guides(
      color = guide_legend(
        nrow = 1,
        override.aes = list(
          fill      = c(adjustcolor(model_colors_2[1], alpha.f = 0.3), adjustcolor(model_colors_2[2], alpha.f = 0.3), NA),
          linetype  = c("solid", "solid", "blank"),
          shape     = c(NA, NA, 19),
          linewidth = c(0.6, 0.6, 0)
        )
      ),
      fill    = "none",
      linetype= "none"
    ) +
    labs(x = "Date", y = "Daily deaths (7-day average)") +
    theme_publication

  # export via tikz
  # tikz(tex_file, width = 10.5, height = 7.2, standAlone = TRUE,
  # NEW: # 6 column x 9 rows
  # tikz(tex_file, width = 13.0, height = 13.0, standAlone = TRUE, 
  tikz(tex_file, width = 14.0, height = 14.5, standAlone = TRUE, engine = "luatex",
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
  # system(paste("pdflatex -halt-on-error -interaction=nonstopmode -output-directory", 
  #              shQuote(out_dir), shQuote(tex_file)))
  system(paste(
    "lualatex -halt-on-error -interaction=nonstopmode -output-directory",
    shQuote(out_dir),
    shQuote(tex_file)
  ))  
  
  aux_base  <- sub("\\.tex$", "", basename(tex_file))
  invisible(lapply(file.path(out_dir, paste0(aux_base, c(".aux", ".log"))), 
                   function(f) if (file.exists(f)) file.remove(f)))

  cat("Saved:", pdf_file, "\n")
}

# Figure 1: Mixed
create_fits_grid("Behavioral (Mixed)", "F1_final_fits_mixed", "order_mixed.csv", "main")

# Figure S4: Exponential
create_fits_grid("Behavioral (Exponential)", "S4_final_fits_exp", "order_exp.csv", "supplement")

# Figure S5: Rational
create_fits_grid("Behavioral (Rational)", "S5_final_fits_rational", "order_rational.csv", "supplement")
