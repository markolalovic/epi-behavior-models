#!/usr/bin/env Rscript
# R/plot_beta_eff_all.R
#
# Produces:
#   - figures/supplement/S9_beta_eff_mixed.pdf
#   - figures/supplement/S10_beta_eff_exp.pdf
#   - figures/supplement/S11_beta_eff_rational.pdf
#
# Each figure:
#   - compares baseline and one behavioral variant
#   - in 5 x 6 grid (30 US states)
#   - showing medians and 90% uncertainty intervals
#

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reticulate)
  library(tikzDevice)
})

project_root <- getwd() # assuming running from root

# load venv and config.py
use_virtualenv(file.path(project_root, ".venv"), required = TRUE)
py_run_string(sprintf("import sys; sys.path.append('%s')", file.path(project_root, "src")))
config <- import("config")

# data paths
beta_path <- file.path(project_root, "data", "plotting", "beta_eff_data.csv")
out_dir <- file.path(project_root, "figures", "supplement")

###
## --- settings ---
#
font_eps <- 1
FONT_SIZES <- list(base=8+font_eps, axis_title=9+font_eps, strip_text=8+font_eps, legend_text=8+font_eps)

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

beta_df <- read.csv(beta_path, stringsAsFactors = FALSE) %>%
  mutate(date = as.Date(date))

# NOTE: keeping alphabetical ordering from config.LOCATIONS
loc_codes <- unlist(config$LOCATIONS)
beta_df$location <- factor(beta_df$location, levels = loc_codes)


###
## --- plotting ---
#
plot_beta_variant <- function(behavioral_display_name, filename_prefix) {
  
  model_order_2 <- c("Baseline", behavioral_display_name)
  
  py_models <- config$MODELS
  cols_map  <- setNames(vapply(py_models, function(m) m$color, character(1)), 
                        vapply(py_models, function(m) m$display_name, character(1)))
  model_colors_2 <- cols_map[model_order_2]

  # build the plot
  p_beta <- ggplot(beta_df %>% filter(model_name %in% model_order_2)) +
    geom_ribbon(
      aes(x = date, ymin = lower, ymax = upper, fill = model_name),
      alpha = 0.40, show.legend = TRUE
    ) +
    geom_line(
      aes(x = date, y = median, color = model_name),
      linewidth = 0.6, show.legend = TRUE
    ) +
    facet_wrap(~ location, ncol = 6) +
    scale_color_manual(name = NULL, values = model_colors_2, breaks = model_order_2) +
    scale_fill_manual(name = NULL, values = model_colors_2, breaks = model_order_2) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    guides(
      color = guide_legend(
        nrow = 1,
        override.aes = list(
          fill      = c(adjustcolor(model_colors_2[1], alpha.f = 0.3),
                         adjustcolor(model_colors_2[2], alpha.f = 0.3)),
          linewidth = c(0.6, 0.6)
        )
      ),
      fill = "none"
    ) +
    labs(x = "Date", y = "Effective transmission rate") +
    theme_publication

  tex_path <- file.path(out_dir, paste0(filename_prefix, ".tex"))
  pdf_path <- file.path(out_dir, paste0(filename_prefix, ".pdf"))

  tikz(tex_path, width = 11.5, height = 8.5, standAlone = TRUE,
    packages = c(
      "\\usepackage{amsmath}", "\\usepackage{tikz}",
      "\\usepackage[active,tightpage,psfixbb]{preview}",
      "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
      "\\renewcommand{\\familydefault}{\\sfdefault}"
    )
  )
  print(p_beta)
  dev.off()

  system(paste("pdflatex -interaction=nonstopmode -output-directory", 
               shQuote(out_dir), shQuote(tex_path)))
  
  aux_base <- sub("\\.tex$", "", basename(tex_path))
  invisible(lapply(file.path(out_dir, paste0(aux_base, c(".aux", ".log"))), 
                   function(f) if (file.exists(f)) file.remove(f)))

  cat("Saved Figure:", pdf_path, "\n")
}

# S9: Mixed
plot_beta_variant("Behavioral (Mixed)", "S9_beta_eff_mixed")

# S10: Exponential
plot_beta_variant("Behavioral (Exponential)", "S10_beta_eff_exp")

# S11: Rational
plot_beta_variant("Behavioral (Rational)", "S11_beta_eff_rational")
