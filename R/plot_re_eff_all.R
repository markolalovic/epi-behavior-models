#!/usr/bin/env Rscript
# R/plot_re_eff_all.R
#
# Produces:
#   - figures/supplement/S12_r_effective_mixed.pdf
#   - figures/supplement/S13_r_effective_exp.pdf
#   - figures/supplement/S14_r_effective_rational.pdf
#
# Each figure:
#   - compares baseline and one behavior variant
#   - 5 x 6 grid (for 30 US states)
#   - R_e(t) medians and 90% pointwise credible intervals
#

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reticulate)
  library(tikzDevice)
})

###
## --- settings ---
#
project_root <- getwd() # assuming running from root

# load venv and config.py
use_virtualenv(file.path(project_root, ".venv"), required = TRUE)
py_run_string(sprintf("import sys; sys.path.append('%s')", file.path(project_root, "src")))
config <- import("config")

# data paths
re_path <- file.path(project_root, "data", "plotting", "re_effective_data.csv")

out_dir <- file.path(project_root, "figures", "supplement")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

font_eps <- 1 
FONT_SIZES <- list(base=8+font_eps, axis_title=9+font_eps, 
  strip_text=8+font_eps, legend_text=8+font_eps)

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

# load data 
re_df <- read.csv(re_path, stringsAsFactors = FALSE) %>%
  mutate(date = as.Date(date))

# NOTE: keeping alphabetical ordering as in config.LOCATIONS
loc_codes <- unlist(config$LOCATIONS)
re_df$location <- factor(re_df$location, levels = loc_codes)

###
## --- plotting ---
#
plot_re_variant <- function(behavioral_display_name, filename_prefix) {
  
  model_order_2 <- c("Baseline", behavioral_display_name)
  
  py_models <- config$MODELS
  cols_map  <- setNames(vapply(py_models, function(m) m$color, character(1)), 
                        vapply(py_models, function(m) m$display_name, character(1)))
  model_colors_2 <- cols_map[model_order_2]

  p_re <- ggplot(re_df %>% filter(model_name %in% model_order_2)) +
    # dashed reference line at R_e = 1
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.4) +
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
    # shared y-axis
    scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2)) +
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
    labs(x = "Date", y = "Effective reproduction number") +
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
  print(p_re)
  dev.off()

  system(paste("pdflatex -halt-on-error -interaction=nonstopmode -output-directory", 
               shQuote(out_dir), shQuote(tex_path)))
  
  aux_base <- sub("\\.tex$", "", basename(tex_path))
  invisible(lapply(file.path(out_dir, paste0(aux_base, c(".aux", ".log"))), 
                   function(f) if (file.exists(f)) file.remove(f)))

  cat("Saved:", pdf_path, "\n")
}


# S12: Mixed
plot_re_variant("Behavioral (Mixed)", "S12_r_effective_mixed")

# S13: Exponential
plot_re_variant("Behavioral (Exponential)", "S13_r_effective_exp")

# S14: Rational
plot_re_variant("Behavioral (Rational)", "S14_r_effective_rational")

