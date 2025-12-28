#!/usr/bin/env Rscript
# R/plot_synthetic_results.R
#
# Produces: 
# - figures/main/F3_synthetic_experiment_mixed.pdf
# - figures/main/S2_synthetic_experiment_exp.pdf
# - figures/main/S3_synthetic_experiment_rational.pdf
#
# Each is a 3-panel figure comparing baseline with behavior (variant) 
# results of synthetic experiments:
# - Panel A: fits
# - Panel B: R_0 boxplots
# - Panel C: final sizes boxplots
# 
# With ground-truth values for R_0 and final epidemic sizes.
#

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reticulate)
  library(patchwork)
  library(tikzDevice)
  library(ggh4x)  
})

project_root <- getwd() # assuming running from root directory

# load venv and config.py
use_virtualenv(file.path(project_root, ".venv"), required = TRUE)
py_run_string(sprintf("import sys; sys.path.append('%s')", file.path(project_root, "src")))
config <- import("config")

# data paths
plotting_dir <- file.path(project_root, "data", "plotting")
true_r0_path <- file.path(plotting_dir, "synthetic_true_r0.csv")
true_fs_path <- file.path(plotting_dir, "synthetic_true_final_sizes.csv")

create_synthetic_plot <- function(variant_code, out_filename, out_subdir) {
  
  # map display names mapping
  display_name <- switch(variant_code, 
    "mixed"    = "Behavioral (Mixed)", 
    "exp"      = "Behavioral (Exponential)", 
    "rational" = "Behavioral (Rational)"
  )
  
  # load plotting data
  fits_path <- file.path(plotting_dir, paste0("fits_data_synthetic_ground_truth_", variant_code, ".csv"))
  r0_path   <- file.path(plotting_dir, paste0("r0_boxplot_data_synthetic_ground_truth_", variant_code, ".csv"))
  fs_path   <- file.path(plotting_dir, paste0("final_size_boxplot_data_synthetic_ground_truth_", variant_code, ".csv"))
  
  # load ground-truth values 
  true_r0_df <- read.csv(true_r0_path)
  TRUE_R0    <- true_r0_df[[paste0("true_r0_", variant_code)]][1]
  
  true_fs_all <- read.csv(true_fs_path)
  true_final_size_df <- true_fs_all %>%
    dplyr::select(zeta_true, true_final_size = paste0("true_final_size_", variant_code))

  # settings
  model_order_2 <- c("Baseline", display_name)
  py_models     <- config$MODELS
  cols_map      <- setNames(vapply(py_models, function(m) m$color, character(1)), 
                            vapply(py_models, function(m) m$display_name, character(1)))
  model_colors_2 <- cols_map[model_order_2]
  MODEL_SHAPES_2 <- setNames(c(16, 17), model_order_2)

  font_eps <- 2
  FONT_SIZES <- list(base=9+font_eps, axis_title=10+font_eps, strip_text=10+font_eps, legend_text=9+font_eps)
  theme_publication <- theme_bw(base_size = FONT_SIZES$base) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          axis.title = element_text(size = FONT_SIZES$axis_title),
          strip.background = element_rect(fill = "gray90", colour = "black", linewidth = 0.5),
          strip.text = element_text(size = FONT_SIZES$strip_text, face = "bold"),
          strip.placement = "outside", legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = FONT_SIZES$legend_text),
          legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.7, "cm"),
          plot.margin = margin(6, 6, 6, 6))

  ###
  ## --- Panel A: fits ---
  #
  fits_df <- read.csv(fits_path) %>% mutate(date = as.Date(date))

  # map names to display_names
  fits_df$model_name <- gsub("\\(Mix\\)", "(Mixed)", fits_df$model_name)
  fits_df$model_name <- gsub("\\(Exp\\)", "(Exponential)", fits_df$model_name)
  fits_df$model_name <- gsub("\\(Rat\\)", "(Rational)", fits_df$model_name)
  
  fits_df <- fits_df %>% filter(model_name %in% model_order_2)
  z_order <- sort(unique(fits_df$zeta_true))
  fits_df$facet_zeta <- factor(fits_df$zeta_true, levels = z_order, labels = sprintf("$\\zeta = %.3f$", z_order))
  fits_df$model_name <- factor(fits_df$model_name, levels = model_order_2)

  legend_breaks_A <- c(model_order_2, "Synthetic data")
  color_values_A  <- c(model_colors_2, "Synthetic data" = "black")
  fill_values_A     <- c(model_colors_2, "Synthetic data" = NA)
  linetype_values_A <- c(setNames(rep("solid", length(model_order_2)), model_order_2), "Synthetic data" = "blank")

  p_fits <- ggplot() +
    geom_ribbon(data = fits_df, aes(x = date, ymin = lower_90, 
                                    ymax = upper_90, fill = model_name), 
                                    alpha = 0.40, show.legend = TRUE) +
    geom_line(data = fits_df, aes(x = date, y = median, color = model_name, 
      linetype = model_name), linewidth = 0.8, show.legend = TRUE) +
    geom_point(data = fits_df, aes(x = date, y = obs, color = "Synthetic data"), 
      size = 0.3, shape = 16, show.legend = TRUE) +
    facet_wrap(~ facet_zeta, scales = "free_y", ncol = 3) +
    scale_color_manual(name = NULL, values = color_values_A, 
      breaks = legend_breaks_A, drop = FALSE) +
    scale_fill_manual(name = NULL, values = fill_values_A, 
      breaks = legend_breaks_A, drop = FALSE) +
    scale_linetype_manual(name = NULL, values = linetype_values_A, 
      breaks = legend_breaks_A, drop = FALSE) +
    guides(color = guide_legend(nrow = 1, override.aes = list(
        fill = c(adjustcolor(model_colors_2[1], alpha.f = 0.3), 
                 adjustcolor(model_colors_2[2], alpha.f = 0.3), NA),
        linetype  = c("solid", "solid", "blank"), 
        shape = c(NA, NA, 16), linewidth = c(0.8, 0.8, 0))),
      fill = "none", linetype = "none") +
    labs(x = "Date", y = "Daily deaths (7-day average)") + theme_publication

  ###
  ## --- Panel B: R_0 boxplots ---
  #
  r0_df <- read.csv(r0_path)  
  # map display names
  r0_df$model_name <- gsub("\\(Mix\\)", "(Mixed)", r0_df$model_name)
  r0_df$model_name <- gsub("\\(Exp\\)", "(Exponential)", r0_df$model_name)
  r0_df$model_name <- gsub("\\(Rat\\)", "(Rational)", r0_df$model_name)
  r0_df <- r0_df %>% filter(model_name %in% model_order_2)
  r0_df$facet_zeta <- factor(r0_df$zeta_true, levels = z_order, 
    labels = sprintf("$\\zeta = %.3f$", z_order))
  r0_df$model_name <- factor(r0_df$model_name, levels = model_order_2)

  p_r0 <- ggplot(r0_df, aes(x = model_name)) +
    geom_boxplot(aes(ymin=ymin, lower=lower, 
                     middle=middle, upper=upper, 
                     ymax=ymax, fill=model_name, 
                     color=model_name),
                 stat="identity", alpha=0.5, width=0.7, linewidth=0.7) +
    geom_point(aes(y=middle, shape=model_name, color=model_name), size=2.6, stroke=0.25) +
    geom_hline(aes(yintercept=TRUE_R0, color="True basic reproduction number", 
      linetype="True basic reproduction number"), linewidth=0.7, inherit.aes=FALSE) +
    facet_grid(. ~ facet_zeta, scales="free_x", space="free_x") +
    scale_color_manual(name=NULL, values=c("Baseline"="black", 
      setNames("black", display_name), "True basic reproduction number"="red"),
                       breaks=c(model_order_2, "True basic reproduction number"), drop=FALSE) +
    scale_fill_manual(name=NULL, values=model_colors_2, breaks=model_order_2, drop=FALSE) +
    scale_shape_manual(name=NULL, values=MODEL_SHAPES_2, breaks=model_order_2) +
    scale_linetype_manual(name=NULL, values=c("Baseline"="solid", 
      "Behavioral"="solid", "True basic reproduction number"="dashed"), drop=FALSE) +
    guides(fill="none", shape="none", linetype="none", 
      color=guide_legend(keywidth=unit(1.5, "cm"), keyheight=unit(0.8, "cm"),
      override.aes = list(linetype=c("solid","solid","dashed"),
                          fill=c(adjustcolor(model_colors_2[1], 0.5), adjustcolor(model_colors_2[2], 0.5), NA),
                          shape=c(16,17,NA)))) +
    labs(x=NULL, y="Basic reproduction number") + theme_publication + 
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

  ###
  ## --- Panel C: final sizes boxplots ---
  #
  fs_df <- read.csv(fs_path)
  fs_df$model_name <- gsub("\\(Mix\\)", "(Mixed)", fs_df$model_name)
  fs_df$model_name <- gsub("\\(Exp\\)", "(Exponential)", fs_df$model_name)
  fs_df$model_name <- gsub("\\(Rat\\)", "(Rational)", fs_df$model_name)
  fs_df <- fs_df %>% filter(model_name %in% model_order_2)
  fs_df$facet_zeta <- factor(fs_df$zeta_true, levels=z_order, labels=sprintf("$\\zeta = %.3f$", z_order))
  fs_df$model_name <- factor(fs_df$model_name, levels=model_order_2)
  true_fs_plot_df <- true_final_size_df %>% mutate(facet_zeta = factor(zeta_true, levels=z_order, labels=sprintf("$\\zeta = %.3f$", z_order)))

  p_fs <- ggplot(fs_df, aes(x = model_name)) +
    geom_boxplot(aes(ymin=ymin, lower=lower, middle=middle, upper=upper, ymax=ymax, fill=model_name, color=model_name),
                 stat="identity", alpha=0.5, width=0.7, linewidth=0.7) +
    geom_point(aes(y=middle, shape=model_name, color=model_name), size=2.6, stroke=0.25) +
    geom_hline(data=true_fs_plot_df, aes(yintercept=true_final_size, color="True final epidemic size", linetype="True final epidemic size"),
               linewidth=0.7, inherit.aes=FALSE) +
    facet_grid(. ~ facet_zeta, scales="free_x", space="free_x") +
    scale_fill_manual(name=NULL, values=model_colors_2, breaks=model_order_2, drop=FALSE) +
    scale_shape_manual(name=NULL, values=MODEL_SHAPES_2, breaks=model_order_2) +
    scale_color_manual(name=NULL, values=c("Baseline"="black", setNames("black", display_name), "True final epidemic size"="red"),
                       breaks=c(model_order_2, "True final epidemic size"), drop=FALSE) +
    scale_linetype_manual(name=NULL, values=c("Baseline"="solid", "Behavioral"="solid", "True final epidemic size"="dashed"), drop=FALSE) +
    guides(fill="none", shape="none", linetype="none", color=guide_legend(keywidth=unit(1.5, "cm"), keyheight=unit(0.8, "cm"),
      override.aes = list(linetype=c("solid","solid","dashed"),
                          fill=c(adjustcolor(model_colors_2[1], 0.5), adjustcolor(model_colors_2[2], 0.5), NA),
                          shape=c(16,17,NA)))) +
    labs(x=NULL, y="Final epidemic size") + theme_publication + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

  ###
  ## -- combine panels and export ---
  #
  combined <- (p_fits / p_r0 / p_fs) + plot_layout(heights = c(2, 1, 1.2)) + 
              plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16, face = "bold"))

  out_dir_full <- file.path(project_root, "figures", out_subdir)
  dir.create(out_dir_full, recursive = TRUE, showWarnings = FALSE)
  tex_file <- file.path(out_dir_full, paste0(out_filename, ".tex"))
  
  tikz(tex_file, width = 11.5, height = 12.5, standAlone = TRUE,
       packages = c("\\usepackage{amsmath}", "\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
                    "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\renewcommand{\\familydefault}{\\sfdefault}"))
  print(combined)
  dev.off()

  system(paste("pdflatex -halt-on-error -interaction=nonstopmode -output-directory", shQuote(out_dir_full), shQuote(tex_file)))
  invisible(lapply(file.path(out_dir_full, paste0(out_filename, c(".aux", ".log"))), function(f) if (file.exists(f)) file.remove(f)))
  cat("Saved:", file.path(out_dir_full, paste0(out_filename, ".pdf")), "\n")
}

# Main Figure 3 (Mixed)
create_synthetic_plot("mixed", "F3_synthetic_experiment_mixed", "main")

# Supplement Figure S2 (Exponential)
create_synthetic_plot("exp", "S2_synthetic_experiment_exp", "supplement")

# Supplement Figure S3 (Rational)
create_synthetic_plot("rational", "S3_synthetic_experiment_rational", "supplement")

