#!/usr/bin/env Rscript
# R/plot_r0_and_final_size_boxplots.R
#
# Produces:
# - figures/main/F2_r0_final_size_mixed.pdf
# - figures/supplement/S6_r0_final_size_exp.pdf
# - figures/supplement/S7_r0_final_size_rational.pdf
#
# Each is a two-panel figure comparing 
# baseline with behavior (variant) across 30 locations:
#
# - Panel A: R_0 boxplots
# - Panel B: Final epidemic size boxplots
#

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reticulate)
  library(patchwork)
  library(tikzDevice)
})

project_root <- getwd() # assuming running from root directory

# load venv and config.py
use_virtualenv(file.path(project_root, ".venv"), required = TRUE)
py_run_string(sprintf("import sys; import os; sys.path.append(os.path.join(os.getcwd(), 'src'))"))
config <- import("config")

# data paths
r0_path    <- file.path(project_root, "data", "plotting", "r0_boxplot_data.csv")
final_path <- file.path(project_root, "data", "plotting", "final_size_boxplot_data.csv")

# theme setup
font_eps <- 1
FONT_SIZES <- list(
    base=8+font_eps, 
    axis_title=9+font_eps, 
    strip_text=9+font_eps, 
    legend_text=8+font_eps)

theme_publication <- theme_bw(base_size = FONT_SIZES$base) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = FONT_SIZES$axis_title),
    strip.background = element_rect(fill = "gray90", colour = "black", linewidth = 0.4),
    strip.text = element_text(size = FONT_SIZES$strip_text, face = "bold"),
    legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(size = FONT_SIZES$legend_text),
    legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.6, "cm"),
    plot.margin = margin(6, 6, 6, 6)
  )

# boxplot panel builder
create_boxplot_figure <- function(behavior_model_name, behavior_color, 
  output_filename, order_filename, out_subdir) {
  
  order_path <- file.path(project_root, "results", order_filename)
  out_dir    <- file.path(project_root, "figures", out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  tex_file <- file.path(out_dir, paste0(output_filename, ".tex"))
  pdf_file <- file.path(out_dir, paste0(output_filename, ".pdf"))

  model_order_2 <- c("Baseline", behavior_model_name)
  model_fills_2 <- c("Baseline" = "#0072B2", setNames(behavior_color, behavior_model_name))
  model_shapes_2 <- c("Baseline" = 16, setNames(17, behavior_model_name))

  # load order and data
  loc_order <- read.csv(order_path)$Location
  r0_df     <- read.csv(r0_path)
  final_df  <- read.csv(final_path)

  # naming mapping
  target_suffix <- gsub("Behavioral ", "", behavior_model_name)
  r0_df$model_name    <- gsub("\\(Exp\\)", paste0("(", target_suffix, ")"), r0_df$model_name)
  final_df$model_name <- gsub("\\(Exp\\)", paste0("(", target_suffix, ")"), final_df$model_name)

  # prepare dataframe
  prep_df <- function(df) {
    df %>% 
      filter(model_name %in% model_order_2) %>%
      mutate(model_name = factor(model_name, levels = model_order_2),
             location = factor(location, levels = loc_order))
  }
  r0_df <- prep_df(r0_df); final_df <- prep_df(final_df)

  # vertical alternating stripes
  even_stripes <- data.frame(idx = seq_along(loc_order)) %>% filter(idx %% 2 == 0)

  # make panel function
  make_panel <- function(df, y_lab, show_x = TRUE) {
    dodge <- position_dodge(width = 0.6)
    p <- ggplot() +
      geom_rect(data=even_stripes, aes(xmin=idx-0.5, xmax=idx+0.5, ymin=-Inf, ymax=Inf), fill="grey95", colour=NA) +
      geom_boxplot(data=df, aes(x=location, ymin=ymin, lower=lower, middle=middle, upper=upper, ymax=ymax, fill=model_name),
                   stat="identity", width=0.5, alpha=0.6, linewidth=0.5, position=dodge, colour="black") +
      geom_point(data=df, aes(x=location, y=middle, shape=model_name), position=dodge, size=1.8, stroke=0.25, colour="black") +
      scale_fill_manual(values=model_fills_2) + scale_shape_manual(values=model_shapes_2) +
      labs(x=NULL, y=y_lab) + theme_publication +
      theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), 
            panel.grid.major.y = element_line(linewidth=0.3, colour="grey85"))
    if(!show_x) p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
    return(p)
  }

  p_a <- make_panel(r0_df, "Basic reproduction number", show_x = FALSE)
  p_b <- make_panel(final_df, "Final epidemic size", show_x = TRUE)

  # combine it
  combined <- (p_a / p_b) + plot_layout(heights=c(1,1), guides='collect') + plot_annotation(tag_levels='A') &
              theme(legend.position="bottom", plot.tag=element_text(size=14, face="bold"), panel.spacing=unit(0.1, "lines"))

  tikz(tex_file, width=11, height=7.5, standAlone=TRUE,
       packages = c("\\usepackage{amsmath}", "\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
                    "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\renewcommand{\\familydefault}{\\sfdefault}"))
  print(combined)
  dev.off()

  system(paste("pdflatex -interaction=nonstopmode -output-directory", shQuote(out_dir), shQuote(tex_file)))
  invisible(lapply(file.path(out_dir, paste0(output_filename, c(".aux", ".log"))), function(f) if (file.exists(f)) file.remove(f)))
  cat("Saved:", pdf_file, "\n")
}

# Main Figure 2
create_boxplot_figure("Behavioral (Mixed)", "#E69F00", "F2_r0_final_size_mixed", "order_2_mixed.csv", "main")

# Supplement Figure S6
create_boxplot_figure("Behavioral (Exponential)", "#D55E00", "S6_r0_final_size_exp", "order_2_exp.csv", "supplement")

# Supplement Figure S7
create_boxplot_figure("Behavioral (Rational)", "#CC79A7", "S7_r0_final_size_rational", "order_2_rational.csv", "supplement")
