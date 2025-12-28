#!/usr/bin/env Rscript
# R/plot_bad_locations.R
#
# Produces single figure:
#   - Grid of 21 excluded states (7 x 3)
#   - Observed 7-day averaged daily deaths only (no fits)
#   - Saves it to `figures/supplement/bad_locations.pdf`
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
py_run_string("import sys; import os; sys.path.append(os.path.join(os.getcwd(), 'src'))")
config <- import("config")

# data paths
smoothed_path <- file.path(project_root, "data", "processed", "smoothed_mortality.csv")
out_dir       <- file.path(project_root, "figures", "supplement")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tex_file <- file.path(out_dir, "S1_bad_locations.tex")
pdf_file <- file.path(out_dir, "S1_bad_locations.pdf")

# 21 excluded states from config.LOCATIONS_BAD
bad_locs <- unlist(config$LOCATIONS_BAD)


# theme
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
    legend.position    = "none",
    plot.margin        = margin(6, 6, 6, 6)
  )


# load smoothed mortality data 
smoothed <- read.csv(smoothed_path, stringsAsFactors = FALSE)
smoothed$date <- as.Date(smoothed$date)

# keep only date and bad location columns
bad_long <- smoothed %>%
  select(date, all_of(bad_locs)) %>%
  pivot_longer(
    cols      = -date,
    names_to  = "state",
    values_to = "deaths_7d"
  )

# restrict to the main analysis window
start_date <- as.Date("2020-03-01")
end_date   <- as.Date("2020-07-01")

bad_long <- bad_long %>%
  filter(date >= start_date, date <= end_date)

# order alphabetically
bad_long$state <- factor(bad_long$state, levels = sort(bad_locs))


# plot
p_bad <- ggplot(bad_long, aes(x = date, y = deaths_7d)) +
  geom_point(size = 0.6) +
  facet_wrap(~ state, ncol = 6, scales = "free_y") +
  labs(
    x = "Date",
    y = "Daily deaths (7-day average)"
  ) +
  theme_publication


# export via tikz and pdflatex
tikz(
  tex_file,
  width  = 11.5,
  height = 8.5,
  standAlone = TRUE,
  packages = c(
    "\\usepackage{amsmath}",
    "\\usepackage{tikz}",
    "\\usepackage[active,tightpage,psfixbb]{preview}",
    "\\PreviewEnvironment{pgfpicture}",
    "\\setlength\\PreviewBorder{0pt}",
    "\\renewcommand{\\familydefault}{\\sfdefault}"
  )
)
print(p_bad)
dev.off()

cmd <- paste(
  "pdflatex -halt-on-error -interaction=nonstopmode",
  "-output-directory", shQuote(out_dir),
  shQuote(tex_file)
)
system(cmd)

# clean up
aux_base  <- sub("\\.tex$", "", basename(tex_file))
aux_files <- file.path(out_dir, paste0(aux_base, c(".aux", ".log")))
invisible(lapply(aux_files, function(f) if (file.exists(f)) file.remove(f)))

cat("Saved:", pdf_file, "\n")
