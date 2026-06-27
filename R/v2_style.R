# v2_style.R
# ---------------------------------------------------------------------------
# Shared styling for the v2 figure outputs. Design goals (v2 brief):
#   * LARGE TEXT             -- readable at print/screen size without zooming
#   * BETTER USE OF SPACE    -- generous panels, bottom legends, trimmed margins,
#                               short axis titles (time-matching detail lives in
#                               the manuscript caption, not on the figure)
#   * NO MAIN TITLES         -- plots carry no title/subtitle; identification is
#                               the caption's job
#   * ACCESSIBILITY          -- Okabe-Ito colour-blind-safe palette for the three
#                               plots, with REDUNDANT linetype on the fitted lines
#                               and shape on the points, so colour is never the
#                               only channel; high-contrast annotation labels.
#
# v1 stays untouched: the v1 scripts and output/figures/*.png are conserved; the
# v2 scripts source this file and write to output/figures_v2/.
# ---------------------------------------------------------------------------

suppressMessages(library(ggplot2))

# --- output location (kept separate from v1 output/figures) ------------------
FIG_V2_DIR <- "output/figures_v2"
dir.create(FIG_V2_DIR, showWarnings = FALSE, recursive = TRUE)
fig_v2 <- function(name) file.path(FIG_V2_DIR, name)

# --- readable plot/treatment labels -----------------------------------------
# In the data the plot factor is raw ("attenuated","no_fungus_control","virulent").
# These remain a within-plot nuisance block (never compared), but readable labels
# and a fixed order help the reader.
PLOT_RECODE <- c(no_fungus_control = "No-fungus control",
                 attenuated        = "Formerly attenuated strain",
                 virulent          = "Virulent strain")
PLOT_ORDER  <- c("No-fungus control", "Formerly attenuated strain", "Virulent strain")
recode_plot <- function(x)
  factor(unname(PLOT_RECODE[as.character(x)]), levels = PLOT_ORDER)

# --- colour-blind-safe encoding (Okabe-Ito) + redundant linetype -------------
PLOT_COLOURS <- c("No-fungus control"          = "#0072B2",   # blue
                  "Formerly attenuated strain" = "#D55E00",    # vermillion
                  "Virulent strain"            = "#009E73")    # bluish green
PLOT_LINETYPES <- c("No-fungus control"          = "solid",
                    "Formerly attenuated strain" = "21",       # dashed
                    "Virulent strain"            = "4121")     # dot-dash

PLOT_LEGEND <- "Plot (within-plot nuisance)"
scale_plot_colour   <- function(name = PLOT_LEGEND)
  scale_colour_manual(values = PLOT_COLOURS, name = name)
scale_plot_linetype <- function(name = PLOT_LEGEND)
  scale_linetype_manual(values = PLOT_LINETYPES, name = name)

# --- sampling month -> shape (a second redundant, colour-blind-safe channel) --
MONTH_SHAPES <- c(may = 17, july = 16, september = 15)   # triangle / circle / square
scale_month_shape <- function(name = "Sampling month")
  scale_shape_manual(values = MONTH_SHAPES, name = name)

# --- short, typeset axis label for the index --------------------------------
LAB_ANI <- expression("Disease-influence index ("*ANI[disease]*")")

# --- the v2 theme ------------------------------------------------------------
# base_size is deliberately large; per-element rel() sizes keep the hierarchy.
theme_v2 <- function(base_size = 16, legend = "bottom") {
  theme_bw(base_size = base_size) +
    theme(
      plot.title        = element_blank(),                       # no main titles
      plot.subtitle     = element_blank(),
      axis.title        = element_text(size = rel(1.00), colour = "grey15"),
      axis.text         = element_text(size = rel(0.82), colour = "grey25"),
      strip.text        = element_text(size = rel(0.92), face = "bold",
                                       colour = "grey10",
                                       margin = margin(5, 5, 5, 5)),
      strip.background   = element_rect(fill = "grey92", colour = NA),
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = "grey90"),
      panel.spacing      = unit(1.1, "lines"),
      legend.position    = legend,
      legend.box         = "vertical",     # stack the (month) and (plot) legends
      legend.box.just    = "center",
      legend.title       = element_text(size = rel(0.88), face = "bold"),
      legend.text        = element_text(size = rel(0.82)),
      legend.key         = element_blank(),
      legend.background  = element_blank(),
      plot.margin        = margin(10, 16, 8, 10)
    )
}

# --- high-contrast corner annotation (slope / P) -----------------------------
# A white, lightly-bordered label so the text stays legible over any background;
# placed in the top-RIGHT whitespace of each panel (where the data are sparse).
annot_layer <- function(annot, size = 4.0) {
  geom_label(data = annot, aes(label = label), x = Inf, y = Inf,
             hjust = 1.04, vjust = 1.12, size = size, colour = "grey10",
             fill = "white", alpha = 0.78, label.size = 0.2,
             label.padding = unit(0.28, "lines"), lineheight = 0.95,
             inherit.aes = FALSE)
}
