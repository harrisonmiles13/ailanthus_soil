suppressMessages(library(ggplot2))

FIG_V2_DIR <- "output/figures_v2"
dir.create(FIG_V2_DIR, showWarnings = FALSE, recursive = TRUE)
fig_v2 <- function(name) file.path(FIG_V2_DIR, name)

PLOT_RECODE <- c(no_fungus_control = "No-fungus control",
                 attenuated        = "Formerly attenuated strain",
                 virulent          = "Virulent strain")
PLOT_ORDER  <- c("No-fungus control", "Formerly attenuated strain", "Virulent strain")
recode_plot <- function(x)
  factor(unname(PLOT_RECODE[as.character(x)]), levels = PLOT_ORDER)

PLOT_COLOURS <- c("No-fungus control"          = "#0072B2",
                  "Formerly attenuated strain" = "#D55E00",
                  "Virulent strain"            = "#009E73")
PLOT_LINETYPES <- c("No-fungus control"          = "solid",
                    "Formerly attenuated strain" = "21",
                    "Virulent strain"            = "4121")

PLOT_LEGEND <- "Plot (within-plot nuisance)"
scale_plot_colour   <- function(name = PLOT_LEGEND)
  scale_colour_manual(values = PLOT_COLOURS, name = name)
scale_plot_linetype <- function(name = PLOT_LEGEND)
  scale_linetype_manual(values = PLOT_LINETYPES, name = name)

MONTH_SHAPES <- c(may = 17, july = 16, september = 15)
scale_month_shape <- function(name = "Sampling month", labels = waiver())
  scale_shape_manual(values = MONTH_SHAPES, name = name, labels = labels)

LAB_ANI <- expression("Disease-influence index ("*ANI[disease]*")")

theme_v2 <- function(base_size = 16, legend = "bottom") {
  theme_bw(base_size = base_size) +
    theme(
      plot.title        = element_blank(),
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
      legend.box         = "vertical",
      legend.box.just    = "center",
      legend.title       = element_text(size = rel(0.88), face = "bold"),
      legend.text        = element_text(size = rel(0.82)),
      legend.key         = element_blank(),
      legend.background  = element_blank(),
      plot.margin        = margin(10, 16, 8, 10)
    )
}

annot_layer <- function(annot, size = 4.0) {
  geom_label(data = annot, aes(label = label), x = Inf, y = Inf,
             hjust = 1.04, vjust = 1.12, size = size, colour = "grey10",
             fill = "white", alpha = 0.78, label.size = 0.2,
             label.padding = unit(0.28, "lines"), lineheight = 0.95,
             inherit.aes = FALSE)
}
