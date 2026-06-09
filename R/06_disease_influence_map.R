# 06_disease_influence_map.R
# ---------------------------------------------------------------------------
# Spatial maps of the disease-weighted Ailanthus Neighbourhood Index
# (ANI_disease) at the representative kernel (alpha = 2, beta = 0.15).
#
# Remakes the archived 09c "disease pressure" figure
# (archive/outputs/09c_disease_influence_all_plots.png), but for the CURRENT
# index used in the manuscript. The old figure mapped a home-grown index
# (score x BA / dist^2); this maps the published Gomez-Aparicio & Canham (2008)
# form modified to AUDPC-weight each tree:
#
#   ANI_disease(g, t) = sum_i  AUDPC_adj_i(t) * dbh_i^alpha * exp(-beta * d_ig)
#
# Timepoints match the manuscript's time-matching:
#   May       = pre-inoculation        -> field = 0 (no disease accumulated yet)
#   July      = 2 months post-inoc.    -> audpc_adj_2mai
#   September = 4 months post-inoc.    -> audpc_adj_4mai
#
# Layout: 3 plots (rows) x 3 timepoints (cols). Background = ANI_disease field;
# points = Ailanthus stems (size = DBH, fill = disease score at that timepoint);
# cyan diamonds = soil-sample sites (labelled with their ANI_disease value);
# dashed ring = 10 m reference; crosshairs = plot centre.
#
# Output: output/figures/disease_influence_map.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
  library(ggnewscale); library(ggrepel)
})

DATA_FILE  <- "bioassay_data_v_final.xlsx"
DIST_FLOOR <- 0.5
ALPHA      <- 2
BETA       <- 0.15
GRID_STEP  <- 0.2          # m
OUT_FIG    <- "output/figures"
if (!dir.exists(OUT_FIG)) dir.create(OUT_FIG, recursive = TRUE)

# display labels, ordered control -> attenuated -> virulent (as in the original)
TREEPLOT_TO_LABEL <- c(no_fungus_control = "No-fungus control",
                       attenuated        = "Attenuated strain",
                       virulent          = "Virulent strain")
TREAT_TO_LABEL    <- c(vnaa140_control = "No-fungus control",
                       vnaa140_2019    = "Attenuated strain",
                       vnaa140_2023    = "Virulent strain")
PLOT_LEVELS <- c("No-fungus control", "Attenuated strain", "Virulent strain")

polar_to_xy <- function(distance, bearing) {
  bearing  <- ifelse(is.na(bearing),  0, bearing)
  distance <- ifelse(is.na(distance), 0, distance)
  rad <- bearing * pi / 180
  list(x = distance * sin(rad), y = distance * cos(rad))
}

# timepoint -> tree disease-score column + AUDPC column (NA = pre-inoc, zero field)
tp <- tibble(
  label     = c("Pre-inoculation (May)",
                "2 mo. post-inoculation (July)",
                "4 mo. post-inoculation (September)"),
  score_col = c("disease_initial", "disease_2_mai", "disease_4_mai"),
  audpc_col = c(NA, "audpc_adj_2mai", "audpc_adj_4mai")
)
tp$label <- factor(tp$label, levels = tp$label)

# --- trees ------------------------------------------------------------------
trees <- read_excel(DATA_FILE, "tree_data") %>%
  mutate(tx = polar_to_xy(Distance_m, Bearing)$x,
         ty = polar_to_xy(Distance_m, Bearing)$y,
         plot_label = factor(unname(TREEPLOT_TO_LABEL[Plot]), levels = PLOT_LEVELS))

# --- soil-sample positions (distinct per plot) ------------------------------
soil <- read_excel(DATA_FILE, "bioassay_primary") %>%
  filter(treatment != "neg_control") %>%
  mutate(plot_label = factor(unname(TREAT_TO_LABEL[treatment]), levels = PLOT_LEVELS),
         sx = polar_to_xy(distance, bearing)$x,
         sy = polar_to_xy(distance, bearing)$y) %>%
  distinct(plot_label, sx, sy)

# --- index field on a grid, and at a point ----------------------------------
ani_field <- function(stems, gx, gy, audpc_col) {
  if (is.na(audpc_col)) return(rep(0, length(gx)))
  w  <- stems[[audpc_col]] * stems$dbh_cm^ALPHA
  dx <- outer(gx, stems$tx, `-`)
  dy <- outer(gy, stems$ty, `-`)
  d  <- pmax(sqrt(dx^2 + dy^2), DIST_FLOOR)
  as.numeric(exp(-BETA * d) %*% w)
}
ani_point <- function(stems, sx, sy, audpc_col) {
  if (is.na(audpc_col)) return(0)
  d <- pmax(sqrt((stems$tx - sx)^2 + (stems$ty - sy)^2), DIST_FLOOR)
  sum(stems[[audpc_col]] * stems$dbh_cm^ALPHA * exp(-BETA * d))
}

# --- map extent (fit all stems + soil + margin) -----------------------------
map_half <- ceiling(max(abs(c(trees$tx, trees$ty, soil$sx, soil$sy)))) + 1
grid_xy  <- expand.grid(x = seq(-map_half, map_half, by = GRID_STEP),
                        y = seq(-map_half, map_half, by = GRID_STEP))

# --- assemble field, stem and soil layers across plot x timepoint -----------
field_df <- list(); stem_df <- list(); soil_df <- list()
for (tpl in PLOT_LEVELS) {
  stems_p <- trees[trees$plot_label == tpl, ]
  soil_p  <- soil[soil$plot_label == tpl, ]
  for (i in seq_len(nrow(tp))) {
    ac <- tp$audpc_col[i]; sc <- tp$score_col[i]; lab <- tp$label[i]

    field_df[[length(field_df) + 1]] <- transform(
      grid_xy, ani = ani_field(stems_p, grid_xy$x, grid_xy$y, ac),
      plot_label = tpl, tpt = lab)

    sm <- stems_p[[sc]]; sm[is.na(sm)] <- 1            # untracked -> asymptomatic
    stem_df[[length(stem_df) + 1]] <- data.frame(
      x = stems_p$tx, y = stems_p$ty, dbh_cm = stems_p$dbh_cm,
      disease_score = sm, plot_label = tpl, tpt = lab)

    soil_df[[length(soil_df) + 1]] <- data.frame(
      x = soil_p$sx, y = soil_p$sy,
      ani = vapply(seq_len(nrow(soil_p)),
                   function(k) ani_point(stems_p, soil_p$sx[k], soil_p$sy[k], ac),
                   numeric(1)),
      plot_label = tpl, tpt = lab)
  }
}
field_all <- bind_rows(field_df); stem_all <- bind_rows(stem_df); soil_all <- bind_rows(soil_df)
field_all$plot_label <- factor(field_all$plot_label, levels = PLOT_LEVELS)
stem_all$plot_label  <- factor(stem_all$plot_label,  levels = PLOT_LEVELS)
soil_all$plot_label  <- factor(soil_all$plot_label,  levels = PLOT_LEVELS)

# 10 m reference ring (one copy per facet via recycling over panels)
theta   <- seq(0, 2 * pi, length.out = 181)
ring    <- data.frame(x = 10 * cos(theta), y = 10 * sin(theta))

cat(sprintf("Field max ANI_disease = %.0f | map half-width = %d m | grid step = %.2f m\n",
            max(field_all$ani), map_half, GRID_STEP))

# --- plot -------------------------------------------------------------------
p <- ggplot() +
  geom_raster(data = field_all, aes(x, y, fill = ani)) +
  scale_fill_viridis_c(option = "magma", name = expression(ANI[disease]),
                       limits = c(0, max(field_all$ani))) +
  geom_path(data = ring, aes(x, y), colour = "grey75",
            linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.2) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.2) +
  ggnewscale::new_scale_fill() +
  geom_point(data = stem_all,
             aes(x, y, size = dbh_cm, fill = disease_score),
             shape = 21, colour = "grey20", stroke = 0.25, alpha = 0.95) +
  scale_fill_gradient(low = "white", high = "#08306b",
                      limits = c(1, 6), breaks = 1:6,
                      name = "Tree disease\nscore (1–6)") +
  scale_size_continuous(range = c(1, 5), breaks = c(2, 5, 10, 20),
                        name = "DBH (cm)") +
  geom_point(data = soil_all, aes(x, y), shape = 23, size = 2.6,
             fill = "#17becf", colour = "black", stroke = 0.4) +
  geom_label_repel(data = soil_all, aes(x, y, label = round(ani)),
                   size = 2.4, colour = "black", fill = "white", alpha = 1,
                   label.padding = grid::unit(0.1, "lines"), label.size = 0.15,
                   segment.colour = "grey25", segment.size = 0.25,
                   min.segment.length = 0, box.padding = 0.45,
                   point.padding = 0.3, force = 2, max.overlaps = Inf, seed = 1) +
  facet_grid(plot_label ~ tpt) +
  coord_fixed(xlim = c(-map_half, map_half), ylim = c(-map_half, map_half)) +
  labs(
    x = "East of plot centre (m)", y = "North of plot centre (m)"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey90"),
        legend.position = "right")

fname <- file.path(OUT_FIG, "disease_influence_map.png")
ggsave(fname, p, width = 12, height = 10.5, dpi = 200)
cat("Wrote", fname, "\n")
