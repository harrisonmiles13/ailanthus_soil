# 09_satellite_map.R
# ---------------------------------------------------------------------------
# Satellite-imagery site map: the three Ailanthus stands, every mapped stem,
# and every soil-collection point.
#
# Trees (tree_data) and soil samples (bioassay_primary) are recorded as POLAR
# offsets -- compass bearing + distance in metres -- from a stand centre. The
# only absolute positions in the study are the three stand centres:
#
#   Formerly attenuated strain   37.31776, -76.88626
#   Virulent strain              37.31766, -76.88659
#   No-fungus control            37.31725, -76.88538
#
# This script projects those offsets onto WGS84 and writes two versions of
# each panel:
#
#   output/figures/satellite_map_overview.png         on Esri World Imagery
#   output/figures/satellite_map_<stand>.png
#   output/figures/satellite_map_overview_blank.png    transparent background,
#   output/figures/satellite_map_<stand>_blank.png     no imagery at all
#
# The "_blank" figures carry only the geometry (stems, soil points, sampling
# rings, a scale bar) on a transparent background, with no basemap fetch and
# no dependency on maptiles/tidyterra. They exist so the point cloud can be
# hand-registered, in image-editing or GIS software, against an independent
# photo of the site from the study period -- useful whenever a basemap
# provider's cached imagery is a different vintage than the vegetation being
# described (e.g. it shows bare ground for a period when the plots were
# forested, or vice versa). Toggle RENDER_BASEMAP_MAPS / RENDER_BLANK_MAPS
# below to skip either set.
#
# The bearing convention matches polar_to_xy() in R/06_disease_influence_map.R
# (x = d*sin(theta) east, y = d*cos(theta) north, theta clockwise from north).
# DECLINATION_DEG rotates every bearing before projection and defaults to 0,
# i.e. bearings are treated as TRUE north exactly as the ANI analyses treat
# them. If the field compass was uncorrected magnetic, set it to the local
# declination (about -11 for Williamsburg VA in 2023) to swing the whole point
# cloud onto true north.
#
# NOTE: get_tiles() downloads imagery, so this script needs a network
# connection the first time it runs (tiles are then cached under CACHE_DIR).
# ---------------------------------------------------------------------------

RENDER_BASEMAP_MAPS <- TRUE     # satellite_map_*.png       (needs maptiles, tidyterra, network)
RENDER_BLANK_MAPS   <- TRUE     # satellite_map_*_blank.png (geometry only, no extra deps)

need <- c("readxl", "dplyr", "ggplot2", "tibble", "sf")
if (RENDER_BASEMAP_MAPS) need <- c(need, "maptiles", "tidyterra")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing))
  stop("Install the missing packages first:\n  install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))", call. = FALSE)

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(ggplot2)
})

DATA_FILE       <- "bioassay_data_v_final.xlsx"
OUT_FIG         <- "output/figures"
CACHE_DIR       <- "output/spatial/tile_cache"
DECLINATION_DEG <- 0
SAMPLING_RADIUS <- 10          # m, radius of the soil-sampling design
TILE_PROVIDER   <- "Esri.WorldImagery"

# every path here is relative to the repository root
if (!file.exists(DATA_FILE))
  stop("Cannot find ", DATA_FILE, " in the working directory (", getwd(), ").\n",
       "Run this from the repository root, e.g.\n",
       "  Rscript R/09_satellite_map.R", call. = FALSE)

dir.create(OUT_FIG,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

STAND_CENTRES <- tibble::tribble(
  ~stand,               ~lat,      ~lon,
  "no_fungus_control", 37.31725, -76.88538,
  "attenuated",        37.31776, -76.88626,
  "virulent",          37.31766, -76.88659
)

TREAT_TO_STAND <- c(vnaa140_control = "no_fungus_control",
                    vnaa140_2019    = "attenuated",
                    vnaa140_2023    = "virulent")
STAND_LABEL <- c(no_fungus_control = "No-fungus control",
                 attenuated        = "Formerly attenuated strain",
                 virulent          = "Virulent strain")
STAND_LEVELS <- unname(STAND_LABEL)

# --- polar offset -> local metres -> WGS84 ----------------------------------
polar_to_xy <- function(distance, bearing) {
  bearing  <- ifelse(is.na(bearing),  0, bearing)
  distance <- ifelse(is.na(distance), 0, distance)
  rad <- (bearing + DECLINATION_DEG) * pi / 180
  list(x = distance * sin(rad), y = distance * cos(rad))
}

# metres per degree of latitude / longitude at a given latitude (WGS84 series)
m_per_deg <- function(lat) {
  phi <- lat * pi / 180
  list(lat = 111132.92 - 559.82 * cos(2 * phi) + 1.175 * cos(4 * phi) -
              0.0023 * cos(6 * phi),
       lon = 111412.84 * cos(phi) - 93.5 * cos(3 * phi) + 0.118 * cos(5 * phi))
}

offset_to_lonlat <- function(lat0, lon0, east, north) {
  m <- m_per_deg(lat0)
  list(lon = lon0 + east / m$lon, lat = lat0 + north / m$lat)
}

georef <- function(df, stand_col = "stand") {
  centres <- STAND_CENTRES[match(df[[stand_col]], STAND_CENTRES$stand), ]
  ll <- offset_to_lonlat(centres$lat, centres$lon, df$east_m, df$north_m)
  df$lon <- ll$lon
  df$lat <- ll$lat
  df
}

# --- trees ------------------------------------------------------------------
trees <- read_excel(DATA_FILE, "tree_data") %>%
  mutate(stand  = Plot,
         east_m = polar_to_xy(Distance_m, Bearing)$x,
         north_m = polar_to_xy(Distance_m, Bearing)$y,
         # untracked stems are asymptomatic, as in 06_disease_influence_map.R
         disease_sep = ifelse(is.na(disease_4_mai), 1, disease_4_mai)) %>%
  georef() %>%
  mutate(stand_label = factor(unname(STAND_LABEL[stand]), levels = STAND_LEVELS))

# --- soil-collection points (one row per distinct position) -----------------
soil <- read_excel(DATA_FILE, "bioassay_primary") %>%
  filter(treatment != "neg_control") %>%      # off-site soils have no coordinates
  mutate(stand  = unname(TREAT_TO_STAND[treatment]),
         east_m = polar_to_xy(distance, bearing)$x,
         north_m = polar_to_xy(distance, bearing)$y) %>%
  distinct(stand, bearing, distance, east_m, north_m) %>%
  georef() %>%
  mutate(stand_label = factor(unname(STAND_LABEL[stand]), levels = STAND_LEVELS))

centres <- STAND_CENTRES %>%
  mutate(stand_label = factor(unname(STAND_LABEL[stand]), levels = STAND_LEVELS))

cat(sprintf("%d stems, %d soil points, %d stands\n",
            nrow(trees), nrow(soil), nrow(centres)))

# --- sf geometry ------------------------------------------------------------
as_pts <- function(df) sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
trees_sf   <- as_pts(trees)
soil_sf    <- as_pts(soil)
centres_sf <- as_pts(centres)

rings_sf <- centres_sf %>%
  sf::st_transform(3857) %>%
  sf::st_buffer(SAMPLING_RADIUS / cos(mean(STAND_CENTRES$lat) * pi / 180)) %>%
  sf::st_transform(4326)          # metres are Web-Mercator-inflated; undo by latitude

# --- shared ground extent ----------------------------------------------------
# `pad` is the margin in ground metres added around the features. Returns the
# WGS84 extent (as an sfc polygon), its bbox, and the width:height of the
# ground it covers (so a saved figure can be sized to fill its panel).
compute_extent <- function(tr, so, ce, pad) {
  feats  <- c(sf::st_geometry(tr), sf::st_geometry(so), sf::st_geometry(ce))
  extent <- sf::st_as_sfc(sf::st_bbox(feats)) %>%
    sf::st_transform(3857) %>%
    # pad is in ground metres; Web Mercator inflates them by 1/cos(latitude)
    sf::st_buffer(pad / cos(mean(STAND_CENTRES$lat) * pi / 180)) %>%
    sf::st_transform(4326)
  bb  <- sf::st_bbox(extent)
  mid <- m_per_deg(mean(c(bb[["ymin"]], bb[["ymax"]])))
  ground_aspect <- ((bb[["xmax"]] - bb[["xmin"]]) * mid$lon) /
                   ((bb[["ymax"]] - bb[["ymin"]]) * mid$lat)
  list(extent = extent, bb = bb, ground_aspect = ground_aspect)
}

# --- shared data layers: stems, soil points, sampling ring -------------------
# Disease score rides `colour`, not `fill`, for the same reason as
# 06_disease_influence_map_v2.R: on the imagery panel the basemap already owns
# a fill scale and ggnewscale is not a dependency here; kept the same way on
# the blank panel so the two stay visually consistent. A fixed grey ring on
# top restores the marker definition that shape 21 would have given.
feature_layers <- function(tr, so, ce, ri, point_range, ring_colour, ring = TRUE) {
  layers <- list()
  if (ring)
    layers <- c(layers, list(geom_sf(data = ri, fill = NA, colour = ring_colour,
                                     linewidth = 0.4, linetype = "22")))
  c(layers, list(
    geom_sf(data = tr, aes(size = dbh_cm, colour = disease_sep),
            shape = 16, alpha = 0.95),
    geom_sf(data = tr, aes(size = dbh_cm), shape = 21, fill = NA,
            colour = "grey15", stroke = 0.3, show.legend = FALSE),
    geom_sf(data = so, shape = 23, size = 3.2, fill = "#17becf",
            colour = "black", stroke = 0.6),
    geom_sf(data = ce, shape = 3, size = 3.4, colour = ring_colour, stroke = 0.9),
    scale_colour_gradient(low = "white", high = "#08306b",
                          limits = c(1, 6), breaks = 1:6,
                          name = "Stem disease score\n(September, 1-6)"),
    scale_size_continuous(range = point_range, breaks = c(2, 5, 10, 20, 30),
                          name = "Stem DBH (cm)"),
    guides(colour = guide_colourbar(order = 1, barheight = grid::unit(6, "lines")),
           size   = guide_legend(order = 2))
  ))
}

# --- imagery panel ------------------------------------------------------------
map_panel <- function(tr, so, ce, ri, pad, zoom, point_range, ring = TRUE) {
  ext <- compute_extent(tr, so, ce, pad)
  bb  <- ext$bb

  # imagery is not published at the same maximum zoom everywhere; step down
  # until the provider actually serves tiles for this extent.
  tiles <- NULL
  for (z in seq(zoom, 17)) {
    tiles <- tryCatch(
      maptiles::get_tiles(ext$extent, provider = TILE_PROVIDER, zoom = z,
                          crop = TRUE, cachedir = CACHE_DIR,
                          forceDownload = FALSE),
      error = function(e) NULL)
    if (!is.null(tiles)) { if (z != zoom) cat("  (fell back to zoom", z, ")\n"); break }
  }
  if (is.null(tiles))
    stop("Could not fetch ", TILE_PROVIDER, " tiles -- check the network connection.",
         call. = FALSE)

  p <- ggplot() +
    tidyterra::geom_spatraster_rgb(data = tiles, maxcell = 5e6) +
    feature_layers(tr, so, ce, ri, point_range, ring_colour = "white", ring = ring) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE, crs = 4326) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_line(colour = alpha("white", 0.18), linewidth = 0.25),
          axis.text = element_text(size = rel(0.7), colour = "grey30"),
          legend.position = "right",
          plot.margin = margin(6, 8, 6, 6))

  structure(p, ground_aspect = ext$ground_aspect)
}

# --- blank panel: geometry only, transparent background, no basemap ----------
# For hand-registering this point cloud against an independent site photo
# (e.g. taken near the study dates) in image-editing or GIS software. Needs
# only sf/ggplot2 -- no maptiles, no tidyterra, no network access.
blank_panel <- function(tr, so, ce, ri, pad, point_range, scale_bar_m = 5) {
  ext <- compute_extent(tr, so, ce, pad)
  bb  <- ext$bb
  mid <- m_per_deg(mean(c(bb[["ymin"]], bb[["ymax"]])))

  # a short reference line of known ground length in the bottom-left margin,
  # for scaling (and, from its bearing, rotating) this layer onto a photo
  # that has no coordinate reference of its own
  bar_x0 <- bb[["xmin"]] + 0.10 * (bb[["xmax"]] - bb[["xmin"]])
  bar_y0 <- bb[["ymin"]] + 0.06 * (bb[["ymax"]] - bb[["ymin"]])
  bar    <- data.frame(x = c(bar_x0, bar_x0 + scale_bar_m / mid$lon), y = bar_y0)

  p <- ggplot() +
    feature_layers(tr, so, ce, ri, point_range, ring_colour = "grey40") +
    geom_line(data = bar, aes(x, y), linewidth = 1, colour = "grey15",
             inherit.aes = FALSE) +
    geom_text(data = bar[1, ], aes(x, y, label = paste(scale_bar_m, "m")),
             colour = "grey15", size = 3.2, vjust = 2.1, hjust = 0,
             inherit.aes = FALSE) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE, crs = 4326) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 13) +
    theme(legend.position = "right",
          legend.background = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          plot.margin = margin(6, 8, 6, 6))

  structure(p, ground_aspect = ext$ground_aspect)
}

# Save at a height that matches the panel's own aspect, so the figure is not
# padded with empty bands. legend_in / extra_in approximate the non-panel
# furniture (colourbar + size legend, axis text, caption).
save_panel <- function(p, file, width, legend_in, extra_in = 1.0, transparent = FALSE) {
  height <- (width - legend_in) / attr(p, "ground_aspect") + extra_in
  ggsave(file, p, width = width, height = height, dpi = 300,
        bg = if (transparent) "transparent" else "white")
  cat(sprintf("Wrote %s  (%.1f x %.1f in @ 300 dpi)\n", file, width, height))
}

# --- overview: all three stands ---------------------------------------------
if (RENDER_BASEMAP_MAPS) {
  ov <- map_panel(trees_sf, soil_sf, centres_sf, rings_sf,
                  pad = 15, zoom = 20, point_range = c(0.7, 3.2)) +
    geom_label(data = centres, aes(x = lon, y = lat + 0.00018, label = stand_label),
               colour = "white", fill = "grey10", alpha = 0.62, size = 3.6,
               fontface = "bold", label.size = 0, label.r = grid::unit(0.12, "lines"),
               label.padding = grid::unit(0.22, "lines"), inherit.aes = FALSE) +
    labs(caption = sprintf(
      "Esri World Imagery. %d mapped Ailanthus stems, %d soil-collection points (stand centre + 10 m N/E/S/W).",
      nrow(trees_sf), nrow(soil_sf)))

  save_panel(ov, file.path(OUT_FIG, "satellite_map_overview.png"),
             width = 12, legend_in = 3.0)
}

if (RENDER_BLANK_MAPS) {
  ovb <- blank_panel(trees_sf, soil_sf, centres_sf, rings_sf,
                     pad = 15, point_range = c(0.7, 3.2)) +
    geom_text(data = centres, aes(x = lon, y = lat + 0.00018, label = stand_label),
              colour = "grey10", size = 3.6, fontface = "bold", inherit.aes = FALSE) +
    labs(caption = sprintf(
      "No basemap -- register against an independent site photo. %d mapped Ailanthus stems, %d soil-collection points.",
      nrow(trees_sf), nrow(soil_sf)))

  save_panel(ovb, file.path(OUT_FIG, "satellite_map_overview_blank.png"),
             width = 12, legend_in = 3.0, transparent = TRUE)
}

# --- one zoomed panel per stand ---------------------------------------------
for (s in STAND_CENTRES$stand) {
  lab <- unname(STAND_LABEL[s])
  keep <- function(x) x[x$stand_label == lab, ]
  tr <- keep(trees_sf); so <- keep(soil_sf); ce <- keep(centres_sf); ri <- keep(rings_sf)

  if (RENDER_BASEMAP_MAPS) {
    p <- map_panel(tr, so, ce, ri, pad = 4, zoom = 20, point_range = c(1.4, 6)) +
      labs(caption = sprintf("%s -- Esri World Imagery; dashed ring = %d m soil-sampling radius.",
                             lab, SAMPLING_RADIUS))
    save_panel(p, file.path(OUT_FIG, sprintf("satellite_map_%s.png", s)),
               width = 10, legend_in = 3.0)
  }

  if (RENDER_BLANK_MAPS) {
    pb <- blank_panel(tr, so, ce, ri, pad = 4, point_range = c(1.4, 6)) +
      labs(caption = sprintf(
        "%s -- no basemap; dashed ring = %d m soil-sampling radius.", lab, SAMPLING_RADIUS))
    save_panel(pb, file.path(OUT_FIG, sprintf("satellite_map_%s_blank.png", s)),
               width = 10, legend_in = 3.0, transparent = TRUE)
  }
}
