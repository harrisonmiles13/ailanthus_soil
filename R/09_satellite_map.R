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
# This script projects those offsets onto WGS84, drops them on Esri World
# Imagery, and writes:
#
#   output/figures/satellite_map_overview.png    all three stands together
#   output/figures/satellite_map_<stand>.png     one zoomed panel per stand
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

need <- c("readxl", "dplyr", "ggplot2", "tibble", "sf", "maptiles", "tidyterra")
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
# Okabe-Ito, matching PLOT_COLOURS in R/v2_style.R
STAND_COLOUR <- c("No-fungus control"          = "#0072B2",
                  "Formerly attenuated strain" = "#D55E00",
                  "Virulent strain"            = "#009E73")

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

# stand outline = convex hull of the mapped stems; ring = 10 m sampling radius
hulls_sf <- do.call(rbind, lapply(STAND_LEVELS, function(lab) {
  stems <- sf::st_geometry(trees_sf[trees_sf$stand_label == lab, ])
  sf::st_sf(stand_label = factor(lab, levels = STAND_LEVELS),
            geometry    = sf::st_convex_hull(sf::st_combine(stems)))
}))

rings_sf <- centres_sf %>%
  sf::st_transform(3857) %>%
  sf::st_buffer(SAMPLING_RADIUS / cos(mean(STAND_CENTRES$lat) * pi / 180)) %>%
  sf::st_transform(4326)          # metres are Web-Mercator-inflated; undo by latitude

# --- one map panel ----------------------------------------------------------
# `pad` is the margin in metres added around the features before fetching tiles.
map_panel <- function(tr, so, ce, hu, ri, pad, zoom, point_range, ring = TRUE) {
  feats  <- c(sf::st_geometry(tr), sf::st_geometry(so), sf::st_geometry(ce))
  extent <- sf::st_as_sfc(sf::st_bbox(feats)) %>%
    sf::st_transform(3857) %>%
    # pad is in ground metres; Web Mercator inflates them by 1/cos(latitude)
    sf::st_buffer(pad / cos(mean(STAND_CENTRES$lat) * pi / 180)) %>%
    sf::st_transform(4326)

  # imagery is not published at the same maximum zoom everywhere; step down
  # until the provider actually serves tiles for this extent.
  tiles <- NULL
  for (z in seq(zoom, 17)) {
    tiles <- tryCatch(
      maptiles::get_tiles(extent, provider = TILE_PROVIDER, zoom = z,
                          crop = TRUE, cachedir = CACHE_DIR,
                          forceDownload = FALSE),
      error = function(e) NULL)
    if (!is.null(tiles)) { if (z != zoom) cat("  (fell back to zoom", z, ")\n"); break }
  }
  if (is.null(tiles))
    stop("Could not fetch ", TILE_PROVIDER, " tiles -- check the network connection.",
         call. = FALSE)
  bb <- sf::st_bbox(extent)

  p <- ggplot() +
    tidyterra::geom_spatraster_rgb(data = tiles, maxcell = 5e6) +
    geom_sf(data = hu, aes(colour = stand_label), fill = NA,
            linewidth = 0.7, show.legend = FALSE)

  if (ring)
    p <- p + geom_sf(data = ri, fill = NA, colour = "white",
                     linewidth = 0.4, linetype = "22")

  p +
    geom_sf(data = tr, aes(size = dbh_cm, fill = disease_sep),
            shape = 21, colour = "grey15", stroke = 0.25, alpha = 0.95) +
    geom_sf(data = so, shape = 23, size = 3.2, fill = "#17becf",
            colour = "black", stroke = 0.6) +
    geom_sf(data = ce, shape = 3, size = 3.4, colour = "white", stroke = 0.9) +
    scale_fill_gradient(low = "white", high = "#08306b",
                        limits = c(1, 6), breaks = 1:6,
                        name = "Stem disease score\n(September, 1-6)") +
    scale_colour_manual(values = STAND_COLOUR) +
    scale_size_continuous(range = point_range, breaks = c(2, 5, 10, 20, 30),
                          name = "Stem DBH (cm)") +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE, crs = 4326) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colourbar(order = 1, barheight = grid::unit(6, "lines")),
           size = guide_legend(order = 2)) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_line(colour = alpha("white", 0.18), linewidth = 0.25),
          axis.text = element_text(size = rel(0.7), colour = "grey30"),
          legend.position = "right",
          plot.margin = margin(6, 8, 6, 6))
}

# --- overview: all three stands ---------------------------------------------
ov <- map_panel(trees_sf, soil_sf, centres_sf, hulls_sf, rings_sf,
                pad = 35, zoom = 20, point_range = c(0.7, 3.2)) +
  geom_sf_text(data = centres_sf, aes(label = stand_label),
               colour = "white", size = 3.6, fontface = "bold",
               nudge_y = 0.00016) +
  labs(caption = sprintf(
    "Esri World Imagery. %d mapped Ailanthus stems, %d soil-collection points (stand centre + 10 m N/E/S/W).",
    nrow(trees_sf), nrow(soil_sf)))

f <- file.path(OUT_FIG, "satellite_map_overview.png")
ggsave(f, ov, width = 12, height = 8.5, dpi = 300)
cat("Wrote", f, "\n")

# --- one zoomed panel per stand ---------------------------------------------
for (s in STAND_CENTRES$stand) {
  lab <- unname(STAND_LABEL[s])
  keep <- function(x) x[x$stand_label == lab, ]
  p <- map_panel(keep(trees_sf), keep(soil_sf), keep(centres_sf),
                 keep(hulls_sf), keep(rings_sf),
                 pad = 12, zoom = 20, point_range = c(1.4, 6)) +
    labs(caption = sprintf("%s -- Esri World Imagery; dashed ring = %d m soil-sampling radius.",
                           lab, SAMPLING_RADIUS))
  f <- file.path(OUT_FIG, sprintf("satellite_map_%s.png", s))
  ggsave(f, p, width = 9, height = 8, dpi = 300)
  cat("Wrote", f, "\n")
}
