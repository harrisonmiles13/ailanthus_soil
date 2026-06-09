# 01_disease_influence.R
# ---------------------------------------------------------------------------
# Ailanthus Neighbourhood Index of DISEASE load for the soil bioassay.
#
# This is the Gomez-Aparicio & Canham (2008, J. Ecology 96:447-458) neighbourhood
# index (their Eqn 3), MODIFIED so that each Ailanthus tree's contribution is
# weighted by its disease load (AUDPC) rather than counting every tree equally.
# A tree that never developed disease contributes ZERO, so the index measures the
# spatial field of *diseased* Ailanthus to which each soil sample is exposed.
#
#   Published (Gomez-Aparicio & Canham 2008, Eqn 3):
#       ANI_j = sum_i  DBH_i^alpha * exp( -beta * d_ij )
#
#   This script (disease-weighted modification):
#       ANI_disease_j(t) = sum_i  AUDPC_i(t) * DBH_i^alpha * exp( -beta * d_ij )
#
#     AUDPC_i(t) = audpc_adj_*  (baseline-subtracted area under the disease
#                  progress curve; healthy tree => 0 contribution)
#     DBH_i^alpha = tree-size weight on the MAGNITUDE of influence
#                   (alpha = 0 -> size-blind disease density;
#                    alpha = 1 -> linear in diameter;
#                    alpha = 2 -> proportional to basal area)
#     exp(-beta*d) = Gomez-Aparicio distance-decay kernel; mean reach ~= 1/beta m
#     d_ij = Euclidean distance (m) between soil sample j and tree i, floored
#            at 0.5 m to cap the weight of a near-coincident tree
#
# NOTE on the difference from our earlier home-grown index: there, tree size set
# the decay LENGTH (exp(-lambda*d/S_i)); here, following the published model, size
# sets the MAGNITUDE (DBH_i^alpha) and the decay is a plain exp(-beta*d). alpha is
# FIXED (primary = 2, basal-area-proportional) and beta is SWEPT for sensitivity,
# because our limited disease gradient cannot support full MLE of the kernel the
# way the published 25-m, many-stand dataset could.
#
# Tree and soil-sample positions are recorded in polar form (distance + bearing
# from the plot centre) and converted to a shared Cartesian frame per plot.
#
# Output: output/tables/ani_disease.csv              (one row per soil sample)
#         output/tables/ani_disease_sensitivity.csv  (alpha x beta grid, 4mai)
# ---------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)

DATA_FILE   <- "bioassay_data_v_final.xlsx"
DIST_FLOOR  <- 0.5          # metres; minimum tree-to-sample distance
ALPHA       <- 2            # primary size exponent (DBH^2 ~ basal area)
BETA        <- 0.15         # primary decay constant (mean reach ~= 1/beta m)
AUDPC_COLS  <- c("audpc_adj_1mai", "audpc_adj_2mai",
                 "audpc_adj_3mai", "audpc_adj_4mai")

# bioassay treatment label  ->  tree_data Plot label
# (recovered from the shared pooled_sample key; neg_control has no Ailanthus trees)
TREATMENT_TO_PLOT <- c(
  vnaa140_2019    = "attenuated",
  vnaa140_2023    = "virulent",
  vnaa140_control = "no_fungus_control"
)

# --- polar (distance m, bearing deg from North) -> Cartesian (x = East, y = North)
polar_to_xy <- function(distance, bearing) {
  bearing  <- ifelse(is.na(bearing),  0, bearing)
  distance <- ifelse(is.na(distance), 0, distance)
  rad <- bearing * pi / 180
  list(x = distance * sin(rad), y = distance * cos(rad))
}

# === Load & prepare trees ===================================================
trees <- read_excel(DATA_FILE, sheet = "tree_data")

trees <- trees %>%
  mutate(
    BA_cm2 = dbh_cm^2 * 0.785,                 # basal area, pi*(dbh/2)^2 (cm^2)
    tx     = polar_to_xy(Distance_m, Bearing)$x,
    ty     = polar_to_xy(Distance_m, Bearing)$y
  )

cat(sprintf(
  "Trees: %d across %d plots | dbh %.1f..%.1f cm | max dist-from-centre %.1f m\n",
  nrow(trees), n_distinct(trees$Plot), min(trees$dbh_cm), max(trees$dbh_cm),
  max(sqrt(trees$tx^2 + trees$ty^2))))

# === Load & prepare soil samples ============================================
# Analysis is restricted to the three Ailanthus plots; neg_control (no Ailanthus)
# and the ailanthus_presence variable are excluded by design.
soil <- read_excel(DATA_FILE, sheet = "bioassay_primary") %>%
  filter(treatment != "neg_control") %>%
  mutate(
    plot = unname(TREATMENT_TO_PLOT[treatment]),
    sx   = polar_to_xy(distance, bearing)$x,
    sy   = polar_to_xy(distance, bearing)$y
  )

# === Disease-weighted neighbourhood index ===================================
# ANI_disease at one soil sample (sx, sy) in a given plot, for one AUDPC column.
# Sum runs over EVERY Ailanthus tree in the plot: the stem census itself is
# bounded (all stems <~11 m from plot centre), so it supplies the neighbourhood
# extent and no explicit radius cap is needed (cf. the published 25-m window).
ani_at <- function(plot, sx, sy, audpc_col, alpha = ALPHA, beta = BETA) {
  if (is.na(plot)) return(0)                       # neg_control: no Ailanthus
  pt <- trees[trees$Plot == plot, ]
  d  <- sqrt((pt$tx - sx)^2 + (pt$ty - sy)^2)
  d  <- pmax(d, DIST_FLOOR)
  sum(pt[[audpc_col]] * pt$dbh_cm^alpha * exp(-beta * d))
}

# Comparator: total basal area within 5 m (a model-free density proxy) --------
ba_within <- function(plot, sx, sy, radius = 5) {
  if (is.na(plot)) return(0)
  pt <- trees[trees$Plot == plot, ]
  d  <- sqrt((pt$tx - sx)^2 + (pt$ty - sy)^2)
  sum(pt$BA_cm2[d <= radius])
}

# === Compute index for every soil sample & timepoint ========================
# One column per AUDPC timepoint, at the primary (alpha, beta).
for (col in AUDPC_COLS) {
  out_col <- sub("audpc_adj", "ani_disease", col)
  soil[[out_col]] <- mapply(
    function(p, x, y) ani_at(p, x, y, col, ALPHA, BETA),
    soil$plot, soil$sx, soil$sy
  )
}
soil$ba_sum_5m <- mapply(function(p, x, y) ba_within(p, x, y, 5),
                         soil$plot, soil$sx, soil$sy)

ani_out <- soil %>%
  select(pooled_sample, treatment, plot, month,
         distance, bearing, sx, sy, ba_sum_5m,
         starts_with("ani_disease_"))

write.csv(ani_out, "output/tables/ani_disease.csv", row.names = FALSE)

# Pooled sanity check only. By-treatment / plot-mean summaries are intentionally
# omitted: the 3 plots are confounded with the 3 treatments, so a per-treatment
# table invites a pseudoreplicated between-plot reading. Treatment is carried in
# the output for joining, not for comparison.
cat(sprintf(
  "\nANI_disease (4mai, alpha=%g beta=%g) across the %d Ailanthus samples: mean %.1f | sd %.1f | range %.1f..%.1f\n",
  ALPHA, BETA, nrow(ani_out), mean(ani_out$ani_disease_4mai),
  sd(ani_out$ani_disease_4mai),
  min(ani_out$ani_disease_4mai),
  max(ani_out$ani_disease_4mai)))

# === alpha x beta sensitivity (primary index = 4mai) ========================
# Reports, for each kernel choice, the index spread and how closely it RANKS the
# 45 samples relative to the primary index (Spearman) and to the model-free BA
# density. A high rank correlation means the downstream regression is insensitive
# to the exact kernel constants.
alpha_grid <- c(0, 1, 2)
beta_grid  <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.50)

ail      <- !is.na(soil$plot)
ref_vals <- ani_out$ani_disease_4mai            # primary (alpha=2, beta=0.15)

sensitivity <- map_dfr(alpha_grid, function(a) {
  map_dfr(beta_grid, function(b) {
    vals <- mapply(function(p, x, y) ani_at(p, x, y, "audpc_adj_4mai", a, b),
                   soil$plot, soil$sx, soil$sy)
    tibble(
      alpha          = a,
      beta           = b,
      mean_reach_m   = 1 / b,
      mean_index     = mean(vals[ail]),
      sd_index       = sd(vals[ail]),
      rho_primary    = suppressWarnings(cor(vals[ail], ref_vals[ail],
                                            method = "spearman")),
      rho_ba_5m      = suppressWarnings(cor(vals[ail], soil$ba_sum_5m[ail],
                                            method = "spearman"))
    )
  })
})

write.csv(sensitivity, "output/tables/ani_disease_sensitivity.csv",
          row.names = FALSE)

cat("\nalpha x beta sensitivity (4mai index; rho = Spearman rank corr):\n")
as.data.frame(sensitivity) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>% print(row.names = FALSE)

cat("\nWrote output/tables/ani_disease.csv (",
    nrow(ani_out), " samples )\n", sep = "")
