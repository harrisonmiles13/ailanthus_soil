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

TREATMENT_TO_PLOT <- c(
  vnaa140_2019    = "attenuated",
  vnaa140_2023    = "virulent",
  vnaa140_control = "no_fungus_control"
)

polar_to_xy <- function(distance, bearing) {
  bearing  <- ifelse(is.na(bearing),  0, bearing)
  distance <- ifelse(is.na(distance), 0, distance)
  rad <- bearing * pi / 180
  list(x = distance * sin(rad), y = distance * cos(rad))
}

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

soil <- read_excel(DATA_FILE, sheet = "bioassay_primary") %>%
  filter(treatment != "neg_control") %>%
  mutate(
    plot = unname(TREATMENT_TO_PLOT[treatment]),
    sx   = polar_to_xy(distance, bearing)$x,
    sy   = polar_to_xy(distance, bearing)$y
  )

ani_at <- function(plot, sx, sy, audpc_col, alpha = ALPHA, beta = BETA) {
  if (is.na(plot)) return(0)                      
  pt <- trees[trees$Plot == plot, ]
  d  <- sqrt((pt$tx - sx)^2 + (pt$ty - sy)^2)
  d  <- pmax(d, DIST_FLOOR)
  sum(pt[[audpc_col]] * pt$dbh_cm^alpha * exp(-beta * d))
}

ba_within <- function(plot, sx, sy, radius = 5) {
  if (is.na(plot)) return(0)
  pt <- trees[trees$Plot == plot, ]
  d  <- sqrt((pt$tx - sx)^2 + (pt$ty - sy)^2)
  sum(pt$BA_cm2[d <= radius])
}

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

cat(sprintf(
  "\nANI_disease (4mai, alpha=%g beta=%g) across the %d Ailanthus samples: mean %.1f | sd %.1f | range %.1f..%.1f\n",
  ALPHA, BETA, nrow(ani_out), mean(ani_out$ani_disease_4mai),
  sd(ani_out$ani_disease_4mai),
  min(ani_out$ani_disease_4mai),
  max(ani_out$ani_disease_4mai)))

alpha_grid <- c(0, 1, 2)
beta_grid  <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.50)

ail      <- !is.na(soil$plot)
ref_vals <- ani_out$ani_disease_4mai            

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
