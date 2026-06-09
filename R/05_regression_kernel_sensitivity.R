# 05_regression_kernel_sensitivity.R
# ---------------------------------------------------------------------------
# Carries the alpha x beta kernel sweep THROUGH the within-plot regression, so
# robustness is shown at the level of the inference (slope / P / dAICc), not just
# the index ranking (which 01_disease_influence.R already covers via Spearman).
#
# For every kernel in the grid the time-matched index is rebuilt, standardised
# per-SD, and all 12 species x response models are refit exactly as in
# 02_bioassay_regression.R. The conclusion's stability across the grid is then
# summarised.
#
# Self-contained on purpose: it mirrors the index formula of 01 (ani_at) and the
# fit logic of 02 (aicc / qaicc / fit_one) rather than sourcing them (sourcing
# those scripts would re-run their file writes). The alpha = 2, beta = 0.15 cell
# MUST reproduce the committed output/tables/bioassay_regression.csv -- this is
# asserted at the end as a sanity anchor.
#
# Output: output/tables/bioassay_regression_sensitivity.csv  (18 kernels x 12 models)
# ---------------------------------------------------------------------------

suppressMessages({ library(readxl); library(dplyr) })

DATA_FILE  <- "bioassay_data_v_final.xlsx"
DIST_FLOOR <- 0.5
SEEDS      <- c(bes = 10, jg = 100, wsr = 40, yar = 40)
SPECIES    <- names(SEEDS)
ALPHA_GRID <- c(0, 1, 2)
BETA_GRID  <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.50)

TREATMENT_TO_PLOT <- c(vnaa140_2019    = "attenuated",
                       vnaa140_2023    = "virulent",
                       vnaa140_control = "no_fungus_control")

# --- geometry (mirrors 01_disease_influence.R) ------------------------------
polar_to_xy <- function(distance, bearing) {
  bearing  <- ifelse(is.na(bearing),  0, bearing)
  distance <- ifelse(is.na(distance), 0, distance)
  rad <- bearing * pi / 180
  list(x = distance * sin(rad), y = distance * cos(rad))
}

trees <- read_excel(DATA_FILE, "tree_data") %>%
  mutate(tx = polar_to_xy(Distance_m, Bearing)$x,
         ty = polar_to_xy(Distance_m, Bearing)$y)

# bioassay_primary carries BOTH sample positions and the per-tray responses
soil <- read_excel(DATA_FILE, "bioassay_primary") %>%
  filter(treatment != "neg_control") %>%
  mutate(plot = unname(TREATMENT_TO_PLOT[treatment]),
         sx   = polar_to_xy(distance, bearing)$x,
         sy   = polar_to_xy(distance, bearing)$y)

# --- index at one sample for one AUDPC column (mirrors 01:ani_at) -----------
ani_at <- function(plot, sx, sy, audpc_col, alpha, beta) {
  pt <- trees[trees$Plot == plot, ]
  d  <- sqrt((pt$tx - sx)^2 + (pt$ty - sy)^2)
  d  <- pmax(d, DIST_FLOOR)
  sum(pt[[audpc_col]] * pt$dbh_cm^alpha * exp(-beta * d))
}

# --- information criteria (mirrors 02_bioassay_regression.R) -----------------
aicc <- function(m) {
  ll <- logLik(m); k <- attr(ll, "df"); n <- nobs(m)
  as.numeric(-2 * ll + 2 * k + (2 * k * (k + 1)) / (n - k - 1))
}
qaicc <- function(m, chat, n) {
  ll <- logLik(m); k <- attr(ll, "df") + 1
  as.numeric(-2 * as.numeric(ll) / chat + 2 * k + (2 * k * (k + 1)) / (n - k - 1))
}

# --- per-species fit (mirrors 02:fit_one; dat passed in explicitly) ---------
fit_one <- function(dat, y_full, label, species_name, family) {
  if (family == "gaussian") {
    df <- data.frame(y = y_full, z_ani = dat$z_ani, plot = dat$plot)
    df <- df[!is.na(df$y), ]
    m0 <- lm(y ~ plot,         data = df)
    m1 <- lm(y ~ z_ani + plot, data = df)
    co <- summary(m1)$coefficients["z_ani", ]
    slope <- co[["Estimate"]]; se <- co[["Std. Error"]]
    stat  <- co[["t value"]];  p  <- co[["Pr(>|t|)"]]
    disp  <- NA_real_
    dic   <- aicc(m0) - aicc(m1)
  } else {                                     # germination: binomial counts
    n_seed <- SEEDS[[species_name]]
    germ   <- round(y_full * n_seed)
    df <- data.frame(germ = germ, fail = n_seed - germ,
                     z_ani = dat$z_ani, plot = dat$plot)
    df <- df[!is.na(df$germ), ]
    n  <- nrow(df)
    b0 <- glm(cbind(germ, fail) ~ plot,         family = binomial, data = df)
    b1 <- glm(cbind(germ, fail) ~ z_ani + plot, family = binomial, data = df)
    disp <- sum(residuals(b1, "pearson")^2) / df.residual(b1)
    chat <- max(disp, 1)
    dic  <- qaicc(b0, chat, n) - qaicc(b1, chat, n)
    q1 <- glm(cbind(germ, fail) ~ z_ani + plot, family = quasibinomial, data = df)
    co <- summary(q1)$coefficients["z_ani", ]
    slope <- co[["Estimate"]]; se <- co[["Std. Error"]]
    stat  <- co[["t value"]];  p  <- co[["Pr(>|t|)"]]
    family <- "quasibinomial"
  }
  data.frame(species = species_name, response = label, family = family,
             n = nrow(df), slope_per_SD = slope, se = se, statistic = stat,
             p = p, dAICc = dic, dispersion = disp)
}

# === sweep ==================================================================
res <- data.frame()
for (a in ALPHA_GRID) for (b in BETA_GRID) {
  i2 <- mapply(function(p, x, y) ani_at(p, x, y, "audpc_adj_2mai", a, b),
               soil$plot, soil$sx, soil$sy)
  i4 <- mapply(function(p, x, y) ani_at(p, x, y, "audpc_adj_4mai", a, b),
               soil$plot, soil$sx, soil$sy)
  ani_tm <- ifelse(soil$month == "may", 0,
                   ifelse(soil$month == "july", i2, i4))     # sep -> 4mai
  dat <- soil
  dat$plot  <- factor(dat$plot)
  dat$z_ani <- as.numeric(scale(ani_tm))
  for (s in SPECIES) {
    res <- rbind(res,
      cbind(alpha = a, beta = b,
            fit_one(dat, dat[[paste0(s, "_total_avg")]], "length",      s, "gaussian")),
      cbind(alpha = a, beta = b,
            fit_one(dat, dat[[paste0(s, "_ratio_avg")]], "ratio",       s, "gaussian")),
      cbind(alpha = a, beta = b,
            fit_one(dat, dat[[paste0(s, "_germ_perc")]], "germination", s, "binomial")))
  }
}

write.csv(res, "output/tables/bioassay_regression_sensitivity.csv", row.names = FALSE)

# === summary across the 18 kernels, per species x response ==================
summ <- res %>%
  group_by(species, response) %>%
  summarise(p_min = min(p), p_max = max(p),
            dAICc_min = min(dAICc), dAICc_max = max(dAICc),
            favours_disease_all = all(dAICc > 0), .groups = "drop") %>%
  as.data.frame()
cat("\n=== Kernel sweep: range of each model's result across 18 kernels ===\n")
print(summ, row.names = FALSE, digits = 3)

cat(sprintf("\nGlobal minimum P across all 12 models x 18 kernels = %.4f (any < 0.05? %s)\n",
            min(res$p), ifelse(min(res$p) < 0.05, "YES", "no")))
cat("Models with dAICc > 0 in EVERY kernel:\n")
print(summ[summ$favours_disease_all, c("species", "response", "dAICc_min", "dAICc_max")],
      row.names = FALSE, digits = 3)

# === sanity anchor: alpha = 2, beta = 0.15 must match the committed primary ==
prim <- res[res$alpha == 2 & res$beta == 0.15,
            c("species", "response", "slope_per_SD", "p", "dAICc")]
ref  <- read.csv("output/tables/bioassay_regression.csv")
m <- merge(prim, ref, by = c("species", "response"), suffixes = c(".new", ".ref"))
stopifnot(nrow(m) == nrow(prim))
cat(sprintf(
  "\nSanity anchor (alpha=2, beta=0.15 cell vs bioassay_regression.csv):\n  max|slope| %.2e | max|p| %.2e | max|dAICc| %.2e\n",
  max(abs(m$slope_per_SD.new - m$slope_per_SD.ref)),
  max(abs(m$p.new - m$p.ref)),
  max(abs(m$dAICc.new - m$dAICc.ref))))

cat("\nWrote output/tables/bioassay_regression_sensitivity.csv (",
    nrow(res), " rows )\n", sep = "")
