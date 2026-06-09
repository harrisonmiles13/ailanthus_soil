# 02_bioassay_regression.R
# ---------------------------------------------------------------------------
# Per-species regression of bioassay response on the TIME-RESOLVED Ailanthus
# disease-neighbourhood index ANI_disease (output/tables/ani_disease.csv, the
# Gomez-Aparicio & Canham 2008 neighbourhood model weighted by tree AUDPC).
#
# One question, asked separately for each of the 4 species:
#   does the disease-neighbourhood load a soil-sample point has accumulated
#   predict the response (total length, root:shoot ratio, germination) measured
#   in soil from that point at that time?
#
# Two time axes, correctly matched:
#
#   * AUDPC timepoint  = accumulated disease ON THE TREES, columns 1mai..4mai
#     ("months after inoculation"), monotonically increasing.
#   * sampling month   = when soil was collected (may / july / september).
#
# Soil collected later in the season has sat under a larger accumulated tree
# disease load. Samplings are 2 months apart, AUDPC readings 1 month apart, so
# the only even mapping that ends on the maximal measured timepoint is
# {0, 2mai, 4mai}:
#
#       may       = pre-inoculation     -> ANI_disease 0  (no disease yet)
#       july      = 2 months post-inoc  -> ani_disease_2mai
#       september = 4 months post-inoc  -> ani_disease_4mai
#
# The predictor genuinely varies across the 3 visits to a point (0 in may,
# moderate in july, highest in september). That time variation IS the disease
# signal, so the 45 point x month rows are NOT collapsed: each row pairs a
# month's response with the disease-neighbourhood load accumulated by that month.
#
# Design facts that fix the form of the regression:
#
#  * Each species is grown in its own cells of soil (same source) and is
#    analysed entirely on its own. No cross-species step of any kind: no
#    pooling, no FDR, no comparison.
#
#  * `plot` is kept in every regression purely to hold the comparison WITHIN
#    plots. Disease load is mostly between-plot (the 3 plots are the 3
#    confounded inoculation treatments); a bare response~ANI slope would be that
#    pseudoreplicated between-plot contrast. With `plot` partialled out, the
#    slope is the confound-free within-plot effect -- the analogue of the
#    neighbourhood coefficient c in the published model. Plot's own coefficients
#    are nuisance and are never interpreted or reported.
#
#  * Error model matched to each response, following the published analyses
#    (Gomez-Aparicio & Canham 2008, Appendix S2: contingency/binomial for
#    survival; ANOVA/normal for growth & biomass):
#
#      - germination : binomial counts. Each pooled_sample x month is ONE tray
#        of N sown seeds (N = 10 bes / 100 jg / 40 wsr / 40 yar), so the per-row
#        response is cbind(germinated, ungerminated) with a fixed, known
#        denominator. The trays are already the binomial units, so no
#        month-summing is needed -- this does NOT reinflate any denominator or
#        re-import month pseudoreplication. These seed counts are strongly
#        OVERDISPERSED (Pearson c-hat ~ 3-13, reported per row), as germination
#        of natural seed lots usually is, so a pure binomial would understate
#        the SEs by ~sqrt(c-hat). Inference therefore uses QUASIBINOMIAL (SEs
#        scaled by sqrt(c-hat)); point estimates of the slope are unchanged.
#      - length, root:shoot ratio : Gaussian lm (continuous growth responses).
#
#  * Model selection per species/response follows the neighbourhood-analysis
#    convention: the neighbourhood model (response ~ z_ani + plot) is compared
#    by a small-sample information criterion against the no-neighbourhood null
#    (response ~ plot). dAICc = IC(null) - IC(full); dAICc > 0 favours a real
#    disease-influence effect. For Gaussian responses IC = AICc; for the
#    overdispersed germination counts IC = QAICc (quasi-AICc, the same AICc with
#    the deviance divided by the model's c-hat and one extra parameter counted
#    for estimating c-hat -- Burnham & Anderson's overdispersion correction). A
#    familiar slope p-value is reported alongside, but dAICc is the framework's
#    selection statistic.
#
# Two honest caveats, carried in the interpretation rather than "fixed" here:
#
#  * The 45 rows are repeated measures of 15 points (3 months each), so the
#    fixed-effect SEs are mildly optimistic. The textbook fix `(1|point)` is the
#    mixed model that is deliberately NOT used (n is too small to spend the df).
#  * Time-matched influence is collinear with calendar season, so the slope
#    conflates accumulated disease with any common seasonal trend. The
#    no_fungus_control plot (Ailanthus, no Verticillium -> ~0 disease all
#    season) is the natural seasonal reference, but using it would be a
#    prohibited between-plot comparison.
#
# Scope: the three Ailanthus plots only. No neg_control, no ailanthus_presence.
# ---------------------------------------------------------------------------

suppressMessages({
  library(readxl)
  library(readr)
})

SEEDS   <- c(bes = 10, jg = 100, wsr = 40, yar = 40)
SPECIES <- names(SEEDS)

ani  <- read_csv("output/tables/ani_disease.csv", show_col_types = FALSE)
resp <- read_excel("bioassay_data_v_final.xlsx", "bioassay_primary")

# === Time-matched disease-neighbourhood index per sampling month ============
# may = pre-inoculation (0); july = 2mai load; september = 4mai load.
stopifnot(all(ani$month %in% c("may", "july", "september")))
ani$ani_tm <- 0                                            # may: no disease yet
is_jul <- ani$month == "july"
is_sep <- ani$month == "september"
ani$ani_tm[is_jul] <- ani$ani_disease_2mai[is_jul]
ani$ani_tm[is_sep] <- ani$ani_disease_4mai[is_sep]

# === Assemble the 45 point x month rows =====================================
keep_resp <- c("pooled_sample", paste0(SPECIES, "_total_avg"),
               paste0(SPECIES, "_ratio_avg"), paste0(SPECIES, "_germ_perc"))
dat <- merge(
  ani[, c("pooled_sample", "plot", "month", "distance", "bearing", "ani_tm")],
  resp[, keep_resp], by = "pooled_sample")
dat$plot  <- factor(dat$plot)
dat$z_ani <- as.numeric(scale(dat$ani_tm))   # per-SD units; +plot makes it within-plot

stopifnot(nrow(dat) == 45)
cat(sprintf("Regression unit: %d point x month rows (15 points x 3 months) | predictor = time-matched within-plot ANI_disease (per SD)\n\n",
            nrow(dat)))

# === Information criteria (small-sample corrected) ==========================
# AICc for the Gaussian models; QAICc for the overdispersed binomial counts.
# Both take the deviance/logLik from a likelihood fit; QAICc divides by c-hat
# (overdispersion from the GLOBAL model) and counts one extra parameter for it.
aicc <- function(m) {
  ll <- logLik(m)
  k  <- attr(ll, "df")          # coefficients (+ sigma for lm)
  n  <- nobs(m)
  as.numeric(-2 * ll + 2 * k + (2 * k * (k + 1)) / (n - k - 1))
}
qaicc <- function(m, chat, n) {
  ll <- logLik(m)
  k  <- attr(ll, "df") + 1      # +1 for estimating c-hat (Burnham & Anderson)
  as.numeric(-2 * as.numeric(ll) / chat + 2 * k + (2 * k * (k + 1)) / (n - k - 1))
}

# === Per-species fit: response ~ time-matched ANI (within plot) =============
# Gaussian for length/ratio; binomial counts (quasibinomial inference) for
# germination. The z_ani slope is the within-plot neighbourhood effect (plot
# partialled out). dAICc compares the neighbourhood model against the plot-only
# null (AICc for Gaussian, QAICc for the overdispersed germination counts).
fit_one <- function(y_full, label, species_name, family) {
  if (family == "gaussian") {
    df <- data.frame(y = y_full, z_ani = dat$z_ani, plot = dat$plot)
    df <- df[!is.na(df$y), ]
    m0 <- lm(y ~ plot,          data = df)
    m1 <- lm(y ~ z_ani + plot,  data = df)
    co <- summary(m1)$coefficients["z_ani", ]
    slope <- co[["Estimate"]]; se <- co[["Std. Error"]]
    stat  <- co[["t value"]];  p  <- co[["Pr(>|t|)"]]
    disp  <- NA_real_
    dic   <- aicc(m0) - aicc(m1)
    unit  <- "resp units / SD"
  } else {                                   # germination: binomial counts
    n_seed <- SEEDS[[species_name]]
    germ   <- round(y_full * n_seed)
    df <- data.frame(germ = germ, fail = n_seed - germ,
                     z_ani = dat$z_ani, plot = dat$plot)
    df <- df[!is.na(df$germ), ]
    n  <- nrow(df)
    # likelihood (binomial) fits supply deviance + logLik for QAICc;
    # c-hat from the global model corrects the overdispersion.
    b0 <- glm(cbind(germ, fail) ~ plot,         family = binomial, data = df)
    b1 <- glm(cbind(germ, fail) ~ z_ani + plot, family = binomial, data = df)
    disp <- sum(residuals(b1, "pearson")^2) / df.residual(b1)   # Pearson c-hat
    chat <- max(disp, 1)                       # never reward underdispersion
    dic  <- qaicc(b0, chat, n) - qaicc(b1, chat, n)
    # quasibinomial supplies the dispersion-corrected slope inference
    q1 <- glm(cbind(germ, fail) ~ z_ani + plot, family = quasibinomial, data = df)
    co <- summary(q1)$coefficients["z_ani", ]
    slope <- co[["Estimate"]]; se <- co[["Std. Error"]]
    stat  <- co[["t value"]];  p  <- co[["Pr(>|t|)"]]
    family <- "quasibinomial"
    unit   <- "log-odds / SD"
  }
  data.frame(species = species_name, response = label, family = family,
             n = nrow(df), slope_per_SD = slope, se = se, statistic = stat,
             p = p, dAICc = dic, dispersion = disp, effect_unit = unit)
}

res <- data.frame()
for (s in SPECIES) {
  res <- rbind(res,
    fit_one(dat[[paste0(s, "_total_avg")]], "length",      s, "gaussian"),
    fit_one(dat[[paste0(s, "_ratio_avg")]], "ratio",       s, "gaussian"),
    fit_one(dat[[paste0(s, "_germ_perc")]], "germination", s, "binomial"))
}

for (r in c("length", "ratio", "germination")) {
  cat(sprintf("\n=== %s ~ time-matched ANI_disease (per SD), each species alone, within plot ===\n",
              tools::toTitleCase(r)))
  print(res[res$response == r, ], row.names = FALSE, digits = 3)
}

write.csv(res, "output/tables/bioassay_regression.csv", row.names = FALSE)
cat("\nWrote output/tables/bioassay_regression.csv\n")
cat("Predictor = time-matched ANI_disease (may=0, july=2mai, sep=4mai) | gaussian lm (length,ratio) + quasibinomial glm (germ, overdispersed) | plot = within-plot nuisance | dAICc = (Q)AICc(null) - (Q)AICc(neighbourhood) | no FDR, no cross-species step\n")
