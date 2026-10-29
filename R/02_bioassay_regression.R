suppressMessages({
  library(readxl)
  library(readr)
})

SEEDS   <- c(bes = 10, jg = 100, wsr = 40, yar = 40)
SPECIES <- names(SEEDS)

ani  <- read_csv("output/tables/ani_disease.csv", show_col_types = FALSE)
resp <- read_excel("bioassay_data_v_final.xlsx", "bioassay_primary")

stopifnot(all(ani$month %in% c("may", "july", "september")))
ani$ani_tm <- 0
is_jul <- ani$month == "july"
is_sep <- ani$month == "september"
ani$ani_tm[is_jul] <- ani$ani_disease_2mai[is_jul]
ani$ani_tm[is_sep] <- ani$ani_disease_4mai[is_sep]

keep_resp <- c("pooled_sample", paste0(SPECIES, "_total_avg"),
               paste0(SPECIES, "_ratio_avg"), paste0(SPECIES, "_germ_perc"))
dat <- merge(
  ani[, c("pooled_sample", "plot", "month", "distance", "bearing", "ani_tm")],
  resp[, keep_resp], by = "pooled_sample")
dat$plot  <- factor(dat$plot)
dat$z_ani <- as.numeric(scale(dat$ani_tm))

stopifnot(nrow(dat) == 45)
cat(sprintf("Regression unit: %d point x month rows (15 points x 3 months) | predictor = time-matched within-plot ANI_disease (per SD)\n\n",
            nrow(dat)))

aicc <- function(m) {
  ll <- logLik(m)
  k  <- attr(ll, "df")
  n  <- nobs(m)
  as.numeric(-2 * ll + 2 * k + (2 * k * (k + 1)) / (n - k - 1))
}
qaicc <- function(m, chat, n) {
  ll <- logLik(m)
  k  <- attr(ll, "df") + 1
  as.numeric(-2 * as.numeric(ll) / chat + 2 * k + (2 * k * (k + 1)) / (n - k - 1))
}

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
  } else {
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
