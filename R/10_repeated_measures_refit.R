suppressMessages({
  library(readxl); library(readr); library(dplyr)
  library(lme4); library(lmerTest); library(glmmTMB)
})

OUT <- "output/r2"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

DATA_FILE  <- "bioassay_data_v_final.xlsx"
DIST_FLOOR <- 0.5
SEEDS      <- c(bes = 10, jg = 100, wsr = 40, yar = 40)
SPECIES    <- names(SEEDS)
RESPONSES  <- c(length = "total_avg", ratio = "ratio_avg", germination = "germ_perc")
ALPHA_GRID <- c(0, 1, 2)
BETA_GRID  <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.50)
TREATMENT_TO_PLOT <- c(vnaa140_2019    = "attenuated",
                       vnaa140_2023    = "virulent",
                       vnaa140_control = "no_fungus_control")

ani  <- read_csv("output/tables/ani_disease.csv", show_col_types = FALSE)
resp <- read_excel(DATA_FILE, "bioassay_primary")

ani$ani_tm <- 0
ani$ani_tm[ani$month == "july"]      <- ani$ani_disease_2mai[ani$month == "july"]
ani$ani_tm[ani$month == "september"] <- ani$ani_disease_4mai[ani$month == "september"]

respcols <- as.vector(outer(SPECIES, RESPONSES, paste, sep = "_"))
dat <- merge(ani[, c("pooled_sample", "plot", "month", "distance", "bearing", "ani_tm")],
             resp[, c("pooled_sample", respcols)], by = "pooled_sample")

station_id <- function(plot, distance, bearing)
  paste(plot, ifelse(distance == 0, "center", paste0("b", bearing)), sep = "_")
dat$station <- factor(station_id(dat$plot, dat$distance, dat$bearing))
dat$plot    <- factor(dat$plot)
dat$month   <- factor(dat$month, levels = c("may", "july", "september"))
dat$z_ani   <- as.numeric(scale(dat$ani_tm))
dat$obs     <- factor(seq_len(nrow(dat)))

stopifnot(nrow(dat) == 45, nlevels(dat$station) == 15,
          all(table(dat$station) == 3))

r2_of <- function(f) summary(lm(f, data = dat))$r.squared
confound <- data.frame(
  model = c("z_ani ~ plot", "z_ani ~ month", "z_ani ~ plot + month",
            "z_ani ~ plot * month"),
  R2    = c(r2_of(z_ani ~ plot), r2_of(z_ani ~ month),
            r2_of(z_ani ~ plot + month), r2_of(z_ani ~ plot * month)))
confound$residual_var_share <- 1 - confound$R2
write.csv(confound, file.path(OUT, "index_confounding.csv"), row.names = FALSE)

cell_sd <- dat %>% group_by(plot, month) %>%
  summarise(n = n(), ani_mean = mean(ani_tm), ani_sd = sd(ani_tm), .groups = "drop")
write.csv(cell_sd, file.path(OUT, "index_by_plot_month.csv"), row.names = FALSE)

aicc_ll <- function(ll, n) {
  k <- attr(ll, "df")
  as.numeric(-2 * ll + 2 * k + (2 * k * (k + 1)) / (n - k - 1))
}

capture <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr,
    warning = function(c) { w <<- c(w, trimws(conditionMessage(c))); invokeRestart("muffleWarning") },
    message = function(c) { w <<- c(w, trimws(conditionMessage(c))); invokeRestart("muffleMessage") })
  list(value = val, notes = paste(unique(w), collapse = " | "))
}

fixed_rhs <- function(spec) if (spec == "B") "plot * month" else "plot + month"

fit_gauss <- function(df, spec) {
  rhs <- fixed_rhs(spec)
  f0  <- as.formula(paste("y ~", rhs, "+ (1 | station)"))
  f1  <- as.formula(paste("y ~", rhs, "+ z_ani + (1 | station)"))
  m0  <- capture(lmer(f0, data = df, REML = FALSE))
  m1  <- capture(lmer(f1, data = df, REML = FALSE))
  mr  <- capture(lmerTest::lmer(f1, data = df, REML = TRUE))
  co  <- summary(mr$value, ddf = "Satterthwaite")$coefficients["z_ani", ]
  vc  <- as.data.frame(VarCorr(mr$value))
  n   <- nrow(df)
  data.frame(
    family = "gaussian LMM", n = n,
    slope_per_SD = co[["Estimate"]], se = co[["Std. Error"]],
    df = co[["df"]], statistic = co[["t value"]], p = co[["Pr(>|t|)"]],
    dAICc = aicc_ll(logLik(m0$value), n) - aicc_ll(logLik(m1$value), n),
    station_var = vc$vcov[vc$grp == "station"], resid_var = vc$vcov[vc$grp == "Residual"],
    singular = isSingular(mr$value) || isSingular(m1$value),
    notes = paste(c(m0$notes, m1$notes, mr$notes)[nzchar(c(m0$notes, m1$notes, mr$notes))],
                  collapse = " | "))
}

fit_betabin <- function(df, spec) {
  rhs <- fixed_rhs(spec)
  f0  <- as.formula(paste("cbind(germ, fail) ~", rhs, "+ (1 | station)"))
  f1  <- as.formula(paste("cbind(germ, fail) ~", rhs, "+ z_ani + (1 | station)"))
  m0  <- capture(glmmTMB(f0, family = betabinomial, data = df))
  m1  <- capture(glmmTMB(f1, family = betabinomial, data = df))
  co  <- summary(m1$value)$coefficients$cond["z_ani", ]
  sv  <- VarCorr(m1$value)$cond$station[1, 1]
  n   <- nrow(df)
  conv_ok <- isTRUE(m1$value$sdr$pdHess) && isTRUE(m0$value$sdr$pdHess)
  data.frame(
    family = "beta-binomial GLMM", n = n,
    slope_per_SD = co[["Estimate"]], se = co[["Std. Error"]], df = NA_real_,
    statistic = co[["z value"]], p = co[["Pr(>|z|)"]],
    dAICc = aicc_ll(logLik(m0$value), n) - aicc_ll(logLik(m1$value), n),
    station_var = sv, resid_var = sigma(m1$value),
    singular = sv < 1e-6 || !conv_ok,
    notes = paste(c(m0$notes, m1$notes, if (!conv_ok) "non-PD Hessian")
                  [nzchar(c(m0$notes, m1$notes, if (!conv_ok) "non-PD Hessian"))],
                  collapse = " | "))
}

fit_olre <- function(df, spec) {
  rhs <- fixed_rhs(spec)
  f0  <- as.formula(paste("cbind(germ, fail) ~", rhs, "+ (1 | station) + (1 | obs)"))
  f1  <- as.formula(paste("cbind(germ, fail) ~", rhs, "+ z_ani + (1 | station) + (1 | obs)"))
  ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  m0  <- capture(glmer(f0, family = binomial, data = df, control = ctl))
  m1  <- capture(glmer(f1, family = binomial, data = df, control = ctl))
  co  <- summary(m1$value)$coefficients["z_ani", ]
  vc  <- as.data.frame(VarCorr(m1$value))
  n   <- nrow(df)
  data.frame(
    family = "binomial GLMM + OLRE", n = n,
    slope_per_SD = co[["Estimate"]], se = co[["Std. Error"]], df = NA_real_,
    statistic = co[["z value"]], p = co[["Pr(>|z|)"]],
    dAICc = aicc_ll(logLik(m0$value), n) - aicc_ll(logLik(m1$value), n),
    station_var = vc$vcov[vc$grp == "station"], resid_var = vc$vcov[vc$grp == "obs"],
    singular = isSingular(m1$value),
    notes = paste(c(m0$notes, m1$notes)[nzchar(c(m0$notes, m1$notes))], collapse = " | "))
}

fit_all <- function(d, spec = "A", germ_fun = fit_betabin) {
  out <- list()
  for (s in SPECIES) for (r in names(RESPONSES)) {
    y <- d[[paste0(s, "_", RESPONSES[[r]])]]
    df <- d[!is.na(y), ]; y <- y[!is.na(y)]
    if (r == "germination") {
      df$germ <- round(y * SEEDS[[s]]); df$fail <- SEEDS[[s]] - df$germ
      res <- germ_fun(df, spec)
    } else {
      df$y <- y
      res <- fit_gauss(df, spec)
    }
    out[[length(out) + 1]] <- cbind(species = s, response = r, spec = spec, res)
  }
  do.call(rbind, out)
}

prim  <- fit_all(dat, "A")
specB <- fit_all(dat, "B")
specC <- fit_all(dat, "A", fit_olre)
specC <- specC[specC$response == "germination", ]
specC$spec <- "C"

write.csv(prim, file.path(OUT, "regression_r2_primary.csv"), row.names = FALSE)
write.csv(rbind(specB, specC), file.path(OUT, "regression_r2_checks.csv"), row.names = FALSE)

r1 <- read.csv("output/tables/bioassay_regression.csv")
cmp <- merge(
  r1[, c("species", "response", "n", "slope_per_SD", "se", "p", "dAICc")],
  prim[, c("species", "response", "slope_per_SD", "se", "p", "dAICc", "singular")],
  by = c("species", "response"), suffixes = c("_R1", "_R2"))
cmpB <- specB[, c("species", "response", "slope_per_SD", "p", "dAICc")]
names(cmpB)[3:5] <- paste0(names(cmpB)[3:5], "_B")
cmp <- merge(cmp, cmpB, by = c("species", "response"))
cmp$sig_R1 <- cmp$p_R1 < 0.05
cmp$sig_R2 <- cmp$p_R2 < 0.05
cmp$sig_B  <- cmp$p_B  < 0.05
cmp$favours_disease_R1 <- cmp$dAICc_R1 > 0
cmp$favours_disease_R2 <- cmp$dAICc_R2 > 0
ord <- order(match(cmp$species, SPECIES), match(cmp$response, names(RESPONSES)))
cmp <- cmp[ord, ]
write.csv(cmp, file.path(OUT, "comparison_R1_vs_R2.csv"), row.names = FALSE)

polar_to_xy <- function(distance, bearing) {
  bearing  <- ifelse(is.na(bearing),  0, bearing)
  distance <- ifelse(is.na(distance), 0, distance)
  rad <- bearing * pi / 180
  list(x = distance * sin(rad), y = distance * cos(rad))
}
trees <- read_excel(DATA_FILE, "tree_data") %>%
  mutate(tx = polar_to_xy(Distance_m, Bearing)$x, ty = polar_to_xy(Distance_m, Bearing)$y)
pos <- dat %>% select(pooled_sample, plot, month, distance, bearing) %>%
  mutate(plot = as.character(plot),
         sx = polar_to_xy(distance, bearing)$x, sy = polar_to_xy(distance, bearing)$y)
ani_at <- function(plot, sx, sy, col, a, b) {
  pt <- trees[trees$Plot == plot, ]
  d  <- pmax(sqrt((pt$tx - sx)^2 + (pt$ty - sy)^2), DIST_FLOOR)
  sum(pt[[col]] * pt$dbh_cm^a * exp(-b * d))
}

sweep <- list()
for (a in ALPHA_GRID) for (b in BETA_GRID) {
  i2 <- mapply(function(p, x, y) ani_at(p, x, y, "audpc_adj_2mai", a, b), pos$plot, pos$sx, pos$sy)
  i4 <- mapply(function(p, x, y) ani_at(p, x, y, "audpc_adj_4mai", a, b), pos$plot, pos$sx, pos$sy)
  tm <- ifelse(pos$month == "may", 0, ifelse(pos$month == "july", i2, i4))
  dk <- dat; dk$z_ani <- as.numeric(scale(tm))
  sweep[[length(sweep) + 1]] <- cbind(alpha = a, beta = b, fit_all(dk, "A"))
}
sweep <- do.call(rbind, sweep)
write.csv(sweep, file.path(OUT, "sensitivity_r2.csv"), row.names = FALSE)

anc <- merge(sweep[sweep$alpha == 2 & sweep$beta == 0.15, c("species", "response", "slope_per_SD", "p")],
             prim[, c("species", "response", "slope_per_SD", "p")], by = c("species", "response"))
stopifnot(max(abs(anc$slope_per_SD.x - anc$slope_per_SD.y) /
              pmax(abs(anc$slope_per_SD.y), 1e-8)) < 1e-3)

sweep_summ <- sweep %>% group_by(species, response) %>%
  summarise(p_min = min(p), p_max = max(p),
            dAICc_min = min(dAICc), dAICc_max = max(dAICc),
            n_sig = sum(p < 0.05), n_favour = sum(dAICc > 0),
            n_singular = sum(singular), .groups = "drop") %>%
  arrange(match(species, SPECIES), match(response, names(RESPONSES)))
write.csv(sweep_summ, file.path(OUT, "sensitivity_r2_summary.csv"), row.names = FALSE)

show <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 3)
cat("\n=== Index confounding (share of z_ani variance explained) ===\n"); show(confound)
cat("\n=== Index by plot x month ===\n"); show(cell_sd)
cat("\n=== PRIMARY (A): y ~ plot + month + z_ani + (1|station) ===\n")
show(prim[, c("species", "response", "family", "n", "slope_per_SD", "se", "statistic",
              "p", "dAICc", "station_var", "singular")])
cat("\n=== CHECK B: plot x month fixed ===\n")
show(specB[, c("species", "response", "slope_per_SD", "se", "p", "dAICc", "singular")])
cat("\n=== CHECK C: germination, binomial + OLRE ===\n")
show(specC[, c("species", "slope_per_SD", "se", "p", "dAICc", "singular")])
cat("\n=== R1 vs R2 ===\n")
show(cmp[, c("species", "response", "slope_per_SD_R1", "p_R1", "dAICc_R1",
             "slope_per_SD_R2", "p_R2", "dAICc_R2", "p_B")])
cat("\n=== Kernel sweep (primary spec, 18 kernels) ===\n"); show(sweep_summ)
cat(sprintf("\nSweep: %d fits | min p = %.4f | fits with p < 0.05: %d\n",
            nrow(sweep), min(sweep$p), sum(sweep$p < 0.05)))
notes <- unique(c(prim$notes, specB$notes, specC$notes))
notes <- notes[nzchar(notes)]
if (length(notes)) { cat("\nFit notes:\n"); cat(paste0("  - ", notes), sep = "\n") }

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("\nWrote", OUT, "\n")
