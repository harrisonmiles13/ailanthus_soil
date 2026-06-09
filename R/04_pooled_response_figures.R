# 04_pooled_response_figures.R
# ---------------------------------------------------------------------------
# EXPLORATORY single-line variant of 03. The bioassay was sown and grown all at
# once, ex situ, so from a seedling's view a soil sample is just a disease
# exposure -- this script drops the per-plot fit AND the month/plot grouping and
# asks the POOLED question: across ALL soil samples, does the time-matched
# ANI_disease predict the response? One pooled fitted line per panel; points are
# shaped by plot (visual only -- plot is NOT in the model, no month shown).
#
#   length, ratio : lm(response ~ z_ani)                       # NO plot/month
#   germination   : glm(cbind(germ,fail) ~ z_ani, quasibinomial)
#                   -- same overdispersed seed-count model as 02/03, just with
#                      the plot block dropped; fitted curve is the pooled logistic.
#   line          : single pooled fit + 95% band over all samples in the panel
#
# READ WITH CARE -- this is deliberately the contrast 02/03 avoid. Most of the
# ANI_disease variance is BETWEEN the 3 plots, and the 3 plots ARE the 3
# inoculation treatments, so the pooled slope is dominated by a between-plot
# (effective n = 3), fully-confounded contrast: it conflates Verticillium with
# anything else differing between those physical plots (soil, drainage, seed
# bank) and with calendar season. The within-plot fits in 02/03 remain the
# defensible inferential version; these are a descriptive complement only.
#
# Outputs:
#   * output/figures/<species>_pooled.png        (4 figures, 3 panels each)
#   * output/tables/pooled_regression_full.txt   (full per-fit summaries, 12 fits)
# Does NOT touch any 02/03 output.
# ---------------------------------------------------------------------------

suppressMessages({
  library(readxl); library(readr); library(ggplot2)
})

SEEDS   <- c(bes = 10, jg = 100, wsr = 40, yar = 40)
SPECIES <- c(bes = "Black-eyed susan", jg = "Johnsongrass",
             wsr = "White snakeroot",  yar = "Yarrow")
RESP <- list(
  total_avg = c("Total seedling length", "length"),
  ratio_avg = c("Root:shoot ratio",      "ratio"),
  germ_perc = c("Germination (proportion)", "germination"))
RESP_ORDER <- unname(vapply(RESP, `[`, character(1), 1))

ani  <- read_csv("output/tables/ani_disease.csv", show_col_types = FALSE)
resp <- read_excel("bioassay_data_v_final.xlsx", "bioassay_primary")

# time-matched ANI_disease per soil sample: may = 0, july = 2mai, sep = 4mai
stopifnot(all(ani$month %in% c("may", "july", "september")))
ani$ani_tm <- 0
ani$ani_tm[ani$month == "july"]      <- ani$ani_disease_2mai[ani$month == "july"]
ani$ani_tm[ani$month == "september"] <- ani$ani_disease_4mai[ani$month == "september"]

respcols <- as.vector(t(outer(names(SPECIES), names(RESP), paste, sep = "_")))
dat <- merge(ani[, c("pooled_sample", "plot", "month", "ani_tm")],
             resp[, c("pooled_sample", respcols)], by = "pooled_sample")
dat$plot  <- factor(dat$plot)
dat$z_ani <- as.numeric(scale(dat$ani_tm))   # per-SD (same scaling family as 02/03)
mu <- mean(dat$ani_tm); sdv <- sd(dat$ani_tm)

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables",  showWarnings = FALSE, recursive = TRUE)

# === Pooled fit, error model matched to 02/03 (NO plot/month term) ===========
# Returns the z_ani slope row and a single pooled fitted line (+95% band) on the
# response scale: lm CI for Gaussian; Wald-on-logit band for the quasibinomial.
pooled_fit <- function(d, sp_code, rcode) {
  gx <- seq(min(d$ani_tm), max(d$ani_tm), length.out = 100)
  gz <- (gx - mu) / sdv
  if (rcode == "germ_perc") {
    n_seed <- SEEDS[[sp_code]]
    d$germ <- round(d$y * n_seed)
    d$fail <- n_seed - d$germ
    m  <- glm(cbind(germ, fail) ~ z_ani, family = quasibinomial, data = d)
    pr <- predict(m, newdata = data.frame(z_ani = gz), type = "link", se.fit = TRUE)
    line <- data.frame(ani_tm = gx, fit = plogis(pr$fit),
                       lwr = plogis(pr$fit - 1.96 * pr$se.fit),
                       upr = plogis(pr$fit + 1.96 * pr$se.fit))
  } else {
    m  <- lm(y ~ z_ani, data = d)
    pr <- as.data.frame(predict(m, newdata = data.frame(z_ani = gz),
                                interval = "confidence"))
    line <- data.frame(ani_tm = gx, fit = pr$fit, lwr = pr$lwr, upr = pr$upr)
  }
  co <- summary(m)$coefficients["z_ani", ]
  list(slope = co[[1]], se = co[[2]], stat = co[[3]], p = co[[4]], line = line)
}

# --- per-species 3-panel figure: ONE pooled line per response ----------------
make_pooled_panel <- function(sp_code, sp_name) {
  long <- lines <- annot <- data.frame()
  for (rcode in names(RESP)) {
    rlabel <- RESP[[rcode]][1]
    col <- paste0(sp_code, "_", rcode)
    d <- dat[!is.na(dat[[col]]), ]; d$y <- d[[col]]

    f  <- pooled_fit(d, sp_code, rcode)
    rf <- factor(rlabel, levels = RESP_ORDER)
    long  <- rbind(long,  data.frame(ani_tm = d$ani_tm, y = d$y,
                                     plot = as.character(d$plot), response = rf))
    lines <- rbind(lines, data.frame(f$line, response = rf))
    annot <- rbind(annot, data.frame(response = rf,
                     label = sprintf("pooled slope %+.3f/SD\np = %.3f  (n = %d)",
                                     f$slope, f$p, nrow(d))))
  }

  plt <- ggplot(long, aes(ani_tm, y)) +
    geom_ribbon(data = lines, aes(ani_tm, ymin = lwr, ymax = upr),
                inherit.aes = FALSE, fill = "grey70", alpha = 0.5) +
    geom_line(data = lines, aes(ani_tm, fit), inherit.aes = FALSE,
              colour = "black", linewidth = 0.9) +
    geom_point(aes(shape = plot), colour = "grey30", size = 2.4, alpha = 0.8) +
    geom_text(data = annot, aes(label = label), x = -Inf, y = Inf,
              hjust = -0.04, vjust = 1.15, size = 3, colour = "grey20",
              inherit.aes = FALSE) +
    facet_wrap(~ response, scales = "free_y") +
    labs(
      caption  = "Pooled slope is mostly the confounded between-plot/treatment contrast -- see 03 for the within-plot test",
      x = "ANI_disease (time-matched AUDPC: may = 0, july = 2mai, sep = 4mai)", y = NULL,
      shape = "Plot (not modelled)"
    ) +
    theme_bw(base_size = 11)

  ggsave(sprintf("output/figures/%s_pooled.png", sp_code), plt,
         width = 11, height = 4.2, dpi = 150)
  cat(sprintf("  %-18s -> output/figures/%s_pooled.png\n", sp_name, sp_code))
}

cat("Writing 4 per-species single-line (pooled) figures, points shaped by plot:\n")
for (sc in names(SPECIES)) make_pooled_panel(sc, SPECIES[[sc]])

# --- full pooled regression outputs (console + text file) -------------------
ovw <- data.frame()
for (sc in names(SPECIES)) for (rc in names(RESP)) {
  col <- paste0(sc, "_", rc); d <- dat[!is.na(dat[[col]]), ]; d$y <- d[[col]]
  if (rc == "germ_perc") {
    n_seed <- SEEDS[[sc]]; d$germ <- round(d$y * n_seed); d$fail <- n_seed - d$germ
    m  <- glm(cbind(germ, fail) ~ z_ani, family = quasibinomial, data = d)
    co <- summary(m)$coefficients["z_ani", ]; fam <- "quasibinomial"; r2 <- NA_real_
  } else {
    s  <- summary(lm(y ~ z_ani, data = d)); co <- s$coefficients["z_ani", ]
    fam <- "gaussian"; r2 <- s$r.squared
  }
  ovw <- rbind(ovw, data.frame(species = SPECIES[[sc]], response = RESP[[rc]][2],
    family = fam, n = nrow(d), slope_per_SD = co[[1]], se = co[[2]],
    statistic = co[[3]], p = co[[4]], R2 = r2))
}

txt <- "output/tables/pooled_regression_full.txt"
sink(txt, split = TRUE)
cat("Pooled regression -- full outputs: each species alone, NO plot/month term.\n")
cat("x = time-matched ANI_disease (may=0, july=2mai, sep=4mai), per SD.\n")
cat("Error model matched to 02: Gaussian lm for length/ratio; quasibinomial glm\n")
cat("(overdispersed seed counts) for germination -- germ slope is log-odds per SD.\n")
cat("CAVEAT: the pooled slope is mostly the confounded between-plot/treatment contrast\n")
cat("        (effective n = 3 plots); the within-plot (defensible) version is R/02 + R/03.\n")
cat("        These pooled fits are a descriptive complement, not the inferential headline.\n\n")
cat("==== Overview (12 fits; slope per SD of ANI_disease) ====\n")
print(ovw, row.names = FALSE, digits = 4)
cat("\n==== Full summaries ====\n")
for (sc in names(SPECIES)) for (rc in names(RESP)) {
  col <- paste0(sc, "_", rc); d <- dat[!is.na(dat[[col]]), ]; d$y <- d[[col]]
  cat(sprintf("\n---------- %s : %s  (n = %d) ----------\n",
              SPECIES[[sc]], RESP[[rc]][1], nrow(d)))
  if (rc == "germ_perc") {
    n_seed <- SEEDS[[sc]]; d$germ <- round(d$y * n_seed); d$fail <- n_seed - d$germ
    m <- glm(cbind(germ, fail) ~ z_ani, family = quasibinomial, data = d)
    print(summary(m))
    ci <- confint.default(m)["z_ani", ]
    cat(sprintf("95%% Wald CI for z_ani slope (log-odds per SD): [%+.4f, %+.4f]\n",
                ci[[1]], ci[[2]]))
  } else {
    m <- lm(y ~ z_ani, data = d); print(summary(m))
    ci <- confint(m)["z_ani", ]
    cat(sprintf("95%% CI for z_ani slope (per SD): [%+.4f, %+.4f]\n", ci[[1]], ci[[2]]))
  }
}
sink()

cat(sprintf("\nDone. Figures: output/figures/<species>_pooled.png | full stats: %s\n", txt))
cat("Caveat: pooled slope mostly = between-plot/treatment (confounded, eff. n = 3). Within-plot test = 02/03.\n")
