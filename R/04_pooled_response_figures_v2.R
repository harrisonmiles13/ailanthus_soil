# 04_pooled_response_figures_v2.R
# ---------------------------------------------------------------------------
# v2 restyle of 04_pooled_response_figures.R (the EXPLORATORY pooled variant).
# SAME pooled fits and SAME numbers as v1 (lm for length/ratio; quasibinomial
# glm for germination; NO plot/month term; single pooled line + 95% band) --
# only the presentation changes per the v2 brief (large text, generous space,
# no titles, colour-blind-safe palette). The confounding CAVEAT is kept, as a
# caption (not a title), because it is essential to reading these fits.
#
# Writes figures to output/figures_v2/; the full-stats text table is the v1
# script's job, so this v2 script does NOT rewrite output/tables/.
#   * output/figures_v2/<species>_pooled.png   (4 figures, 3 panels each)
# ---------------------------------------------------------------------------

suppressMessages({
  library(readxl); library(readr); library(ggplot2)
})
source("R/v2_style.R")

SEEDS   <- c(bes = 10, jg = 100, wsr = 40, yar = 40)
SPECIES <- c(bes = "Black-eyed Susan", jg = "Johnsongrass",
             wsr = "White snakeroot",  yar = "Yarrow")
RESP <- list(
  total_avg = c("Total seedling length (mm)", "length"),
  ratio_avg = c("Root-to-shoot ratio",        "ratio"),
  germ_perc = c("Germination (proportion)",   "germination"))
RESP_ORDER <- unname(vapply(RESP, `[`, character(1), 1))

ani  <- read_csv("output/tables/ani_disease.csv", show_col_types = FALSE)
resp <- read_excel("bioassay_data_v_final.xlsx", "bioassay_primary")

stopifnot(all(ani$month %in% c("may", "july", "september")))
ani$ani_tm <- 0
ani$ani_tm[ani$month == "july"]      <- ani$ani_disease_2mai[ani$month == "july"]
ani$ani_tm[ani$month == "september"] <- ani$ani_disease_4mai[ani$month == "september"]

respcols <- as.vector(t(outer(names(SPECIES), names(RESP), paste, sep = "_")))
dat <- merge(ani[, c("pooled_sample", "plot", "month", "ani_tm")],
             resp[, c("pooled_sample", respcols)], by = "pooled_sample")
dat$plot  <- factor(dat$plot)
dat$z_ani <- as.numeric(scale(dat$ani_tm))
mu <- mean(dat$ani_tm); sdv <- sd(dat$ani_tm)

# === Pooled fit, error model matched to v1 (NO plot/month term) ==============
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

pfmt <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)

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
                                     plot_lab = recode_plot(d$plot), response = rf))
    lines <- rbind(lines, data.frame(f$line, response = rf))
    annot <- rbind(annot, data.frame(response = rf,
                     label = sprintf("pooled slope %+.3f / SD\n%s  (n = %d)",
                                     f$slope, pfmt(f$p), nrow(d))))
  }

  plt <- ggplot(long, aes(ani_tm, y)) +
    geom_ribbon(data = lines, aes(ani_tm, ymin = lwr, ymax = upr),
                inherit.aes = FALSE, fill = "grey75", alpha = 0.55) +
    geom_line(data = lines, aes(ani_tm, fit), inherit.aes = FALSE,
              colour = "grey15", linewidth = 1.2) +
    geom_point(aes(shape = plot_lab, colour = plot_lab), size = 3.0,
               alpha = 0.9, stroke = 0.6) +
    annot_layer(annot, size = 4.2) +
    facet_wrap(~ response, scales = "free_y") +
    scale_colour_manual(values = PLOT_COLOURS, name = "Plot (not modelled)") +
    scale_shape_manual(values = c("No-fungus control" = 16,
                                  "Attenuated strain" = 17,
                                  "Virulent strain"   = 15),
                       name = "Plot (not modelled)") +
    labs(x = LAB_ANI, y = NULL,
         caption = paste("Pooled slope is mostly the confounded",
                         "between-plot / treatment contrast (effective n = 3);",
                         "see the within-plot test in Figure 5.")) +
    guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
    theme_v2(base_size = 16) +
    theme(plot.caption = element_text(size = rel(0.72), colour = "grey35",
                                      hjust = 0, margin = margin(t = 8)))

  fn <- fig_v2(sprintf("%s_pooled.png", sp_code))
  ggsave(fn, plt, width = 13, height = 6.0, dpi = 200)
  cat(sprintf("  %-18s -> %s\n", sp_name, fn))
}

cat("Writing 4 per-species pooled v2 figures (points coloured/shaped by plot):\n")
for (sc in names(SPECIES)) make_pooled_panel(sc, SPECIES[[sc]])

cat("\nDone. v2 pooled figures in output/figures_v2/.\n")
cat("Caveat unchanged: pooled slope mostly = confounded between-plot/treatment (eff. n = 3).\n")
