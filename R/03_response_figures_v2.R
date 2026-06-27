# 03_response_figures_v2.R
# ---------------------------------------------------------------------------
# v2 restyle of 03_response_figures.R. SAME within-plot fits and SAME numbers
# (Gaussian lm for length/ratio; quasibinomial glm for germination; plot
# partialled out) -- only the presentation changes, per the v2 brief: large
# text, generous space, no titles, colour-blind-safe palette with redundant
# linetype + shape. Reads the v1 outputs, writes to output/figures_v2/; the v1
# script and output/figures/*.png are left untouched.
#
#   * 12 single-panel figures (3 responses x 4 species)  -> <species>_<tag>.png
#   * 4 per-species 3-panel figures (length|ratio|germ)  -> <species>_panel.png
# ---------------------------------------------------------------------------

suppressMessages({
  library(readxl); library(readr); library(ggplot2)
})
source("R/v2_style.R")

SEEDS   <- c(bes = 10, jg = 100, wsr = 40, yar = 40)   # seeds sown per tray (germ denom)
SPECIES <- c(bes = "Black-eyed Susan", jg = "Johnsongrass",
             wsr = "White snakeroot",  yar = "Yarrow")
RESP <- list(                                   # column -> (axis label, file tag)
  total_avg = c("Total seedling length (mm)", "length"),
  ratio_avg = c("Root-to-shoot ratio",        "ratio"),
  germ_perc = c("Germination (proportion)",   "germination"))
RESP_ORDER <- unname(vapply(RESP, `[`, character(1), 1))  # panel order = RESP order

ani  <- read_csv("output/tables/ani_disease.csv", show_col_types = FALSE)
resp <- read_excel("bioassay_data_v_final.xlsx", "bioassay_primary")

# === Time-matched ANI_disease: may = 0, july = 2mai, september = 4mai ========
stopifnot(all(ani$month %in% c("may", "july", "september")))
ani$ani_tm <- 0
ani$ani_tm[ani$month == "july"]      <- ani$ani_disease_2mai[ani$month == "july"]
ani$ani_tm[ani$month == "september"] <- ani$ani_disease_4mai[ani$month == "september"]

respcols <- as.vector(t(outer(names(SPECIES), names(RESP), paste, sep = "_")))
dat <- merge(ani[, c("pooled_sample", "plot", "month", "ani_tm")],
             resp[, c("pooled_sample", respcols)], by = "pooled_sample")
dat$plot  <- factor(dat$plot)
dat$z_ani <- as.numeric(scale(dat$ani_tm))      # common per-SD scaling (matches 02)
mu <- mean(dat$ani_tm); sdv <- sd(dat$ani_tm)   # to map grid x -> z for predict()

# === Within-plot fit, error model matched to the response (identical to v1) ===
fit_within <- function(d, sp_code, rcode) {
  if (rcode == "germ_perc") {
    n_seed <- SEEDS[[sp_code]]
    d$germ <- round(d$y * n_seed)
    d$fail <- n_seed - d$germ
    m       <- glm(cbind(germ, fail) ~ z_ani + plot, family = quasibinomial, data = d)
    predfun <- function(g) predict(m, newdata = g, type = "response")
  } else {
    m       <- lm(y ~ z_ani + plot, data = d)
    predfun <- function(g) predict(m, newdata = g)
  }
  co <- summary(m)$coefficients["z_ani", ]
  grid <- do.call(rbind, lapply(levels(d$plot), function(pl) {
    rng <- range(d$ani_tm[d$plot == pl])
    g <- data.frame(ani_tm = seq(rng[1], rng[2], length.out = 50),
                    plot = factor(pl, levels = levels(d$plot)))
    g$z_ani <- (g$ani_tm - mu) / sdv
    g$fit   <- predfun(g)
    g
  }))
  list(slope = co[["Estimate"]], p = co[["Pr(>|t|)"]], grid = grid)
}

pfmt <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)

# --- 12 single-panel figures -------------------------------------------------
make_fig <- function(sp_code, sp_name, rcode, rlabel, rtag) {
  col <- paste0(sp_code, "_", rcode)
  d <- dat[!is.na(dat[[col]]), ]
  d$y    <- d[[col]]
  d$plot <- droplevels(d$plot)

  f <- fit_within(d, sp_code, rcode)
  d$plot_lab     <- recode_plot(d$plot)
  f$grid$plot_lab <- recode_plot(f$grid$plot)
  annot <- data.frame(label = sprintf("slope %+.3f / SD\n%s  (n = %d)",
                                      f$slope, pfmt(f$p), nrow(d)))

  plt <- ggplot(d, aes(ani_tm, y, colour = plot_lab)) +
    geom_line(data = f$grid, aes(ani_tm, fit, colour = plot_lab, linetype = plot_lab),
              linewidth = 1.2) +
    geom_point(aes(shape = month), size = 3.6, alpha = 0.9, stroke = 0.7) +
    annot_layer(annot, size = 5) +
    scale_plot_colour() + scale_plot_linetype() + scale_month_shape() +
    labs(x = LAB_ANI, y = rlabel) +
    guides(colour = guide_legend(nrow = 1, title.position = "top",
                                 override.aes = list(linewidth = 1.2)),
           linetype = guide_legend(nrow = 1, title.position = "top"),
           shape = guide_legend(nrow = 1, title.position = "top")) +
    theme_v2(base_size = 18) +
    theme(legend.title = element_text(size = rel(0.88), face = "bold", hjust = 0.5))

  fn <- fig_v2(sprintf("%s_%s.png", sp_code, rtag))
  ggsave(fn, plt, width = 9, height = 6.6, dpi = 200)
  cat(sprintf("  %-26s slope %+.3f /SD  p = %.3f  (n = %d)  -> %s\n",
              paste(sp_code, rtag), f$slope, f$p, nrow(d), fn))
}

# --- per-species 3-panel manuscript figure -----------------------------------
make_species_panel <- function(sp_code, sp_name) {
  long <- grids <- annot <- data.frame()
  for (rcode in names(RESP)) {
    rlabel <- RESP[[rcode]][1]
    col <- paste0(sp_code, "_", rcode)
    d <- dat[!is.na(dat[[col]]), ]
    d$y    <- d[[col]]
    d$plot <- droplevels(d$plot)

    f <- fit_within(d, sp_code, rcode)

    rf    <- factor(rlabel, levels = RESP_ORDER)
    long  <- rbind(long,  data.frame(ani_tm = d$ani_tm, y = d$y,
                                     plot_lab = recode_plot(d$plot), month = d$month,
                                     response = rf))
    grids <- rbind(grids, data.frame(ani_tm = f$grid$ani_tm, fit = f$grid$fit,
                                     plot_lab = recode_plot(f$grid$plot), response = rf))
    annot <- rbind(annot, data.frame(response = rf,
                     label = sprintf("slope %+.3f / SD\n%s  (n = %d)",
                                     f$slope, pfmt(f$p), nrow(d))))
  }

  plt <- ggplot(long, aes(ani_tm, y, colour = plot_lab)) +
    geom_line(data = grids, aes(ani_tm, fit, colour = plot_lab, linetype = plot_lab),
              linewidth = 1.1) +
    geom_point(aes(shape = month), size = 3.0, alpha = 0.9, stroke = 0.6) +
    annot_layer(annot, size = 4.2) +
    facet_wrap(~ response, scales = "free_y") +
    scale_plot_colour() + scale_plot_linetype() + scale_month_shape() +
    labs(x = LAB_ANI, y = NULL) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)),
           linetype = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
    theme_v2(base_size = 16)

  fn <- fig_v2(sprintf("%s_panel.png", sp_code))
  ggsave(fn, plt, width = 13, height = 5.8, dpi = 200)
  cat(sprintf("  %-18s -> %s\n", sp_name, fn))
}

cat("Writing 12 single-panel v2 figures (3 responses x 4 species):\n")
for (sc in names(SPECIES)) for (rc in names(RESP))
  make_fig(sc, SPECIES[[sc]], rc, RESP[[rc]][1], RESP[[rc]][2])

cat("\nWriting 4 per-species 3-panel v2 figures (length | ratio | germination):\n")
for (sc in names(SPECIES)) make_species_panel(sc, SPECIES[[sc]])

cat("\nDone. v2 figures in output/figures_v2/.\n")
