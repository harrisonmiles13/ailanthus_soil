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

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

fit_within <- function(d, sp_code, rcode) {
  if (rcode == "germ_perc") {
    n_seed <- SEEDS[[sp_code]]
    d$germ <- round(d$y * n_seed)
    d$fail <- n_seed - d$germ
    m       <- glm(cbind(germ, fail) ~ z_ani + plot, family = quasibinomial, data = d)
    predfun <- function(g) predict(m, newdata = g, type = "response")
    fam     <- "quasibinomial glm; slope = log-odds / SD"
  } else {
    m       <- lm(y ~ z_ani + plot, data = d)
    predfun <- function(g) predict(m, newdata = g)
    fam     <- "Gaussian lm; slope = resp units / SD"
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
  list(slope = co[["Estimate"]], p = co[["Pr(>|t|)"]], grid = grid, family = fam)
}

make_fig <- function(sp_code, sp_name, rcode, rlabel, rtag) {
  col <- paste0(sp_code, "_", rcode)
  d <- dat[!is.na(dat[[col]]), ]
  d$y    <- d[[col]]
  d$plot <- droplevels(d$plot)

  f <- fit_within(d, sp_code, rcode)

  plt <- ggplot(d, aes(ani_tm, y, colour = plot)) +
    geom_line(data = f$grid, aes(ani_tm, fit, colour = plot), linewidth = 0.8) +
    geom_point(aes(shape = month), size = 2.4, alpha = 0.9) +
    labs(
      x = "ANI_disease (time-matched: may = 0, july = 2mai, sep = 4mai)",
      y = rlabel, colour = "Plot (within-plot nuisance)", shape = "Sampling month"
    ) +
    theme_bw(base_size = 11)

  fn <- sprintf("output/figures/%s_%s.png", sp_code, rtag)
  ggsave(fn, plt, width = 8, height = 5.5, dpi = 150)
  cat(sprintf("  %-26s slope %+.3f /SD  p = %.3f  (n = %d)  -> %s\n",
              paste(sp_code, rtag), f$slope, f$p, nrow(d), fn))
}

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
                                     plot = as.character(d$plot), month = d$month,
                                     response = rf))
    grids <- rbind(grids, data.frame(ani_tm = f$grid$ani_tm, fit = f$grid$fit,
                                     plot = as.character(f$grid$plot), response = rf))
    annot <- rbind(annot, data.frame(response = rf,
                     label = sprintf("slope %+.3f/SD\np = %.3f  (n = %d)",
                                     f$slope, f$p, nrow(d))))
  }

  plt <- ggplot(long, aes(ani_tm, y, colour = plot)) +
    geom_line(data = grids, aes(ani_tm, fit, colour = plot), linewidth = 0.8) +
    geom_point(aes(shape = month), size = 2.2, alpha = 0.9) +
    geom_text(data = annot, aes(label = label), x = -Inf, y = Inf,
              hjust = -0.04, vjust = 1.15, size = 3, colour = "grey20",
              inherit.aes = FALSE) +
    facet_wrap(~ response, scales = "free_y") +
    labs(
      x = "ANI_disease (time-matched: may = 0, july = 2mai, sep = 4mai)",
      y = NULL, colour = "Plot (within-plot nuisance)", shape = "Sampling month"
    ) +
    theme_bw(base_size = 11)

  fn <- sprintf("output/figures/%s_panel.png", sp_code)
  ggsave(fn, plt, width = 12, height = 4.6, dpi = 150)
  cat(sprintf("  %-18s -> %s\n", sp_name, fn))
}

cat("Writing 12 single-panel figures (3 responses x 4 species):\n")
for (sc in names(SPECIES)) for (rc in names(RESP))
  make_fig(sc, SPECIES[[sc]], rc, RESP[[rc]][1], RESP[[rc]][2])

cat("\nWriting 4 per-species 3-panel figures (length | ratio | germination):\n")
for (sc in names(SPECIES)) make_species_panel(sc, SPECIES[[sc]])

cat("\nDone. Figures in output/figures/:\n")
cat("  single : <species>_<length|ratio|germination>.png  (12)\n")
cat("  panel  : <species>_panel.png                        (4, for the manuscript)\n")
