# 07_combined_panel.R
# ---------------------------------------------------------------------------
# Single 12-panel manuscript figure (Figure 5): 4 species (rows) x 3 responses
# (columns) of the within-plot fit used in 02/03. Each panel is independent
# (free y) and annotated with its own within-plot slope/p, so the one figure
# replaces the four per-species panels and removes the "as in Figure 1" shorthand.
#
# Reuses the exact within-plot fit of 03_response_figures.R (Gaussian lm for
# length/ratio; quasibinomial glm for germination; plot partialled out).
#
# Output: output/figures/combined_panel.png
# ---------------------------------------------------------------------------

suppressMessages({ library(readxl); library(readr); library(ggplot2) })

SEEDS   <- c(bes = 10, jg = 100, wsr = 40, yar = 40)
SPECIES <- c(bes = "Black-eyed Susan", jg = "Johnsongrass",
             wsr = "White snakeroot",  yar = "Yarrow")
RESP <- list(total_avg = c("Total seedling length (mm)", "length"),
             ratio_avg = c("Root-to-shoot ratio",        "ratio"),
             germ_perc = c("Germination (proportion)",   "germination"))

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

fit_within <- function(d, sp_code, rcode) {
  if (rcode == "germ_perc") {
    n_seed <- SEEDS[[sp_code]]
    d$germ <- round(d$y * n_seed); d$fail <- n_seed - d$germ
    m <- glm(cbind(germ, fail) ~ z_ani + plot, family = quasibinomial, data = d)
    predfun <- function(g) predict(m, newdata = g, type = "response")
  } else {
    m <- lm(y ~ z_ani + plot, data = d)
    predfun <- function(g) predict(m, newdata = g)
  }
  co <- summary(m)$coefficients["z_ani", ]
  grid <- do.call(rbind, lapply(levels(d$plot), function(pl) {
    rng <- range(d$ani_tm[d$plot == pl])
    g <- data.frame(ani_tm = seq(rng[1], rng[2], length.out = 50),
                    plot = factor(pl, levels = levels(d$plot)))
    g$z_ani <- (g$ani_tm - mu) / sdv; g$fit <- predfun(g); g
  }))
  list(slope = co[["Estimate"]], p = co[["Pr(>|t|)"]], grid = grid)
}

# panel order: species (rows) x response (cols), row-major for facet_wrap(ncol=3)
panel_levels <- as.vector(t(outer(unname(SPECIES),
                                  vapply(RESP, `[`, character(1), 1),
                                  function(s, r) paste0(s, "\n", r))))

long <- grids <- annot <- data.frame()
for (sc in names(SPECIES)) {
  for (rc in names(RESP)) {
    rlab <- RESP[[rc]][1]
    col  <- paste0(sc, "_", rc)
    d <- dat[!is.na(dat[[col]]), ]; d$y <- d[[col]]; d$plot <- droplevels(d$plot)
    f <- fit_within(d, sc, rc)
    pl_lab <- factor(paste0(SPECIES[[sc]], "\n", rlab), levels = panel_levels)
    long  <- rbind(long,  data.frame(ani_tm = d$ani_tm, y = d$y,
                  plot = as.character(d$plot), month = d$month, panel = pl_lab))
    grids <- rbind(grids, data.frame(ani_tm = f$grid$ani_tm, fit = f$grid$fit,
                  plot = as.character(f$grid$plot), panel = pl_lab))
    annot <- rbind(annot, data.frame(panel = pl_lab,
                  label = sprintf("slope %+.3f/SD\np = %.3f", f$slope, f$p)))
  }
}

plt <- ggplot(long, aes(ani_tm, y, colour = plot)) +
  geom_line(data = grids, aes(ani_tm, fit, colour = plot), linewidth = 0.7) +
  geom_point(aes(shape = month), size = 1.9, alpha = 0.9) +
  geom_text(data = annot, aes(label = label), x = -Inf, y = Inf,
            hjust = -0.05, vjust = 1.12, size = 2.6, colour = "grey20",
            inherit.aes = FALSE) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  labs(x = "Disease influence index (ANI_disease; time-matched: May = 0, July = 2 mai, Sep = 4 mai)",
       y = NULL, colour = "Plot (within-plot nuisance)", shape = "Sampling month") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8),
        panel.spacing = unit(0.5, "lines"))

ggsave("output/figures/combined_panel.png", plt, width = 9.5, height = 11, dpi = 200)
cat("Wrote output/figures/combined_panel.png\n")
