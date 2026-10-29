suppressMessages({
  library(readxl); library(readr); library(ggplot2)
  library(lme4); library(lmerTest); library(glmmTMB)
})
source("R/v2_style.R")

OUT_DIR <- "output/r2/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SEEDS   <- c(bes = 10, jg = 100, wsr = 40, yar = 40)
SPECIES <- c(bes = "Black-eyed Susan", jg = "Johnsongrass",
             wsr = "White snakeroot",  yar = "Yarrow")
RESP <- list(total_avg = c("Total seedling length (mm)", "length"),
             ratio_avg = c("Root-to-shoot ratio",        "ratio"),
             germ_perc = c("Germination (proportion)",   "germination"))

ani  <- read_csv("output/tables/ani_disease.csv", show_col_types = FALSE)
resp <- read_excel("bioassay_data_v_final.xlsx", "bioassay_primary")
ani$ani_tm <- 0
ani$ani_tm[ani$month == "july"]      <- ani$ani_disease_2mai[ani$month == "july"]
ani$ani_tm[ani$month == "september"] <- ani$ani_disease_4mai[ani$month == "september"]

respcols <- as.vector(t(outer(names(SPECIES), names(RESP), paste, sep = "_")))
dat <- merge(ani[, c("pooled_sample", "plot", "month", "distance", "bearing", "ani_tm")],
             resp[, c("pooled_sample", respcols)], by = "pooled_sample")
dat$station <- factor(paste(dat$plot, ifelse(dat$distance == 0, "center",
                                             paste0("b", dat$bearing)), sep = "_"))
dat$plot  <- factor(dat$plot)
dat$month <- factor(dat$month, levels = c("may", "july", "september"))
dat$z_ani <- as.numeric(scale(dat$ani_tm))
mu <- mean(dat$ani_tm); sdv <- sd(dat$ani_tm)
stopifnot(nrow(dat) == 45, nlevels(dat$station) == 15)

partial_panel <- function(d, sp, rc) {
  if (rc == "germ_perc") {
    n <- SEEDS[[sp]]
    d$germ <- round(d$y * n); d$fail <- n - d$germ
    m   <- glmmTMB(cbind(germ, fail) ~ plot + month + z_ani + (1 | station),
                   family = betabinomial, data = d)
    co  <- summary(m)$coefficients$cond["z_ani", ]
    b <- co[["Estimate"]]; se <- co[["Std. Error"]]; p <- co[["Pr(>|z|)"]]
    q   <- qnorm(0.975)
    eta <- predict(m, type = "link")
    ylk <- log((d$germ + 0.5) / (d$fail + 0.5))
    inv <- plogis
  } else {
    m   <- lmerTest::lmer(y ~ plot + month + z_ani + (1 | station),
                          data = d, REML = TRUE)
    co  <- summary(m, ddf = "Satterthwaite")$coefficients["z_ani", ]
    b <- co[["Estimate"]]; se <- co[["Std. Error"]]; p <- co[["Pr(>|t|)"]]
    q   <- qt(0.975, co[["df"]])
    eta <- fitted(m)
    ylk <- d$y
    inv <- identity
  }
  zbar <- mean(d$z_ani); c0 <- mean(eta)
  pts  <- data.frame(ani_tm = d$ani_tm, month = d$month, plot = d$plot,
                     y = inv((ylk - eta) + b * (d$z_ani - zbar) + c0))
  zg   <- seq(min(d$z_ani), max(d$z_ani), length.out = 100)
  lin  <- c0 + b * (zg - zbar); hw <- q * se * abs(zg - zbar)
  line <- data.frame(ani_tm = zg * sdv + mu, fit = inv(lin),
                     lo = inv(lin - hw), hi = inv(lin + hw))
  list(pts = pts, line = line, b = b, p = p, n = nrow(d))
}

pfmt <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
panel_levels <- as.vector(t(outer(unname(SPECIES),
                                  vapply(RESP, `[`, character(1), 1),
                                  function(s, r) paste0(s, "\n", r))))
RCODE <- c(total_avg = "length", ratio_avg = "ratio", germ_perc = "germination")

pts <- lines <- annot <- data.frame(); chk <- data.frame()
for (sc in names(SPECIES)) for (rc in names(RESP)) {
  col <- paste0(sc, "_", rc)
  d <- dat[!is.na(dat[[col]]), ]; d$y <- d[[col]]
  f <- partial_panel(d, sc, rc)
  pl <- factor(paste0(SPECIES[[sc]], "\n", RESP[[rc]][1]), levels = panel_levels)
  pts   <- rbind(pts,   cbind(f$pts,  panel = pl))
  lines <- rbind(lines, cbind(f$line, panel = pl))
  annot <- rbind(annot, data.frame(panel = pl,
                 label = sprintf("slope %+.3f / SD\n%s", f$b, pfmt(f$p))))
  chk   <- rbind(chk, data.frame(species = sc, response = RCODE[[rc]],
                                 slope_fig = f$b, p_fig = f$p))
}

ref <- read.csv("output/r2/regression_r2_primary.csv")
mm  <- merge(chk, ref, by = c("species", "response"))
stopifnot(nrow(mm) == 12,
          max(abs(mm$slope_fig - mm$slope_per_SD) / pmax(abs(mm$slope_per_SD), 1e-8)) < 1e-4,
          max(abs(mm$p_fig - mm$p)) < 1e-4)

pts$plot_lab <- recode_plot(pts$plot)

plt <- ggplot() +
  geom_ribbon(data = lines, aes(ani_tm, ymin = lo, ymax = hi),
              fill = "grey55", alpha = 0.25) +
  geom_line(data = lines, aes(ani_tm, fit), colour = "grey10", linewidth = 0.9) +
  geom_point(data = pts, aes(ani_tm, y, colour = plot_lab, shape = month),
             size = 2.3, alpha = 0.9, stroke = 0.5) +
  geom_text(data = annot, aes(label = label), x = Inf, y = Inf,
            hjust = 1.05, vjust = 1.12, size = 3.6, colour = "grey10",
            lineheight = 0.95) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  scale_plot_colour("Plot") +
  scale_month_shape(labels = c(may = "May", july = "July", september = "September")) +
  labs(x = LAB_ANI, y = "Partial residual (plot, month and station effects removed)") +
  guides(colour = guide_legend(nrow = 1, order = 1),
         shape  = guide_legend(nrow = 1, order = 2)) +
  theme_v2(base_size = 14) +
  theme(strip.text = element_text(size = rel(0.82), face = "bold", colour = "grey10",
                                  lineheight = 0.95, margin = margin(4, 4, 4, 4)),
        panel.spacing = unit(0.8, "lines"))

fn <- file.path(OUT_DIR, "figure4_partial_residuals.png")
ggsave(fn, plt, width = 10.5, height = 13.2, dpi = 300)
cat("Wrote", fn, "| panel slopes and p values match regression_r2_primary.csv\n")
