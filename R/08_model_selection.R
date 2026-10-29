suppressMessages({ library(readxl); library(readr) })

SEEDS   <- c(bes = 10, jg = 100, wsr = 40, yar = 40)
SPNAME  <- c(bes="Rudbeckia hirta", jg="Sorghum halepense",
             wsr="Ageratina altissima", yar="Achillea millefolium")
RESP    <- c(total_avg="Total length", ratio_avg="Root:shoot ratio", germ_perc="Germination")

ani  <- read_csv("output/tables/ani_disease.csv", show_col_types = FALSE)
resp <- read_excel("bioassay_data_v_final.xlsx", "bioassay_primary")
ani$ani_tm <- 0
ani$ani_tm[ani$month=="july"]      <- ani$ani_disease_2mai[ani$month=="july"]
ani$ani_tm[ani$month=="september"] <- ani$ani_disease_4mai[ani$month=="september"]
respcols <- as.vector(t(outer(names(SEEDS), names(RESP), paste, sep="_")))
dat <- merge(ani[,c("pooled_sample","plot","month","ani_tm")],
             resp[,c("pooled_sample",respcols)], by="pooled_sample")
dat$plot  <- factor(dat$plot)
dat$z_ani <- as.numeric(scale(dat$ani_tm))

aicc  <- function(m){ ll<-logLik(m); k<-attr(ll,"df"); n<-nobs(m); as.numeric(-2*ll+2*k+(2*k*(k+1))/(n-k-1)) }
qaicc <- function(m,chat,n){ ll<-logLik(m); k<-attr(ll,"df")+1; as.numeric(-2*as.numeric(ll)/chat+2*k+(2*k*(k+1))/(n-k-1)) }

rows <- list()
for (sc in names(SEEDS)) for (rc in names(RESP)) {
  col <- paste0(sc,"_",rc); d <- dat[!is.na(dat[[col]]),]; d$y <- d[[col]]; d$plot <- droplevels(d$plot)
  if (rc=="germ_perc") {
    ns <- SEEDS[[sc]]; d$germ <- round(d$y*ns); d$fail <- ns-d$germ; n <- nrow(d)
    b0 <- glm(cbind(germ,fail)~plot, family=binomial, data=d)
    b1 <- glm(cbind(germ,fail)~z_ani+plot, family=binomial, data=d)
    chat <- max(sum(residuals(b1,"pearson")^2)/df.residual(b1), 1)
    ic0 <- qaicc(b0,chat,n); ic1 <- qaicc(b1,chat,n)
    k0 <- attr(logLik(b0),"df")+1; k1 <- attr(logLik(b1),"df")+1; ict <- "QAICc"
  } else {
    n <- nrow(d)
    m0 <- lm(y~plot, data=d); m1 <- lm(y~z_ani+plot, data=d)
    ic0 <- aicc(m0); ic1 <- aicc(m1)
    k0 <- attr(logLik(m0),"df"); k1 <- attr(logLik(m1),"df"); ict <- "AICc"
  }
  ics <- c(ic0, ic1); dd <- ics - min(ics); wt <- exp(-0.5*dd); wt <- wt/sum(wt)
  for (j in 1:2) rows[[length(rows)+1]] <- data.frame(
    species=SPNAME[[sc]], response=RESP[[rc]],
    model=c("null (plot only)","disease-influence")[j],
    ic_type=ict, K=c(k0,k1)[j], n=n, IC=round(ics[j],2),
    dAICc=round(dd[j],2), weight=round(wt[j],3),
    best=ifelse(dd[j]==0,"*",""))
}
tab <- do.call(rbind, rows)
write.csv(tab, "output/tables/model_selection.csv", row.names=FALSE)

cat("\nModel selection (per species x response): null [response ~ plot] vs disease-influence [response ~ z_ani + plot]\n")
cat("IC = AICc (length, ratio) / QAICc (germination); dAICc within each pair; weight = Akaike weight; * = better model\n\n")
print(tab, row.names=FALSE)
ref <- read.csv("output/tables/bioassay_regression.csv")
chk <- merge(
  reshape(tab[,c("species","response","model","IC")], idvar=c("species","response"),
          timevar="model", direction="wide"),
  data.frame(species=unname(SPNAME[ref$species]),
             response=unname(RESP[c(length="total_avg", ratio="ratio_avg", germination="germ_perc")[ref$response]]),
             dAICc_ref=ref$dAICc),
  by=c("species","response"))
stopifnot(nrow(chk) == 12)
chk$dIC <- chk$`IC.null (plot only)` - chk$`IC.disease-influence`
max_diff <- max(abs(chk$dIC - chk$dAICc_ref))
cat(sprintf("\nSanity vs bioassay_regression.csv dAICc: %d of 12 matched | max|diff| = %.3f\n",
            nrow(chk), max_diff))
stopifnot(max_diff < 0.011)
cat("Wrote output/tables/model_selection.csv\n")
