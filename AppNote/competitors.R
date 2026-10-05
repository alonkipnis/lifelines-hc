#!/usr/bin/env Rscript
# MaxCombo (Fleming-Harrington dictionary) and Yang-Prentice adaptive log-rank.
#
# Usage:
#   Rscript competitors.R <input.csv> <output.json>
#
# Input CSV columns: time, event, group
#   group = 0 for control (A), 1 for treatment (B).
#
# MaxCombo uses the four-weight immuno-oncology family
#   FH(rho, gamma) in {(0,0), (0,1), (1,0), (1,1)}
# with a two-sided multivariate-normal p-value as in nph::logrank.maxtest
# (Ristl et al. 2021; Lin et al. 2020).
#
# Yang-Prentice reports the adaptive weighted log-rank test of
# Yang & Prentice (2010), built on the short-term/long-term hazard-ratio
# model of Yang & Prentice (2005), via YPmodel::YPmodel.adlgrk.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript competitors.R <input.csv> <output.json>")
}

suppressPackageStartupMessages({
  library(nph)
  library(YPmodel)
  library(jsonlite)
})

infile  <- args[[1]]
outfile <- args[[2]]
seed    <- if (length(args) >= 3) as.integer(args[[3]]) else 42L

df <- read.csv(infile)
stopifnot(all(c("time", "event", "group") %in% names(df)))

time  <- as.numeric(df$time)
event <- as.integer(df$event)
group <- as.integer(df$group)

set.seed(seed)
mc <- logrank.maxtest(
  time, event, group,
  alternative = "two.sided",
  rho   = c(0, 0, 1, 1),
  gamma = c(0, 1, 0, 1)
)
zmax <- max(abs(mc$tests$z))
fh_labels <- c("FH(0,0)", "FH(0,1)", "FH(1,0)", "FH(1,1)")
components <- lapply(seq_len(nrow(mc$tests)), function(i) {
  list(
    weight    = fh_labels[[i]],
    z         = unname(mc$tests$z[[i]]),
    p_value   = unname(mc$tests$p[[i]])
  )
})

ypdat <- data.frame(V1 = time, V2 = event, V3 = group)
est <- YPmodel.estimate(data = ypdat, interval = 0)
data_yp <- YPmodel.inputData(data = ypdat)
adl <- fun.adlgrk(est$beta, est$r, data_yp)
theta <- as.numeric(exp(est$beta))

out <- list(
  maxcombo = list(
    statistic = unname(zmax),
    p_value   = unname(mc$pmult),
    p_bonferroni = unname(mc$p.Bonf),
    components = components
  ),
  yang_prentice = list(
    statistic = unname(as.numeric(adl$t)),
    p_value   = unname(as.numeric(adl$pval)),
    theta_short = unname(theta[[1]]),
    theta_long  = unname(theta[[2]])
  )
)

write(toJSON(out, auto_unbox = TRUE, digits = 12, pretty = TRUE), outfile)
cat("Wrote", outfile, "\n")
