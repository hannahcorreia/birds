## Compare microsat-only with fingerprinting-only and full estimates

rm(list = ls(all = T))

library(tidyverse)
library(ggplot2)
library(RColorBrewer)

## Results of MCMC analyses
load("../analysis_all/MCMC_resids.rda")
MCMC_resids_all <- MCMC_resids

load("../analysis_microsat/MCMC_resids.rda")
MCMC_resids_microsat <- MCMC_resids

load("../analysis_fingerprint/MCMC_resids.rda")
MCMC_resids_fingerprint <- MCMC_resids

## Original data
dat1 <- read.csv("../paternity_birds.csv", header = TRUE)
## remove birds with nbrood = NA and avgbrood = NA (can't run without avgbrood)
dat <- dat1[!is.na(dat1$nbrood) & !is.na(dat1$avgbrood),]


#### Plot microsat-only, fingerprint-only, and all null models on single plot ####
MCMC_resids_all$DNA <- "all"
MCMC_resids_microsat$DNA <- "microsatellite"
MCMC_resids_fingerprint$DNA <- "fingerprint"

MCMC_resids <- rbind(rbind(MCMC_resids_all, MCMC_resids_microsat), MCMC_resids_fingerprint)


gg1 <- ggplot() +
  geom_point(data = dat, aes(avgbrood, pmult), alpha = 0.5, size = 2.5) +
  # stat_smooth(data = dat, aes(avgbrood, pmult), color = "black", se = FALSE, method = "lm") +
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult, color = DNA), 
              se = FALSE, method = "loess") +
  geom_point(data = MCMC_resids, aes(avgbrood, mean.est.pmult, color = DNA), 
             alpha = 0.5, size = 2.5) +
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult_ucl, color = DNA), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "loess") +
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult_lcl, color = DNA), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "loess") +
  labs(x="Brood size", y="Probability of multiple paternity") +
  lims(y = c(0,1)) +
  scale_x_continuous(breaks = c(seq(1,12, by = 1)) , limits = c(1,12)) +
  scale_color_manual(name = "Technique", 
                     labels = c("Both", "DNA fingerprinting", "Microsatellite"),
                     values = c("#1B9E77", "#D95F02", "#7570B3")) +
  theme_bw(base_size = 16)
ggsave("p_methods.eps", plot = gg1, device = cairo_ps, width = 11, height = 7.5, dpi = 320)



#### Meta-analysis of microsat vs fingerprinting ####
library(metafor)
library(metaviz)

birds_resids <- rbind(MCMC_resids_microsat, MCMC_resids_fingerprint)

birds_singlerun1 <- merge(dat[,c(1:3)], birds_resids, by = c("species", "avgbrood"))
# order by brood size
birds_singlerun <- birds_singlerun1[order(birds_singlerun1$avgbrood, birds_singlerun1$pmult, birds_singlerun1$nbrood),]

### variance p(1-p)/n for fixed p assumption
birds_singlerun$vi <- birds_singlerun$pmult*(1-birds_singlerun$pmult)/birds_singlerun$nbrood

## labels
bird_species <- birds_singlerun[order(birds_singlerun$avgbrood),]$species

# paternity.meta.bayes <- rma(resid_pmult, vi, data = birds_singlerun, slab = bird_species)
# paternity.meta.bayes

## Fit two separate random-effects models within each subset of DNA
meta1 <- rma(resid_pmult, vi, data = birds_singlerun, slab = bird_species, 
             subset = DNA=="microsatellite")
meta2 <- rma(resid_pmult, vi, data = birds_singlerun, slab = bird_species, 
             subset = DNA=="fingerprint")

meta_dat <- data.frame(estimate = c(coef(meta1), coef(meta2)), 
                       stderror = c(meta1$se, meta2$se),
                       DNA = c("microsatellite","fingerprint"), 
                       tau2 = round(c(meta1$tau2, meta2$tau2),3))
meta_dat

## Compare the two estimates (i.e., the estimated average resids) by feeding them back 
## to the rma() function and using the variable to distinguish the two estimates as a 
## moderator. We use a fixed-effects model, because the (residual) heterogeneity within 
## each subset has already been accounted for by fitting random-effects models above
rma(estimate, sei = stderror, mods = ~ DNA, method = "FE", data = meta_dat, digits = 3)

## So, studies using the microsatellite technique have significantly smaller 
## effects than studies using the fingerprinting technique; i.e., 
## pB-p resids for microsat data are significantly (p<0.001) smaller on average 
## than resids for fingerprinting by 0.082 

## Meta-analysis with DNA as a moderator instead...
rma(resid_pmult, vi, mods = ~ DNA, slab = bird_species, data = birds_singlerun, digits = 3)

## Pretty much the same results as the comparison between the two separate meta-analyses; 
## though an important difference to note between the methods
## The mixed-effects meta-regression model fitted above has a single variance component 
## for the amount of residual heterogeneity, which implies that the amount of 
## heterogeneity within each subset is assumed to be the same. The separate 
## random-effects models in the two subsets allows the amount of heterogeneity 
## within each set to be different. We can do the same thing using the following:
birds_singlerun$pop <- 1:nrow(birds_singlerun)
meta.tog <- rma.mv(resid_pmult, vi, mods = ~ DNA, random = ~ DNA | pop, struct = "DIAG", 
                   slab = bird_species, data = birds_singlerun, digits = 3)
meta.tog


pdf(file = paste0("forest_fp-v-ms.pdf"), width = 11, height = 35)
print(
  forest(meta.tog, cex = 1, xlim = c(-2.5,2), slab = bird_species,
         xlab = expression(paste(p[B], " - p"))),
  quote = FALSE,
  text(-2.5, 246, "Species",  pos=4),
  text(2, 246, expression(paste(p[B], " - p [95% CI]")), pos=2)
)
dev.off()

## better plot with metaviz
## split labels to get species and DNA technique
spnames <- gsub("\\.[0-9]+$", "\\1", unlist(attributes(meta.tog$yi)$slab, use.names = F))
birds_togmeta <- data.frame(resids = as.vector(meta.tog$yi), vi = meta.tog$vi, 
                            species = spnames,
                            DNA = meta.tog$mf.r[[1]][,1] )

# viz_forest(x = birds_togmeta[,c("resids", "vi")], 
#            group = birds_togmeta[,c("DNA")], 
#            study_labels = birds_togmeta[,c("species")], 
#            summary_label = c("Summary (Microsat-only)", "Summary (Both techniques)"), 
#            xlab = expression(paste(p[B], " - p [95% CI]")),
#            # col = "Greys",
#            col = c("firebrick", "steelblue4")[birds_togmeta[,"DNA"]],
#            summary_col = c("firebrick", "steelblue4"),
#            variant = "classic")

vp1 <- viz_forest(meta.tog, 
                  group = birds_togmeta[,"DNA"], 
                  study_labels = birds_togmeta[,"species"], 
                  summary_label = c("Summary (Fingerprint-only)", "Summary (Microsat-only)"), 
                  xlab = expression(paste(p[B], " - p [95% CI]")),
                  # col = "Greys",
                  col = c("firebrick", "steelblue4")[as.factor(birds_togmeta[,"DNA"])],
                  summary_col = c("firebrick", "steelblue4"),
                  variant = "classic")
ggsave("forest_fp-v-ms_snazzy.pdf", plot = vp1, device = "pdf", 
       dpi = 320, width = 11, height = 35, units = "in")



## Let's look at the difference in effect size between whether we include all the data or
## just the microsat-only data to estimate pB-p
birds_resids2 <- rbind(MCMC_resids_microsat, MCMC_resids_all)

birds_all2 <- merge(dat[,c(1:3)], birds_resids2, by = c("species", "avgbrood"))
# order by brood size
birds_all <- birds_all2[order(birds_all2$avgbrood, birds_all2$pmult, birds_all2$nbrood),]

### variance p(1-p)/n for fixed p assumption
birds_all$vi <- birds_all$pmult*(1-birds_all$pmult)/birds_all$nbrood

## labels
bird_species2 <- birds_all[order(birds_all$avgbrood),]$species

birds_all$pop <- 1:nrow(birds_all)
meta.tog <- rma.mv(resid_pmult, vi, mods = ~ DNA, random = ~ DNA | pop, struct = "DIAG", 
                   slab = bird_species2, data = birds_all, digits = 3)
meta.tog
## So, using the microsat-only data would not have a significant effect on the pB-p resids;
## i.e., pB-p resids for microsat-only data are slightly smaller (by 0.038) on average 
## than resids for using all the data (microsats plus fingerprinting)  


pdf(file = paste0("forest_ms-v-all.pdf"), width = 11, height = 55)
print(
  forest(meta.tog, cex = 1, xlim = c(-2.5,2), slab = bird_species2,
         xlab = expression(paste(p[B], " - p"))),
  quote = FALSE,
  text(-2.5, 382, "Species",  pos=4),
  text(2, 382, expression(paste(p[B], " - p [95% CI]")), pos=2)
)
dev.off()


spnames <- gsub("\\.[0-9]+$", "\\1", unlist(attributes(meta.tog$yi)$slab, use.names = F))
birds_togmeta <- data.frame(resids = as.vector(meta.tog$yi), vi = meta.tog$vi, 
                            species = spnames,
                            DNA = meta.tog$mf.r[[1]][,1] )
vp1 <- viz_forest(meta.tog, 
                  group = birds_togmeta[,"DNA"], 
                  study_labels = birds_togmeta[,"species"], 
                  summary_label = c("Summary (Both techniques)", "Summary (Microsat-only)"), 
                  xlab = expression(paste(p[B], " - p [95% CI]")),
                  # col = "Greys",
                  col = c("firebrick", "steelblue4")[as.factor(birds_togmeta[,"DNA"])],
                  summary_col = c("firebrick", "steelblue4"),
                  variant = "classic")
ggsave("forest_ms-v-all_snazzy.pdf", plot = vp1, device = "pdf", 
       dpi = 320, width = 11, height = 49, units = "in")


#### NOTE: MS vs ALL probably not a valid comparison test to determine 
#### if FP should be removed from analyses (direct MS-FP comparison better)




#### Hedges' g meta-analysis ####
birdsg <- within(birds_all,
                 {m1i <- mean.est.pmult
                 m2i <- pmult
                 sd1i <- sqrt(pmult*(1-pmult))
                 n1i <- nbrood})

meta.g1 <- rma(measure="SMD", 
               m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i, sd2i = sd1i, n2i = n1i, 
               slab = bird_species2, data = birdsg, subset = DNA=="microsatellite")
meta.g2 <- rma(measure="SMD", 
               m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i, sd2i = sd1i, n2i = n1i, 
               slab = bird_species2, data = birdsg, subset = DNA=="all")

meta_datg <- data.frame(estimate = c(coef(meta.g1), coef(meta.g2)), 
                        stderror = c(meta.g1$se, meta.g2$se),
                        DNA = c("microsatellite","fingerprint"), 
                        tau2 = round(c(meta.g1$tau2, meta.g2$tau2),3))
meta_datg

rma(estimate, sei = stderror, mods = ~ DNA, method = "FE", data = meta_datg, digits = 3)
## Similar to the results for pB-p, using the microsat-only data would not have a 
## significant effect on the Hedges' g effect estimate: Hedges' g for microsat-only data 
## are slightly smaller (by 0.104) on average than Hedges' g for both DNA together





#### Conditional probabilities ####

## Probability of MP = 0 given data used fingerprinting = 
## Probability of MP = 0 AND data used fingerprinting / Probability of data used fingerprinting

pAB <- nrow(birds_resids[birds_resids$DNA=="fingerprint" & birds_resids$pmult==0,])
pA <- nrow(birds_resids[birds_resids$DNA=="fingerprint",])

pAB/pA

## Probability of MP = 0 given data used microsat
pAB <- nrow(birds_resids[birds_resids$DNA=="microsatellite" & birds_resids$pmult==0,])
pA <- nrow(birds_resids[birds_resids$DNA=="microsatellite",])

pAB/pA

## Probability of data used fingerprinting given MP = 0
pAB <- nrow(birds_resids[birds_resids$pmult==0 & birds_resids$DNA=="fingerprint",])
pA <- nrow(birds_resids[birds_resids$pmult==0,])

pAB/pA

## Probability of data used microsatellite given MP = 0
pAB <- nrow(birds_resids[birds_resids$pmult==0 & birds_resids$DNA=="microsatellite",])
pA <- nrow(birds_resids[birds_resids$pmult==0,])

pAB/pA

## The odds of the fingerprint method to produce an MP = 0 observation
nrow(birds_resids[birds_resids$DNA=="fingerprint" & birds_resids$pmult==0,])/
  nrow(birds_resids[birds_resids$DNA=="fingerprint" & !birds_resids$pmult==0,])

## The odds of the microsat method to produce an MP = 0 observation
nrow(birds_resids[birds_resids$DNA=="microsatellite" & birds_resids$pmult==0,])/
  nrow(birds_resids[birds_resids$DNA=="microsatellite" & !birds_resids$pmult==0,])


## Fisher's Exact Test for zero MP versus nonzero MP in each method
dat <- data.frame(
  "Zero" = c(48, 29),
  "Nonzero" = c(108, 136),
  row.names = c("Fingerprint", "Microsat"),
  stringsAsFactors = FALSE
)

fisher.test(dat)
#### See cmh_test.R for the Cochran-Mantel-Haenszel test used for this instead


# ## For MP > 0, difference in MP between methods 
# pmult_ms <- birds_resids[birds_resids$DNA=="microsatellite" & birds_resids$pmult>0,]$pmult
# pmult_fp <- birds_resids[birds_resids$DNA=="fingerprint" & birds_resids$pmult>0,]$pmult
# 
# t.test(pmult_ms, pmult_fp, alternative = "two.sided", var.equal = FALSE)

## For all MP, difference in MP between methods - Ha is we expect microsat to find more MP 
pmult_ms <- birds_resids[birds_resids$DNA=="microsatellite",]$pmult
pmult_fp <- birds_resids[birds_resids$DNA=="fingerprint",]$pmult
# 
# t.test(pmult_ms, pmult_fp, alternative = "two.sided", var.equal = FALSE)

## Should probably do clustered Wilcoxon rank sum test for p and pB-p b/w DNA techniques
birds_resids$DNA_num <- ifelse(birds_resids$DNA=="microsatellite", 1, 2)
# birds_resids$brood_id2 <- ifelse((birds_resids$avgbrood - trunc(birds_resids$avgbrood)) < 0.25, 
#                                  floor(birds_resids$avgbrood), 
#                                  ifelse((birds_resids$avgbrood - trunc(birds_resids$avgbrood)) >= 0.25 &
#                                           (birds_resids$avgbrood - trunc(birds_resids$avgbrood)) < 0.75, 
#                                         floor(birds_resids$avgbrood) + 0.5, 
#                                         ceiling(birds_resids$avgbrood)))
## Clusters matched to same as CMH test for DNA method comparison
birds_resids$brood_id3 <- ifelse((birds_resids$avgbrood - trunc(birds_resids$avgbrood)) < 0.5,
                                 floor(birds_resids$avgbrood), ceiling(birds_resids$avgbrood))
birds_resids$brood_group <- ifelse(birds_resids$brood_id3 >= 7, ">6.5",
                                   ifelse(birds_resids$brood_id3 <=2, "<=2",
                                          as.character(birds_resids$brood_id3)))
## clusWilcox doesn't play nice with characters or factors - convert to numeric
birds_resids$brood_group <- factor(birds_resids$brood_group)
birds_resids$brood_grpnum <- as.integer(factor(birds_resids$brood_group))+1

## (expect microsat to find more MP > 0, therefore larger observed p > 0 in microsat data?)
library(clusrank)
# clusWilcox.test(pmult, cluster = brood_id2, group = DNA_num, data = birds_resids, 
#                 alternative = "greater", method = "ds")
# clusWilcox.test(pmult, cluster = brood_id3, group = DNA_num, data = birds_resids[birds_resids$pmult>0,], 
#                 alternative = "greater", method = "ds")

## (expect microsat to find more MP >0, therefore larger observed pmult in microsat data)
clusWilcox.test(pmult, cluster = brood_grpnum, group = DNA_num, data = birds_resids,
                alternative = "greater", method = "ds", exact = TRUE)

## Difference in residuals between methods 
resid_ms <- birds_resids[birds_resids$DNA=="microsatellite",]$resid_pmult
resid_fp <- birds_resids[birds_resids$DNA=="fingerprint",]$resid_pmult
# 
# t.test(resid_ms, resid_fp, alternative = "two.sided", var.equal = FALSE)

## (expect microsat to find more MP > 0, therefore smaller residuals?)
## i.e., (pB-p)_microsat - (pB-p)_fingerprint < 0
# clusWilcox.test(resid_pmult, cluster = brood_id2, group = DNA_num, data = birds_resids, 
#                 alternative = "less", method = "ds")
clusWilcox.test(resid_pmult, cluster = brood_grpnum, group = DNA_num, data = birds_resids,
                alternative = "less", method = "ds", exact = TRUE)


## Difference in nsires between methods 
sire_ms <- birds_resids[birds_resids$DNA=="microsatellite",]$avgsire
sire_fp <- birds_resids[birds_resids$DNA=="fingerprint",]$avgsire

t.test(sire_ms, sire_fp, alternative = "two.sided", var.equal = FALSE)


## Difference in brood sizes between methods 
broodsz_ms <- birds_resids[birds_resids$DNA=="microsatellite",]$avgbrood
broodsz_fp <- birds_resids[birds_resids$DNA=="fingerprint",]$avgbrood

t.test(broodsz_ms, broodsz_fp, alternative = "two.sided", var.equal = FALSE)


## Difference in brood samples between methods (need original data)
birds_data <- readxl::read_xlsx("/Users/Hannah/Documents/Research/Masamu/paternity/birds/paternity_birds.xlsx")
samples_ms <- birds_data[birds_data$DNA=="microsatellite",]$Nbr
samples_fp <- birds_data[birds_data$DNA=="fingerprint",]$Nbr

t.test(samples_ms, samples_fp, alternative = "two.sided", var.equal = FALSE)



