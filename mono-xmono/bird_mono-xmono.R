## Compare microsat-only with fingerprinting-only and full estimates

rm(list = ls(all = T))

library(tidyverse)
library(ggplot2)
library(RColorBrewer)

## Results of MCMC analyses
load("../analysis_mono/MCMC_resids.rda")
MCMC_resids_mono <- MCMC_resids

load("../analysis_xmono/MCMC_resids.rda")
MCMC_resids_xmono <- MCMC_resids

## Original data
dat1 <- read.csv("../paternity_birds.csv", header = TRUE)
## remove birds with nbrood = NA and avgbrood = NA (can't run without avgbrood)
dat <- dat1[!is.na(dat1$nbrood) & !is.na(dat1$avgbrood),]


#### Plot microsat-only, fingerprint-only, and all null models on single plot ####
MCMC_resids_mono$SocMon <- "mono"
MCMC_resids_xmono$SocMon <- "xmono"

MCMC_resids <- rbind(MCMC_resids_mono, MCMC_resids_xmono)


gg1 <- ggplot() +
  geom_point(data = dat, aes(avgbrood, pmult), alpha = 0.5, size = 2.5) +
  # stat_smooth(data = dat, aes(avgbrood, pmult), color = "black", se = FALSE, method = "lm") +
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult, color = SocMon), 
              se = FALSE, method = "loess") +
  geom_point(data = MCMC_resids, aes(avgbrood, mean.est.pmult, color = SocMon), 
             alpha = 0.5, size = 2.5) +
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult_ucl, color = SocMon), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "loess") +
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult_lcl, color = SocMon), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "loess") +
  labs(x="Brood size", y="Probability of multiple paternity") +
  lims(y = c(0,1)) +
  scale_x_continuous(breaks = c(seq(1,12, by = 1)) , limits = c(1,12)) +
  scale_color_manual(name = "", 
                     labels = c("Socially monogamous", "Not socially monogamous"),
                     values = c("#1B9E77", "#D95F02")) + # "#7570B3"
  theme_bw(base_size = 16)
ggsave("p_monog.eps", plot = gg1, device = cairo_ps, width = 11, height = 7.5, dpi = 320)




#### Conditional probabilities ####
birds_resids <- MCMC_resids

## Probability of MP = 0 given birds were socially monogamous =
## Probability of MP = 0 AND birds were socially monogamous / Probability of birds being socially monogamous 

pAB <- nrow(birds_resids[birds_resids$SocMon=="mono" & birds_resids$pmult==0,])
pA <- nrow(birds_resids[birds_resids$SocMon=="mono",])

pAB/pA

## Probability of MP = 0 given birds were NOT socially monogamous
pAB <- nrow(birds_resids[birds_resids$SocMon=="xmono" & birds_resids$pmult==0,])
pA <- nrow(birds_resids[birds_resids$SocMon=="xmono",])

pAB/pA

## Probability of birds being socially monogamous given MP = 0
pAB <- nrow(birds_resids[birds_resids$pmult==0 & birds_resids$SocMon=="mono",])
pA <- nrow(birds_resids[birds_resids$pmult==0,])

pAB/pA

## Probability of birds being NOT socially monogamous given MP = 0
pAB <- nrow(birds_resids[birds_resids$pmult==0 & birds_resids$SocMon=="xmono",])
pA <- nrow(birds_resids[birds_resids$pmult==0,])

pAB/pA

## The odds of socially monogamous birds to produce an MP = 0 observation
nrow(birds_resids[birds_resids$SocMon=="mono" & birds_resids$pmult==0,])/
  nrow(birds_resids[birds_resids$SocMon=="mono" & !birds_resids$pmult==0,])

## The odds of NON socially monogamous birds to produce an MP = 0 observation
nrow(birds_resids[birds_resids$SocMon=="xmono" & birds_resids$pmult==0,])/
  nrow(birds_resids[birds_resids$SocMon=="xmono" & !birds_resids$pmult==0,])


## Fisher's Exact Test for zero MP versus nonzero MP in each method
dat <- data.frame(
  "Zero" = c(20, 6),
  "Nonzero" = c(66, 24),
  row.names = c("Mono", "XMono"),
  stringsAsFactors = FALSE
)

fisher.test(dat)
#### See cmh_test.R for the Cochran-Mantel-Haenszel test used for this instead



## Clustered Wilcoxon rank sum test for p and pB-p b/w Socially Monog AND NON SocMono
birds_resids$mono_num <- ifelse(birds_resids$SocMon=="xmono", 1, 2)
## Clusters matched to same as CMH test for socmono comparison
birds_resids$brood_id3 <- ifelse((birds_resids$avgbrood - trunc(birds_resids$avgbrood)) < 0.5,
                                 floor(birds_resids$avgbrood), ceiling(birds_resids$avgbrood))
birds_resids$brood_group <- ifelse(birds_resids$brood_id3 >= 7, ">6.5",
                                   ifelse(birds_resids$brood_id3 <=2, "<2.5",
                                          as.character(birds_resids$brood_id3)))
## clusWilcox doesn't play nice with characters or factors - convert to numeric
birds_resids$brood_group <- factor(birds_resids$brood_group)
birds_resids$brood_grpnum <- as.integer(factor(birds_resids$brood_group))+1


library(clusrank)
## (expect non-soc mono to find more MP > 0, therefore larger observed pmult in non-soc mono birds)
## i.e., p_[xmono] - p_[mono] > 0
clusWilcox.test(pmult, cluster = brood_grpnum, group = mono_num, data = birds_resids,
                alternative = "greater", method = "ds", exact = TRUE)
## mean pmult for xmono
mean(birds_resids[birds_resids$mono_num==1,]$pmult)
## mean pmult for mono
mean(birds_resids[birds_resids$mono_num==2,]$pmult)


## (expect non-soc mono to find more MP > 0, therefore smaller residuals)
## i.e., (pB-p)_[xmono] - (pB-p)_[mono] < 0
clusWilcox.test(resid_pmult, cluster = brood_grpnum, group = mono_num, data = birds_resids,
                alternative = "less", method = "ds", exact = TRUE)
## mean pB-p for xmono
mean(birds_resids[birds_resids$mono_num==1,]$resid_pmult)
## mean pB-p for mono
mean(birds_resids[birds_resids$mono_num==2,]$resid_pmult)


## Difference in nsires between methods 
sire_mono <- birds_resids[birds_resids$SocMon=="mono",]$avgsire
sire_xmono <- birds_resids[birds_resids$SocMon=="xmono",]$avgsire

t.test(sire_mono, sire_xmono, alternative = "two.sided", var.equal = FALSE)


## Difference in brood sizes between methods 
broodsz_mono <- birds_resids[birds_resids$SocMon=="mono",]$avgbrood
broodsz_xmono <- birds_resids[birds_resids$SocMon=="xmono",]$avgbrood

t.test(broodsz_mono, broodsz_xmono, alternative = "two.sided", var.equal = FALSE)


## Difference in brood samples between methods (need original data)
birds_data <- readxl::read_xlsx("/Users/Hannah/Documents/Research/Masamu/paternity/birds/paternity_birds.xlsx")
samples_mono <- birds_data[birds_data$SocMon=="mono",]$Nbr
samples_xmono <- birds_data[birds_data$SocMon=="xmono",]$Nbr

t.test(samples_mono, samples_xmono, alternative = "two.sided", var.equal = FALSE)



