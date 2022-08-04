## Simple calculations for birds vs. mammals analyses

library(tidyverse)
library(ggplot2)
library(reshape2)
library(irr)

## both DNA methods
# paternity_birds <- read.csv("/Users/Hannah/Documents/Research/Masamu/paternity/birds/analysis_all/paternity_birds.csv") 
## microsat-only
paternity_birds <- read.csv("/Users/Hannah/Documents/Research/Masamu/paternity/birds/analysis_microsat/paternity_birds_ms.csv") 
## mammals - latest from phlyo analysis (64spp)
paternity_mammals <- read.csv("/Users/Hannah/Documents/Research/Masamu/paternity/birds_vs_mammals/mammals_for_comparison/MCMC_estimates_byspecies/phylo_mammals_fixed.csv")


#### Krippendorff's alpha for MP (birds - microsat only) ####

## Create unique species ID number for each species and obs number for each pmult within species
bird_dat <- paternity_birds %>%
  mutate(species2 = factor(species), species_ID2 = factor(paste0("S",factor(as.numeric(species2))))) %>%
  group_by(species_ID2) %>%
  mutate(obs_ID2 = factor(paste0("P", 1:length(pmult)))) %>%
  ungroup() 

## It doesn't make sense to run Krippendorff's alpha on the species with only one observation
bird_mp <- bird_dat %>% select(species_ID2, obs_ID2, pmult) %>%
  group_by(species_ID2) %>% filter(length(pmult)>1) %>% ungroup() %>%
  pivot_wider(names_from = species_ID2, values_from = pmult, id_cols = obs_ID2, 
              values_fill = NA)

mp_mat <- as.matrix(bird_mp[,-1])

## Krippendorff's alpha for ratio data of only species with >1 observation
kripp.alpha(mp_mat, method = "ratio")

# Krippendorff's alpha
# 
#  Subjects = 11 
#    Raters = 5 
#     alpha = 0.402  ## This is very low (>0.8 considered good agreement)

## Can use custom R code from "Additional file 3" of Zapf et al. 2016 (https://doi.org/10.1186/s12874-016-0200-9)
source("k_alpha.R")
## Custom version (with bootstrap 95% CI)
k_alpha(t(mp_mat), alpha_q = 0.05, nboot = 1000, scaling = "ratio")
## So, species with >1 obs had relatively dissimilar observed MP values, apparently


bird_dat %>% select(species_ID2, obs_ID2, pmult) %>%
  group_by(species_ID2) %>% filter(length(pmult)>1) %>% ungroup() %>%
  ggplot(aes(x = species_ID2, y = pmult, fill = obs_ID2)) +
  geom_dotplot(binaxis = "y", stackdir = "center",
               stackratio=1.5, dotsize=0.8, binwidth = 0.02, alpha = 0.6)

## Difference between min(MP) and max(MP) for these species with >1 obs
bird_dat %>% select(species_ID2, obs_ID2, pmult) %>%
  group_by(species_ID2) %>% filter(length(pmult)>1) %>% 
  summarise(sp_max_diff = max(pmult) - min(pmult)) %>%
  print() %>% ungroup() %>%
  summarise(mean_diff = mean(sp_max_diff), 
            min_diff = min(sp_max_diff), max_diff = max(sp_max_diff))



#### Krippendorff's alpha for brood size (birds - microsat only) ####
bird_bs <- bird_dat %>% select(species_ID2, obs_ID2, avgbrood) %>%
  group_by(species_ID2) %>% filter(length(avgbrood)>1) %>% ungroup() %>%
  # mutate(avgbrood = avgbrood/10) %>%
  pivot_wider(names_from = species_ID2, values_from = avgbrood, id_cols = obs_ID2, 
              values_fill = NA)

bs_mat <- as.matrix(bird_bs[,-1])

## Krippendorff's alpha for ratio data of only species with >1 observation
kripp.alpha(bs_mat, method = "ratio")

## Custom version (with bootstrap 95% CI)
k_alpha(t(bs_mat), alpha_q = 0.05, nboot = 1000, scaling = "ratio")


#### Krippendorff's alpha for nsires (birds - microsat only) ####
bird_sire <- bird_dat %>% select(species_ID2, obs_ID2, avgsire) %>%
  group_by(species_ID2) %>% filter(length(avgsire)>1) %>% ungroup() %>%
  # mutate(avgsire = avgsire/10) %>%
  pivot_wider(names_from = species_ID2, values_from = avgsire, id_cols = obs_ID2, 
              values_fill = NA)

sire_mat <- as.matrix(bird_sire[,-1])

## Krippendorff's alpha for ratio data of only species with >1 observation
kripp.alpha(sire_mat, method = "ratio")

## Custom version (with bootstrap 95% CI)
k_alpha(t(sire_mat), alpha_q = 0.05, nboot = 1000, scaling = "ratio")



#### Krippendorff's alpha for pB-p (birds - microsat only) ####
## both DNA methods
# load("/Users/Hannah/Documents/Research/Masamu/paternity/birds/analysis_all/MCMC_resids.rda")
## microsat-only
load("/Users/Hannah/Documents/Research/Masamu/paternity/birds/analysis_microsat/MCMC_resids.rda")
bird_resids <- MCMC_resids

## Create unique species ID number for each species and obs number for each pmult within species
bird_rdat <- bird_resids %>%
  mutate(species2 = factor(species), species_ID2 = factor(paste0("S",factor(as.numeric(species2))))) %>%
  group_by(species_ID2) %>%
  mutate(obs_ID2 = factor(paste0("P", 1:length(pmult)))) %>%
  ungroup() 
## Create matrix
bird_residdat <- bird_rdat %>% select(species_ID2, obs_ID2, resid_pmult) %>%
  group_by(species_ID2) %>% filter(length(resid_pmult)>1) %>% ungroup() %>%
  pivot_wider(names_from = species_ID2, values_from = resid_pmult, id_cols = obs_ID2, 
              values_fill = NA)

resid_mat <- as.matrix(bird_residdat[,-1])

## Krippendorff's alpha for interval data of only species with >1 observation
kripp.alpha(resid_mat, method = "interval")

## Custom version (with bootstrap 95% CI)
k_alpha(t(resid_mat), alpha_q = 0.05, nboot = 1000, scaling = "interval")




#### Krippendorff's alpha for MP (mammals) ####

## Create unique species ID number for each species and obs number for each pmult within species
mammal_dat <- paternity_mammals %>%
  mutate(species2 = factor(species), species_ID2 = factor(paste0("S",factor(as.numeric(species2))))) %>%
  group_by(species_ID2) %>%
  mutate(obs_ID2 = factor(paste0("P", 1:length(pmult)))) %>%
  ungroup() 

## It doesn't make sense to run Krippendorff's alpha on the species with only one observation
mammal_mp <- mammal_dat %>% select(species_ID2, obs_ID2, pmult) %>%
  group_by(species_ID2) %>% filter(length(pmult)>1) %>% ungroup() %>%
  pivot_wider(names_from = species_ID2, values_from = pmult, id_cols = obs_ID2, 
              values_fill = NA)

mp_mat <- as.matrix(mammal_mp[,-1])

## Krippendorff's alpha for ratio data of only species with >1 observation
kripp.alpha(mp_mat, method = "ratio")

# Krippendorff's alpha
# 
#  Subjects = 3 
#    Raters = 2 
#     alpha = 0.65  ## This is somewhat low (>0.8 considered good agreement)

## Custom version (with bootstrap 95% CI)
k_alpha(t(mp_mat), alpha_q = 0.05, nboot = 1000, scaling = "ratio")


## So, species with >1 obs had somewhat dissimilar observed MP values, apparently
mammal_dat %>% select(species_ID2, obs_ID2, pmult) %>%
  group_by(species_ID2) %>% filter(length(pmult)>1) %>% ungroup() %>%
  ggplot(aes(x = species_ID2, y = pmult, fill = obs_ID2)) +
  geom_dotplot(binaxis = "y", stackdir = "center",
               stackratio=1.5, dotsize=0.8, binwidth = 0.02, alpha = 0.6)

## Difference between min(MP) and max(MP) for these species with >1 obs
mammal_dat %>% select(species_ID2, obs_ID2, pmult) %>%
  group_by(species_ID2) %>% filter(length(pmult)>1) %>% 
  summarise(sp_max_diff = max(pmult) - min(pmult)) %>%
  print() %>% ungroup() %>%
  summarise(mean_diff = mean(sp_max_diff), 
            min_diff = min(sp_max_diff), max_diff = max(sp_max_diff))



#### Krippendorff's alpha for brood size (mammals) ####
mammal_bs <- mammal_dat %>% select(species_ID2, obs_ID2, avgbrood) %>%
  group_by(species_ID2) %>% filter(length(avgbrood)>1) %>% ungroup() %>%
  # mutate(avgbrood = avgbrood/10) %>%
  pivot_wider(names_from = species_ID2, values_from = avgbrood, id_cols = obs_ID2, 
              values_fill = NA)

bs_mat <- as.matrix(mammal_bs[,-1])

## Krippendorff's alpha for ratio data of only species with >1 observation
kripp.alpha(bs_mat, method = "ratio")

## Custom version (with bootstrap 95% CI)
k_alpha(t(bs_mat), alpha_q = 0.05, nboot = 1000, scaling = "ratio")



#### Krippendorff's alpha for nsires (mammals) ####
mammal_sire <- mammal_dat %>% select(species_ID2, obs_ID2, avgsire) %>%
  group_by(species_ID2) %>% filter(length(avgsire)>1) %>% ungroup() %>%
  # mutate(avgsire = avgsire/10) %>%
  pivot_wider(names_from = species_ID2, values_from = avgsire, id_cols = obs_ID2, 
              values_fill = NA)

sire_mat <- as.matrix(mammal_sire[,-1])

## Krippendorff's alpha for ratio data of only species with >1 observation
kripp.alpha(sire_mat, method = "ratio")

## Custom version (with bootstrap 95% CI)
k_alpha(t(sire_mat), alpha_q = 0.05, nboot = 1000, scaling = "ratio")



#### Krippendorff's alpha for pB-p (mammals) ####
load("/Users/Hannah/Documents/Research/Masamu/paternity/birds_vs_mammals/mammals_for_comparison/MCMC_estimates_byspecies/phyloMCMC_resids.rda")
mammal_resids <- MCMC_resids 

## Create unique species ID number for each species and obs number for each pmult within species
mammal_rdat <- mammal_resids %>%
  mutate(species2 = factor(species), species_ID2 = factor(paste0("S",factor(as.numeric(species2))))) %>%
  group_by(species_ID2) %>%
  mutate(obs_ID2 = factor(paste0("P", 1:length(pmult)))) %>%
  ungroup() 
## Create matrix
mammal_residdat <- mammal_rdat %>% select(species_ID2, obs_ID2, resid_pmult) %>%
  group_by(species_ID2) %>% filter(length(resid_pmult)>1) %>% ungroup() %>%
  pivot_wider(names_from = species_ID2, values_from = resid_pmult, id_cols = obs_ID2, 
              values_fill = NA)

resid_mat <- as.matrix(mammal_residdat[,-1])

## Krippendorff's alpha for interval data of only species with >1 observation
kripp.alpha(resid_mat, method = "interval")

## Custom version (with bootstrap 95% CI)
k_alpha(t(resid_mat), alpha_q = 0.05, nboot = 1000, scaling = "interval")




#### Plot pB-p vs avgbrood (mammals & birds) ####
## Residuals for both DNA methods
# load("/Users/Hannah/Documents/Research/Masamu/paternity/birds/analysis_all/MCMC_resids.rda")
## Residuals for microsat-only
load("/Users/Hannah/Documents/Research/Masamu/paternity/birds/analysis_microsat/MCMC_resids.rda")
bird_resids <- MCMC_resids
## Residuals for mammals (latest from phylo, 64spp)
load("/Users/Hannah/Documents/Research/Masamu/paternity/birds_vs_mammals/mammals_for_comparison/MCMC_estimates_byspecies/phyloMCMC_resids.rda")
mammal_resids <- MCMC_resids 

mammal_resids$group <- "Mammals"
bird_resids$group <- "Birds"

resid_dat <- rbind(mammal_resids, bird_resids)


## Difference in MP between birds and mammals 
pmult_b <- bird_resids$pmult
pmult_m <- mammal_resids$pmult

t.test(pmult_b, pmult_m, alternative = "two.sided", var.equal = FALSE)


## Difference in avgsires between birds and mammals 
sires_b <- bird_resids$avgsire
sires_m <- mammal_resids$avgsire

t.test(sires_b, sires_m, alternative = "two.sided", var.equal = FALSE)


## Difference in avgbrood between birds and mammals 
broodsz_b <- bird_resids$avgbrood
broodsz_m <- mammal_resids$avgbrood

t.test(broodsz_b, broodsz_m, alternative = "two.sided", var.equal = FALSE)


## Difference in number of broods sampled between birds and mammals (need original data)
birds_data <- readxl::read_xlsx("/Users/Hannah/Documents/Research/Masamu/paternity/birds/paternity_birds.xlsx")[,-1]
mammals_data <- read.csv("/Users/Hannah/Documents/Research/Masamu/paternity/birds_vs_mammals/mammals_for_comparison/MCMC_estimates_byspecies/phylo_mammals_fixed.csv")

samples_b <- birds_data[birds_data$DNA=="microsatellite",]$Nbr
samples_m <- mammals_data[!is.na(mammals_data$avgsire),]$nbrood

t.test(samples_b, samples_m, alternative = "two.sided", var.equal = FALSE)


## Difference in residuals between birds and mammals 
resid_b <- bird_resids$resid_pmult
resid_m <- mammal_resids$resid_pmult

t.test(resid_b, resid_m, alternative = "two.sided", var.equal = FALSE)


## Fisher's Exact Test for zero MP versus nonzero MP in each taxa
dat <- data.frame(
  "Zero" = c(29, 2),
  "Nonzero" = c(136, 63),
  row.names = c("Birds", "Mammals"),
  stringsAsFactors = FALSE
)

fisher.test(dat)




## Plot pB for mammals and birds on the same plot (just for curiosity)
gp0 <- ggplot() +
  geom_point(data = resid_dat, aes(avgbrood, pmult, color = group), alpha = 0.5, size = 2.5) +
  # stat_smooth(data = resid_dat, aes(avgbrood, pmult), se = FALSE, color = "black", method = "lm") +
  stat_smooth(data = resid_dat, aes(avgbrood, pred_pmult, color = group), 
              alpha = 1, se = FALSE, method = "loess") +
  geom_point(data = resid_dat, aes(avgbrood, mean.est.pmult, color = group), 
             alpha = 0.5, size = 2.5) +
  stat_smooth(data = resid_dat, aes(avgbrood, pred_pmult_ucl, color = group), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "loess") +
  stat_smooth(data = resid_dat, aes(avgbrood, pred_pmult_lcl, color = group), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "loess") +
  labs(x="Brood/litter size", y="Probability of multiple paternity") +
  lims(y = c(0,1)) +
  scale_color_manual(name = "",  values = c("#1B9E77", "#D95F02")) +
  scale_x_continuous(breaks = c(seq(1,12, by = 1)) , limits = c(1,12)) +
  theme_bw(base_size = 16)

ggsave("p_mamm-bird.eps", plot = gp0, device = cairo_ps, width = 11, height = 7.5, dpi = 320)


## Plot of pB-p versus avgbrood
gp1 <- ggplot() +
  # geom_point(data = dat, aes(avgbrood, pmult), alpha = 0.5, size = 2.5) +
  # stat_smooth(data = dat, aes(avgbrood, pmult), color = "black", se = FALSE, method = "lm") +
  stat_smooth(data = resid_dat, aes(avgbrood, resid_pmult, color = group), 
              se = FALSE, method = "lm") +
  geom_point(data = resid_dat, aes(avgbrood, resid_pmult, color = group), 
             alpha = 0.5, size = 2.5) +
  stat_smooth(data = resid_dat, aes(avgbrood, resid_pmult_ucl, color = group), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "lm") +
  stat_smooth(data = resid_dat, aes(avgbrood, resid_pmult_lcl, color = group), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "lm") +
  labs(x="Brood/litter size", y="Probability of multiple paternity") +
  lims(y = c(-1,1)) +
  scale_x_continuous(breaks = c(seq(1,12, by = 1)) , limits = c(1,12)) +
  scale_color_manual(name = "",  #labels = c(),
                     values = c("#1B9E77", "#D95F02")) +
  theme_bw(base_size = 16)

gp1 # The outlier (avgbrood >10) for birds makes me queasy...


# Histogram
hist(resid_dat[resid_dat$group=="Mammals",]$resid_pmult)
hist(resid_dat[resid_dat$group=="Birds",]$resid_pmult)

# LM by group
lmres <- lm(resid_pmult ~ avgbrood*group, data = resid_dat)
summary(lmres)
# plot(lmres)

cor(predict(lmres), residuals(lmres))

predls <- within(lmres$model, pred <- predict(lmres))
AIC(lmres)

gp_ls <- ggplot(data = predls, aes(avgbrood, resid_pmult, color = group)) +
  geom_point(alpha = 0.5) + geom_line(aes(y = pred)) +
  scale_color_manual(name = "",  values = c("#1B9E77", "#D95F02")) +
  labs(title = "LS", x="Brood/litter size", y=expression(paste(p[B], " - p"))) +
  scale_x_continuous(breaks = c(seq(1,12, by = 2)) , limits = c(1,12)) +
  theme_bw(base_size = 16) + theme(legend.position="bottom")


## LS LM by group without outliers
lmres_no <- lm(resid_pmult ~ avgbrood*group, data = resid_dat[resid_dat$avgbrood<10,])
summary(lmres_no)
# plot(lmres_no)

cor(predict(lmres_no), residuals(lmres_no))

predls_no <- within(lmres_no$model, pred <- predict(lmres_no))
AIC(lmres_no)

gp_ls_no <- ggplot(data = predls_no, aes(avgbrood, resid_pmult, color = group)) +
  geom_point(alpha = 0.5) + geom_line(aes(y = pred)) +
  scale_color_manual(name = "",  values = c("#1B9E77", "#D95F02")) +
  labs(title = "LS - 1 outlier removed", x="Brood/litter size", y=expression(paste(p[B], " - p"))) +
  scale_x_continuous(breaks = c(seq(1,12, by = 2)) , limits = c(1,12)) +
  theme_bw(base_size = 16)  + theme(legend.position="bottom")


# Rank LM by group
library(Rfit)
rlmres <- rfit(resid_pmult ~ avgbrood*group, data = resid_dat)
summary(rlmres)

## Residuals plot
plot(rlmres$fitted.values, rlmres$residuals)
cor(rlmres$fitted.values, rlmres$residuals)

predrank <- data.frame(avgbrood = rlmres$x[,2],
                       resid_pmult = rlmres$y,
                       groupMammals = rlmres$x[,3],
                       pred = rlmres$fitted.values)
predrank$group <- ifelse(predrank$groupMammals=="1", "Mammals", "Birds")

gp_r <- ggplot(data = predrank, aes(avgbrood, resid_pmult, color = group)) +
  geom_point(alpha = 0.5) + geom_line(aes(y = pred)) +
  scale_color_manual(name = "",  values = c("#1B9E77", "#D95F02")) +
  labs(title = "Rank", x="Brood/litter size", y=expression(paste(p[B], " - p"))) +
  scale_x_continuous(breaks = c(seq(1,12, by = 2)) , limits = c(1,12)) +
  theme_bw(base_size = 16) + theme(legend.position="bottom")


library(gridExtra)

gg3 <- grid.arrange(gp_ls, gp_r, gp_ls_no, ncol = 3)

ggsave("resids_mamm-bird.eps", plot = gg3, device = cairo_ps, width = 15, height = 5, dpi = 320)


## save rank version of mammals vs birds
gp_r2 <- ggplot(data = predrank, aes(avgbrood, resid_pmult, color = group)) +
  geom_point(alpha = 0.5, size = 2.5) + geom_line(aes(y = pred), size = 1.5) +
  scale_color_manual(name = "",  values = c("#1B9E77", "#D95F02")) +
  labs(x="Brood/litter size", y=expression(paste(p[B], " - p"))) +
  scale_x_continuous(breaks = c(seq(1,12, by = 2)) , limits = c(1,12)) +
  theme_bw(base_size = 16) #+ theme(legend.position="top")

ggsave("rank_mamm-bird.eps", plot = gp_r2, device = cairo_ps, width = 11, height = 7.5, dpi = 320)



## Non-parametric comparison of mammals and birds in binned broodsizes 
#### (rank-sum test for clustered data, Datta and Satten 2005) ####
library(clusrank)

## Create "bins" of avgbrood sizes with whole values, 
## e.g., 0.5 <= bin("1") < 1.5; 1.5 <= bin("2") < 2.5; and so forth
resid_dat$brood_id <- ifelse((resid_dat$avgbrood - trunc(resid_dat$avgbrood)) < 0.5, 
                             floor(resid_dat$avgbrood), ceiling(resid_dat$avgbrood))
## "ds" method (i.e., Datta and Satten, informative clusters) does not behave with
## groups that are factors or characters; need numeric groups
resid_dat$group_num <- ifelse(resid_dat$group=="Birds", 1, 2) ## bird = 1, mammals = 2

## birds vs mammals binned by avgbrood - nonparametric test asks 
## "Are values of pB-p for birds signif diff from those of mammals?"
clusWilcox.test(resid_pmult, cluster = brood_id, group = group_num, data = resid_dat, 
                alternative = "two.sided", method = "ds")

## "Are values of pB-p for birds (group_num=1) signif larger than those of mammals (group_num=2)?"
## i.e. (pB-p)_birds - (pB-p)_mammals > 0
## We would need evidence to suggest that we expect bird resids to be smaller than mamm resids
clusWilcox.test(resid_pmult, cluster = brood_id, group = group_num, data = resid_dat, 
                alternative = "greater", method = "ds", exact = TRUE)


## Consider a different number of bins (each smaller)
# ## e.g., 0.75 <= bin("1") < 1.25; 1.25 <= bin("1.5") < 1.75; and so forth
# resid_dat$brood_id2 <- ifelse((resid_dat$avgbrood - trunc(resid_dat$avgbrood)) < 0.25, 
#                               floor(resid_dat$avgbrood), 
#                               ifelse((resid_dat$avgbrood - trunc(resid_dat$avgbrood)) >= 0.25 &
#                                        (resid_dat$avgbrood - trunc(resid_dat$avgbrood)) < 0.75, 
#                                      floor(resid_dat$avgbrood) + 0.5, 
#                                      ceiling(resid_dat$avgbrood)))

## "ds" method (i.e., Datta and Satten, informative clusters) does not behave with
## groups that are factors or characters; need numeric groups
resid_dat$group_num <- ifelse(resid_dat$group=="Birds", 1, 2) ## bird = 1, mammals = 2

## Clusters matched to same as CMH test for DNA method comparison
resid_dat$brood_id3 <- ifelse((resid_dat$avgbrood - trunc(resid_dat$avgbrood)) < 0.5,
                              floor(resid_dat$avgbrood), ceiling(resid_dat$avgbrood))
resid_dat$brood_group <- ifelse(resid_dat$brood_id3 >= 7, ">6.5",
                                ifelse(resid_dat$brood_id3 <=2, "<=2",
                                       as.character(resid_dat$brood_id3)))
## clusWilcox doesn't play nice with characters or factors - convert to numeric
resid_dat$brood_group <- factor(resid_dat$brood_group)
resid_dat$brood_grpnum <- as.integer(factor(resid_dat$brood_group))+1

## "Are values of pmult for birds (group_num=1) signif smaller than those of mammals (group_num=2)?"
## i.e. pmult_birds - pmult_mammals < 0
clusWilcox.test(pmult, cluster = brood_grpnum, group = group_num, data = resid_dat, 
                alternative = "less", method = "ds", exact = TRUE)

## "Are values of pB-p for birds (group_num=1) signif larger than those of mammals (group_num=2)?"
## i.e. (pB-p)_birds - (pB-p)_mammals > 0
## For smaller bins (larger number of bins @ 19)
clusWilcox.test(resid_pmult, cluster = brood_grpnum, group = group_num, data = resid_dat, 
                alternative = "greater", method = "ds", exact = TRUE)



## Boxplot of residuals by bin
library(wesanderson)
mypal <- c(wes_palette("Darjeeling1")[c(2,4)])

## replace brood group labels with ordered ones
levels(resid_dat$brood_group)[levels(resid_dat$brood_group)=="<=2"] <- "<2.5"
resid_dat$brood_group2 <- ordered(resid_dat$brood_group, 
                                  levels = c("<2.5","3","4","5","6",">6.5"))

ggbp <- ggplot(resid_dat) + 
  geom_boxplot(aes(x = brood_group2, y = resid_pmult, 
                   group = interaction(brood_group2, factor(group_num)), 
                   fill = factor(group_num)), position=position_dodge(0.8)) +
  labs(x = "Brood/clutch/litter size", y = expression(paste(p[B], " - p"))) +
  scale_x_discrete(labels = c(expression("" <= 2.5), "[2.5, 3.5)", "[3.5, 4.5)", 
                              "[4.5, 5.5)", "[5.5, 6.5)", expression("" > 6.5))) +
  scale_fill_manual(name = "Taxa", values = mypal,
                     labels = c("Birds","Mammals")) + 
  theme_bw()

ggsave("boxplot_resids-birds-mammals.eps", plot = ggbp, device = cairo_ps, width = 12, height = 4, dpi = 320)


## rotated version of boxplot - might spread taxa better
ggbpr <- ggplot(resid_dat) + 
  geom_boxplot(aes(x = brood_group2, y = resid_pmult, 
                   group = interaction(brood_group2, factor(group_num)), 
                   fill = factor(group_num)), position=position_dodge(0.8)) +
  labs(x = "Brood/clutch/litter size", y = expression(paste(p[B], " - p"))) +
  scale_y_discrete(labels = c(expression("" <= 2.5), "[2.5, 3.5)", "[3.5, 4.5)", 
                              "[4.5, 5.5)", "[5.5, 6.5)", expression("" > 6.5))) +
  scale_fill_manual(name = "Taxa", values = mypal,
                    labels = c("Birds","Mammals")) + 
  theme_bw() + coord_flip()

ggsave("boxplot_resids-birds-mammals_rotated.eps", plot = ggbpr, device = cairo_ps, 
       width = 8, height = 12, dpi = 320)




#### Wilcoxon rank sum test for ALL four groups of species ####
## Residuals for birds (microsat-only)
load("/Users/Hannah/Documents/Research/Masamu/paternity/birds/analysis_microsat/MCMC_resids.rda")
bird_resids <- MCMC_resids
## Residuals for mammals (latest from phylo, 64spp)
load("/Users/Hannah/Documents/Research/Masamu/paternity/birds_vs_mammals/mammals_for_comparison/MCMC_estimates_byspecies/phyloMCMC_resids.rda")
mammal_resids <- MCMC_resids 
## Residuals for fish
load("/Users/Hannah/Documents/Research/Masamu/paternity/fish/single_run/MCMC_resids.rda")
fish_resids <- MCMC_resids 
## Residuals for reptiles/amphibians
load("/Users/Hannah/Documents/Research/Masamu/paternity/herps/single_run/MCMC_resids.rda")
herp_resids <- MCMC_resids 
## Residuals for inverts
load("/Users/Hannah/Documents/Research/Masamu/paternity/inverts/single_run/MCMC_resids.rda")
inverts_resids <- MCMC_resids 


mammal_resids$group <- "Mammals"
bird_resids$group <- "Birds"
fish_resids$group <- "Fish"
herp_resids$group <- "Herps"
inverts_resids$group <- "Inverts"

resid_dat_all <- rbind(mammal_resids, bird_resids, fish_resids, herp_resids, inverts_resids)

resid_dat_all$group_fac <- as.factor(resid_dat_all$group)
## "ds" method (i.e., Datta and Satten, informative clusters) does not behave with
## groups that are factors or characters; need numeric groups
resid_dat_all$group_num <- as.numeric(resid_dat_all$group_fac) 
# birds = 1, fish = 2, herps = 3, inverts = 4, mammals = 5


## subset to keep only broodsize < 2500 (fish/inverts) and convert to log scale
resid_dat_k <- resid_dat_all[resid_dat_all$avgbrood < 2500,]
resid_dat_k$avgbrood_log <- log(resid_dat_k$avgbrood)

## Consider many number of bins for k on the regular scale (use for birds, herps, mammals only)
## e.g., 0.75 <= bin("1") < 1.25; 1.25 <= bin("1.5") < 1.75; and so forth
resid_dat_k$brood_id1 <- ifelse((resid_dat_k$avgbrood - trunc(resid_dat_k$avgbrood)) < 0.25, 
                                floor(resid_dat_k$avgbrood), 
                                ifelse((resid_dat_k$avgbrood - trunc(resid_dat_k$avgbrood)) >= 0.25 &
                                         (resid_dat_k$avgbrood - trunc(resid_dat_k$avgbrood)) < 0.75, 
                                       floor(resid_dat_k$avgbrood) + 0.5, 
                                       ceiling(resid_dat_k$avgbrood)))

## Consider many number of bins for k on the log scale
## e.g., 0.75 <= bin("1") < 1.25; 1.25 <= bin("1.5") < 1.75; and so forth
resid_dat_k$brood_id2 <- ifelse((resid_dat_k$avgbrood_log - trunc(resid_dat_k$avgbrood_log)) < 0.25, 
                                floor(resid_dat_k$avgbrood_log), 
                                ifelse((resid_dat_k$avgbrood_log - trunc(resid_dat_k$avgbrood_log)) >= 0.25 &
                                         (resid_dat_k$avgbrood_log - trunc(resid_dat_k$avgbrood_log)) < 0.75, 
                                       floor(resid_dat_k$avgbrood_log) + 0.5, 
                                       ceiling(resid_dat_k$avgbrood_log)))


## "Are values of pB-p signif different across major taxa groups?"
## For smaller bins (larger number of bins @ 16)
clusWilcox.test(resid_pmult, cluster = brood_id2, group = group_num, data = resid_dat_k, 
                method = "ds")

## "Are values of pmult signif different across 3 major taxa groups (birds, mammals, herps)?"
## For smaller bins (larger number of bins @ 25)
clusWilcox.test(pmult, cluster = brood_id1, group = group_num, 
                data = resid_dat_k[resid_dat_k$group %in% c("Birds", "Mammals", "Herps"),], 
                method = "ds")


## Boxplot of residuals by bin
library(wesanderson)
mypal2 <- c(wes_palette("Darjeeling1")[c(5,2,4,1)], wes_palette("Darjeeling2")[5])
ggplot(resid_dat_k) + 
  geom_boxplot(aes(x = brood_id2, y = resid_pmult, 
                   group = interaction(brood_id2, group_num), color = factor(group_num))) +
  labs(x = "Log brood/clutch/litter size", y = expression(paste(p[B], " - p"))) +
  scale_color_manual(name = "Taxa", values = mypal2,
                     labels = c("Birds","Fish","Reptiles & amphibians","Invertebrates","Mammals")) + 
  theme_bw()


## Boxplot of observed pmult by bin
mypal2 <- c(wes_palette("Darjeeling1")[c(5,2,4,1)], wes_palette("Darjeeling2")[5])
ggplot(resid_dat_k) + 
  geom_boxplot(aes(x = brood_id2, y = pmult, 
                   group = interaction(brood_id2, group_num), color = factor(group_num))) +
  labs(x = "Log brood/clutch/litter size", y = "Observed probability of multiple paternity (p)") +
  scale_color_manual(name = "Taxa", values = mypal2,
                     labels = c("Birds","Fish","Reptiles & amphibians","Invertebrates","Mammals")) + 
  theme_bw()



## Plot of pB curve over p
ga <- ggplot(resid_dat_k) + 
  geom_point(aes(x = avgbrood_log, y = pmult, group = group_fac, color = group_fac), 
             alpha = 0.6, size = 2.5) +
  ## individual smooths per taxa
  geom_smooth(aes(x = avgbrood_log, y = pred_pmult, group = group_fac, color = group_fac), 
              se = FALSE, method = "loess") +
  # geom_ribbon(aes(x = avgbrood_log, ymin = loess(pred_pmult_lcl ~ avgbrood_log)$fitted,
  #                 ymax = loess(pred_pmult_ucl ~ avgbrood_log)$fitted,
  #                 group = group_fac, fill = group_fac), alpha = 0.4) +
  geom_smooth(aes(x = avgbrood_log, y = pred_pmult_ucl, group = group_fac, color = group_fac),
              se = FALSE, method = "loess", size = 0.5, alpha = 0.6) +
  geom_smooth(aes(x = avgbrood_log, y = pred_pmult_lcl, group = group_fac, color = group_fac),
              se = FALSE, method = "loess", size = 0.5, alpha = 0.6) +
  ## overall smooth for all taxa together
  geom_ribbon(aes(x = avgbrood_log, ymin = loess(pred_pmult_lcl ~ avgbrood_log)$fitted,
                  ymax = loess(pred_pmult_ucl ~ avgbrood_log)$fitted),
              color = "grey20", fill = "grey60", alpha = 0.4) +
  geom_smooth(aes(x = avgbrood_log, y = pred_pmult), color = "grey60", linetype = "dashed",
              se = FALSE, method = "loess", alpha = 1) +
  # ## overall smooth for all taxa together
  # geom_smooth(aes(x = avgbrood_log, y = pred_pmult), color = "blue", linetype = "dashed",
  #             se = FALSE, method = "loess") +
  # geom_smooth(aes(x = avgbrood_log, y = pred_pmult_ucl), color = "blue", linetype = "dashed", 
  #             se = FALSE, method = "loess", size = 0.5, alpha = 0.6) +
  # geom_smooth(aes(x = avgbrood_log, y = pred_pmult_lcl), color = "blue", linetype = "dashed", 
  #             se = FALSE, method = "loess", size = 0.5, alpha = 0.6) +
  labs(x = "Log brood/clutch/litter size", y = "Probability of multiple paternity") +
  scale_color_manual(name = "Taxa", values = mypal2,
                     labels = c("Birds","Fish","Reptiles & amphibians","Invertebrates","Mammals")) + 
  theme_bw()

ggsave("pB_all-taxa.eps", ga, device = cairo_ps, width = 10, height = 6.5, dpi = 320)



#### Plot of pB-p versus avgbrood (predict with loess) - not helpful
resid_dat2 <- resid_dat
resid_dat2$pred_resid_pmult <- NA
resid_dat2$pred_resid_pmult_lcl <- NA
resid_dat2$pred_resid_pmult_ucl <- NA

# mammals - preds
loess_resid <- loess(resid_pmult ~ avgbrood, data = resid_dat2[resid_dat2$group=="Mammals",])
resid_dat2[resid_dat2$group=="Mammals",]$pred_resid_pmult <- 
  predict(loess_resid, newdata = resid_dat2[resid_dat2$group=="Mammals",]$avgbrood)
# mammals - lower limit
loess_resid_pmult_lcl <- loess(resid_pmult_lcl ~ avgbrood, data = resid_dat2[resid_dat2$group=="Mammals",])
resid_dat2[resid_dat2$group=="Mammals",]$pred_resid_pmult_lcl <- 
  predict(loess_resid_pmult_lcl, newdata = resid_dat2[resid_dat2$group=="Mammals",]$avgbrood)
# mammals - upper limit
loess_resid_pmult_ucl <- loess(resid_pmult_ucl ~ avgbrood, data = resid_dat2[resid_dat2$group=="Mammals",])
resid_dat2[resid_dat2$group=="Mammals",]$pred_resid_pmult_ucl <- 
  predict(loess_resid_pmult_ucl, newdata = resid_dat2[resid_dat2$group=="Mammals",]$avgbrood)

# birds - preds
loess_resid <- loess(resid_pmult ~ avgbrood, data = resid_dat2[resid_dat2$group=="Birds",])
resid_dat2[resid_dat2$group=="Birds",]$pred_resid_pmult <- 
  predict(loess_resid, newdata = resid_dat2[resid_dat2$group=="Birds",]$avgbrood)
# birds - lower limit
loess_resid_pmult_lcl <- loess(resid_pmult_lcl ~ avgbrood, data = resid_dat2[resid_dat2$group=="Birds",])
resid_dat2[resid_dat2$group=="Birds",]$pred_resid_pmult_lcl <- 
  predict(loess_resid_pmult_lcl, newdata = resid_dat2[resid_dat2$group=="Birds",]$avgbrood)
# birds - upper limit
loess_resid_pmult_ucl <- loess(resid_pmult_ucl ~ avgbrood, data = resid_dat2[resid_dat2$group=="Birds",])
resid_dat2[resid_dat2$group=="Birds",]$pred_resid_pmult_ucl <- 
  predict(loess_resid_pmult_ucl, newdata = resid_dat2[resid_dat2$group=="Birds",]$avgbrood)


# stat_smooth with method = "loess" and default span
gpl1 <- ggplot() +
  # geom_point(data = dat, aes(avgbrood, pmult), alpha = 0.5, size = 2.5) +
  # stat_smooth(data = dat, aes(avgbrood, pmult), color = "black", se = FALSE, method = "loess") +
  stat_smooth(data = resid_dat2, aes(avgbrood, resid_pmult, color = group), 
              se = FALSE, method = "loess") +
  geom_point(data = resid_dat2, aes(avgbrood, resid_pmult, color = group), 
             alpha = 0.5, size = 2.5) +
  stat_smooth(data = resid_dat2, aes(avgbrood, resid_pmult_ucl, color = group), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "loess") +
  stat_smooth(data = resid_dat2, aes(avgbrood, resid_pmult_lcl, color = group), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "loess") +
  labs(x="Brood/litter size", y=expression(p[B]-p)) +
  lims(y = c(-1,1)) +
  scale_x_continuous(breaks = c(seq(1,12, by = 1)) , limits = c(1,12)) +
  scale_color_manual(name = "",  #labels = c(),
                     values = c("#1B9E77", "#D95F02")) +
  theme_bw(base_size = 16)

gpl1 


# predictions from loess(), then plotted directly with stat_smooth
gpl2 <- ggplot() +
  # geom_point(data = dat, aes(avgbrood, pmult), alpha = 0.5, size = 2.5) +
  # stat_smooth(data = dat, aes(avgbrood, pmult), color = "black", se = FALSE, method = "loess") +
  stat_smooth(data = resid_dat2, aes(avgbrood, pred_resid_pmult, color = group), 
              se = FALSE, method = "loess") +
  geom_point(data = resid_dat2, aes(avgbrood, resid_pmult, color = group), 
             alpha = 0.5, size = 2.5) +
  stat_smooth(data = resid_dat2, aes(avgbrood, pred_resid_pmult_ucl, color = group), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "loess") +
  stat_smooth(data = resid_dat2, aes(avgbrood, pred_resid_pmult_lcl, color = group), 
              alpha = 0.5, se = FALSE, size = 0.5, method = "loess") +
  labs(x="Brood/litter size", y=expression(p[B]-p)) +
  lims(y = c(-1,1)) +
  scale_x_continuous(breaks = c(seq(1,12, by = 1)) , limits = c(1,12)) +
  scale_color_manual(name = "",  #labels = c(),
                     values = c("#1B9E77", "#D95F02")) +
  theme_bw(base_size = 16)

gpl2 

