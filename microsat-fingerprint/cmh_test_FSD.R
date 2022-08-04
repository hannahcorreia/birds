### Cochran-Mantel-Haenszel test for DNA methods
library(tidyverse)
library(vcd)
library(readxl)
paternity_birds <- read.csv("paternity_birds.csv", header = TRUE)
# View(paternity_birds)

## Create category for pmult
paternity_birds$pmult_cat <- ifelse(paternity_birds$pmult==0, "azero", "positive")

## Create category for brood size 
## half-integer with a "7.5+" group
paternity_birds$brood_id3 <- ifelse((paternity_birds$avgbrood - trunc(paternity_birds$avgbrood)) < 0.25,
                                    floor(paternity_birds$avgbrood),
                                    ifelse((paternity_birds$avgbrood - trunc(paternity_birds$avgbrood)) >= 0.25 &
                                             (paternity_birds$avgbrood - trunc(paternity_birds$avgbrood)) < 0.75,
                                           floor(paternity_birds$avgbrood) + 0.5,
                                           ceiling(paternity_birds$avgbrood)))
paternity_birds$brood_group3 <- ifelse(paternity_birds$brood_id3 >= 7.5, "7.5+",
                                       as.character(paternity_birds$brood_id3))
## What are the sizes of each of these groups?
paternity_birds %>% filter(!is.na(pmult)) %>% group_by(brood_group3) %>% summarise(count = n())

## Three-dimensional table 
ptab <- paternity_birds %>%
  filter(!is.na(pmult)) %>%
  count(DNA, pmult_cat, brood_group3) %>%
  group_by(DNA) 

ptab5 <- xtabs(n ~ DNA + pmult_cat + brood_group3, ptab) 

## CMH test
mantelhaen.test(ptab5, exact = TRUE, correct = FALSE)

## Conditional odds ratios
oddsratio(ptab5, 3, log = FALSE)

## Conditional log odds ratio plot
lor <- oddsratio(ptab5,3)
exp(confint(lor)) ## CI
summary(lor)
plot(lor, xlab="Brood size", ylab = "LOR(DNA / pmult)", main = "")


#### Do the same for mammals vs birds
## Residuals for microsat-only
load("/Users/Hannah/Dropbox/Paternity_Birds/MCMC_birds.rda")
bird_resids <- MCMC_resids
## Residuals for mammals (latest from phylo, 63spp)
load("/Users/Hannah/Dropbox/Paternity_Birds/MCMC_mammals.rda")
mammal_resids <- MCMC_resids 

mammal_resids$group <- "Mammals"
bird_resids$group <- "Birds"

resid_dat <- rbind(mammal_resids, bird_resids)

## Create category for pmult
resid_dat$pmult_cat <- ifelse(resid_dat$pmult==0, "azero", "positive")

## half-integer with a "7.5+" group
resid_dat$brood_id3 <- ifelse((resid_dat$avgbrood - trunc(resid_dat$avgbrood)) < 0.25,
                              floor(resid_dat$avgbrood),
                              ifelse((resid_dat$avgbrood - trunc(resid_dat$avgbrood)) >= 0.25 &
                                       (resid_dat$avgbrood - trunc(resid_dat$avgbrood)) < 0.75,
                                     floor(resid_dat$avgbrood) + 0.5,
                                     ceiling(resid_dat$avgbrood)))
resid_dat$brood_group3 <- ifelse(resid_dat$brood_id3 >= 7.5, "7.5+",
                                 as.character(resid_dat$brood_id3))
## What are the sizes of each of these groups?
resid_dat %>% filter(!is.na(pmult)) %>% group_by(brood_group3) %>% summarise(count = n())

## Three-dimensional table 
ptab <- resid_dat %>%
  filter(!is.na(pmult)) %>%
  count(group, pmult_cat, brood_group3) %>%
  group_by(group) 

ptab6 <- xtabs(n ~ group + pmult_cat + brood_group3, ptab) 

## CMH test
mantelhaen.test(ptab6, exact = TRUE, correct = FALSE)

## Marginal odds ratios
AC <- margin.table(ptab6, c(2,1))

oddsratio(ptab6, 3, log = FALSE)

lor <- oddsratio(ptab6,3)
exp(confint(lor)) ## CI
summary(lor)
plot(lor, xlab="Brood/litter size", ylab = "LOR(taxa / pmult)", main = "")


