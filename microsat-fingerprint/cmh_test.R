### Cochran-Mantel-Haenszel test for DNA methods
library(tidyverse)
library(vcd)
library(readxl)
# paternity_birds <- read_excel("~/Documents/Research/Masamu/paternity/birds/paternity_birds.xlsx", 
#                               col_types = c("skip", "text", "text", 
#                                             "text", "text", "numeric", "text", 
#                                             "text", "text", "text", "text", "text", 
#                                             "numeric", "text", "text", "numeric", 
#                                             "numeric", "numeric", "numeric", 
#                                             "numeric", "numeric", "numeric", 
#                                             "text", "text", "text", "text", "text", 
#                                             "text", "numeric"))
paternity_birds <- read.csv("~/Documents/Research/Masamu/paternity/birds/paternity_birds.csv")
# View(paternity_birds)

## Create category for MP
paternity_birds$MP_cat <- ifelse(paternity_birds$pmult==0, "azero", "positive")

## Create category for brood size 
## (half-integer and integer alone create stratums with sample size < 2)
## try "small", "med", "large" instead
paternity_birds$brood_id <- ifelse(paternity_birds$avgbrood <= 3, "small",
                                   ifelse(paternity_birds$avgbrood <= 6, "medium", "large"))

## What are the sizes of each of these groups?
paternity_birds %>% filter(!is.na(pmult)) %>% group_by(brood_id) %>% summarise(count = n())


## Three-dimensional table (as list)
tab3 <- with(paternity_birds[!is.na(paternity_birds$pmult),], table(DNA, MP_cat, brood_id))

## Three-dimensional table (with xtabs)
ptab <- paternity_birds %>%
  filter(!is.na(pmult)) %>%
  count(DNA, MP_cat, brood_id) %>%
  group_by(DNA) 

ptab3 <- xtabs(n ~ DNA + MP_cat + brood_id, ptab) ## same as tab3, but going to use this one

## CMH test
mantelhaen.test(ptab3)



### Could try smaller groupings
## try integer with a "6.5+" group and a "1 to 2.5" group instead 
paternity_birds$brood_id2 <- ifelse((paternity_birds$avgbrood - trunc(paternity_birds$avgbrood)) < 0.5, 
                                    floor(paternity_birds$avgbrood), ceiling(paternity_birds$avgbrood))
paternity_birds$brood_group <- ifelse(paternity_birds$brood_id2 >= 7, ">6.5",
                                      ifelse(paternity_birds$brood_id2 <=2, "<2.5",
                                             as.character(paternity_birds$brood_id2)))
paternity_birds$brood_group <- ordered(paternity_birds$brood_group, 
                                       levels = c("<2.5","3","4","5","6",">6.5"))

## What are the sizes of each of these groups?
paternity_birds %>% filter(!is.na(pmult)) %>% group_by(brood_group) %>% summarise(count = n())

## Three-dimensional table 
ptab <- paternity_birds %>%
  filter(!is.na(pmult)) %>%
  count(DNA, MP_cat, brood_group) %>%
  group_by(DNA) 

ptab4 <- xtabs(n ~ DNA + MP_cat + brood_group, ptab) 

## CMH test
mantelhaen.test(ptab4, exact = TRUE, correct = FALSE)

## Marginal odds ratios
AC <- margin.table(ptab4, c(2,1))

oddsratio(ptab4, 3, log = FALSE)

lor <- oddsratio(ptab4,3)
exp(confint(lor)) ## CI
summary(lor)
plot(lor, xlab="Brood size", ylab = "LOR(DNA / MP)", main = "")


# ## or half-integer with a "7.5+" group
# paternity_birds$brood_id3 <- ifelse((paternity_birds$avgbrood - trunc(paternity_birds$avgbrood)) < 0.25,
#                                     floor(paternity_birds$avgbrood),
#                                     ifelse((paternity_birds$avgbrood - trunc(paternity_birds$avgbrood)) >= 0.25 &
#                                              (paternity_birds$avgbrood - trunc(paternity_birds$avgbrood)) < 0.75,
#                                            floor(paternity_birds$avgbrood) + 0.5,
#                                            ceiling(paternity_birds$avgbrood)))
# paternity_birds$brood_group3 <- ifelse(paternity_birds$brood_id3 >= 7.5, "7.5+",
#                                        as.character(paternity_birds$brood_id3))
# ## What are the sizes of each of these groups?
# paternity_birds %>% filter(!is.na(pmult)) %>% group_by(brood_group3) %>% summarise(count = n())
# 
# ## Three-dimensional table 
# ptab <- paternity_birds %>%
#   filter(!is.na(pmult)) %>%
#   count(DNA, MP_cat, brood_group3) %>%
#   group_by(DNA) 
# 
# ptab5 <- xtabs(n ~ DNA + MP_cat + brood_group3, ptab) 
# 
# ## CMH test
# mantelhaen.test(ptab5, exact = TRUE, correct = FALSE)




#### Do the same for mammals vs birds
## Residuals for microsat-only
load("/Users/Hannah/Documents/Research/Masamu/paternity/birds/analysis_microsat/MCMC_resids.rda")
bird_resids <- MCMC_resids
## Residuals for mammals (latest from phylo, 63spp)
load("/Users/Hannah/Documents/Research/Masamu/paternity/phylogeny/MCMC_estimates_byspecies/phyloMCMC_resids.rda")
mammal_resids <- MCMC_resids 

mammal_resids$group <- "Mammals"
bird_resids$group <- "Birds"

resid_dat <- rbind(mammal_resids, bird_resids)

## Create category for MP
resid_dat$MP_cat <- ifelse(resid_dat$pmult==0, "azero", "positive")

# ## half-integer with a "7.5+" group
# resid_dat$brood_id3 <- ifelse((resid_dat$avgbrood - trunc(resid_dat$avgbrood)) < 0.25,
#                               floor(resid_dat$avgbrood),
#                               ifelse((resid_dat$avgbrood - trunc(resid_dat$avgbrood)) >= 0.25 &
#                                        (resid_dat$avgbrood - trunc(resid_dat$avgbrood)) < 0.75,
#                                      floor(resid_dat$avgbrood) + 0.5,
#                                      ceiling(resid_dat$avgbrood)))
# resid_dat$brood_group3 <- ifelse(resid_dat$brood_id3 >= 7.5, "7.5+",
#                                  as.character(resid_dat$brood_id3))

## whole-integer with a "7+" group
resid_dat$brood_id2 <- ifelse((resid_dat$avgbrood - trunc(resid_dat$avgbrood)) < 0.5,
                              floor(resid_dat$avgbrood), ceiling(resid_dat$avgbrood))
resid_dat$brood_group <- ifelse(resid_dat$brood_id2 >= 7, ">6.5",
                                ifelse(resid_dat$brood_id2 <=2, "<2.5",
                                       as.character(resid_dat$brood_id2)))
resid_dat$brood_group <- ordered(resid_dat$brood_group, 
                                 levels = c("<2.5","3","4","5","6",">6.5"))

## What are the sizes of each of these groups?
resid_dat %>% filter(!is.na(pmult)) %>% group_by(brood_group) %>% summarise(count = n())

## Three-dimensional table 
ptab <- resid_dat %>%
  filter(!is.na(pmult)) %>%
  count(group, MP_cat, brood_group) %>%
  group_by(group) 

ptab6 <- xtabs(n ~ group + MP_cat + brood_group, ptab) 

## CMH test
mantelhaen.test(ptab6, exact = TRUE, correct = FALSE)

## Marginal odds ratios
AC <- margin.table(ptab6, c(2,1))

oddsratio(ptab6, 3, log = FALSE)

lor <- oddsratio(ptab6,3)
exp(confint(lor)) ## CI
summary(lor)
plot(lor, xlab="Brood/litter size", ylab = "LOR(taxa / MP)", main = "")



#### Panel plot of conditional LOR for DNA methods and birds vs mammals comparisons ####
library(gridExtra)
p1 <- plot(oddsratio(ptab4,3), xlab="Brood size", ylab = "LOR", main = "", return_grob = TRUE)
p2 <- plot(oddsratio(ptab6,3), xlab="Brood/litter size", ylab = "LOR", main = "", return_grob = TRUE)

gpl <- grid.arrange(p1, p2, newpage = TRUE)

ggsave("conditionalLOR.eps", plot = gpl, device = cairo_ps, 
       dpi = 320, width = 7, height = 10, units = "in")

## alternative panel plot with labels
library(cowplot)
gpl2 <- plot_grid(p1, p2, labels = c("(a)","(b)"), ncol = 1)
ggsave("conditionalLOR_labelled.eps", plot = gpl2, device = cairo_ps, 
       dpi = 320, width = 7, height = 10, units = "in")

## save individual plots
ggsave("conditionalLOR_DNA.eps", plot = p1, device = cairo_ps, 
       dpi = 320, width = 7, height = 5, units = "in")
ggsave("conditionalLOR_BvM.eps", plot = p2, device = cairo_ps, 
       dpi = 320, width = 7, height = 5, units = "in")


#### ggplot version because Ash wants the common odds on the plot too ####
library(ggplot2)
library(wesanderson)

lor_dna <- oddsratio(ptab4, 3)
lor_dna_ci <- confint(lor_dna)
mtest_dna <- mantelhaen.test(ptab4, exact = TRUE, correct = FALSE)
lor_bvm <- oddsratio(ptab6, 3)
lor_bvm_ci <- confint(lor_bvm)
mtest_bvm <- mantelhaen.test(ptab6, exact = TRUE, correct = FALSE)

lor_dna_dat <- data.frame(brood_group = ordered(names(lor_dna$coefficients),
                                               levels = c("<2.5","3","4","5","6",">6.5")),
                          xAxis = 1:length(names(lor_dna$coefficients)),
                          coeff = unname(lor_dna$coefficients),
                          lcl = unname(lor_dna_ci[,1]),
                          ucl = unname(lor_dna_ci[,2]))

lor_bvm_dat <- data.frame(brood_group = ordered(names(lor_bvm$coefficients),
                                                levels = c("<2.5","3","4","5","6",">6.5")),
                          xAxis = 1:length(names(lor_bvm$coefficients)),
                          coeff = unname(lor_bvm$coefficients),
                          lcl = unname(lor_bvm_ci[,1]),
                          ucl = unname(lor_bvm_ci[,2]))

## plot of DNA comparison conditional odds ratios
lor_dna_p <- ggplot(lor_dna_dat) +
  geom_col(aes(x = brood_group, y = coeff), fill = "grey90", color = "grey70") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey20") +
  ## common LOR error band for all brood groups as colored block (not pretty)
  # geom_ribbon(aes(ymin = log(mtest_dna$conf.int[1]), ymax = log(mtest_dna$conf.int[2]),
  #                 x = xAxis), fill = wes_palette("Darjeeling1")[2], alpha = 0.3) +
  ## common LOR estimate
  geom_hline(yintercept = log(unname(mtest_dna$estimate)), linetype = "longdash", size = 1.2) +
  ## common LOR error band for all brood groups
  geom_hline(yintercept = log(mtest_dna$conf.int[1]), linetype = "solid",
             color = wes_palette("Darjeeling1")[2], size = 1.2) +
  geom_hline(yintercept = log(mtest_dna$conf.int[2]), linetype = "solid",
             color = wes_palette("Darjeeling1")[2], size = 1.2) +
  geom_point(aes(x = brood_group, y = coeff, group = 1), 
             color = wes_palette("FantasticFox1")[5], size = 4) +
  geom_line(aes(x = brood_group, y = coeff, group = 1), 
            color = wes_palette("FantasticFox1")[5], size = 1) +
  ## conditional LOR error bars for each brood group
  # geom_errorbar(aes(x = brood_group, ymax = ucl, ymin = lcl), 
  #               color = wes_palette("FantasticFox1")[5],
  #               size = 1, width = 0.0) +
  # geom_dl(aes(x = 14.2, y = (log(unname(mtest_dna$estimate)) + 0.5), label = "CLOR"), 
  #         method = list(dl.trans(x = x + .2), "last.points")) +
  labs(x = "Brood size", y = "Log odds ratio") +
  scale_y_continuous(breaks = seq(-2, 4, 1)) +
  scale_x_discrete(labels = c(expression("" <= 2.5), "[2.5, 3.5)", "[3.5, 4.5)", 
                                          "[4.5, 5.5)", "[5.5, 6.5)", expression("" > 6.5))) +
  # scale_x_continuous(limits = c(min(lor_dna_dat$xAxis) - 0.5, max(lor_dna_dat$xAxis) + 1),
  #                    breaks = (min(lor_dna_dat$xAxis)):(max(lor_dna_dat$xAxis) + 1), 
  #                    labels = c(lor_dna_dat$brood_group,'')) +
  annotate("text", x = 7, y = log(unname(mtest_dna$estimate)), hjust = -0.1,
           label = paste0("CLOR = ", round(log(unname(mtest_dna$estimate)), 2))) +
  theme_classic() + theme(panel.border = element_rect(color='grey50', fill = NA),
                          text = element_text(size = 12),
                          plot.margin = unit(c(1, 6, 1, 1), "lines")) + 
  coord_cartesian(clip = "off")

ggsave("conditionalLOR_DNA_gg.eps", plot = lor_dna_p, device = cairo_ps, 
       dpi = 320, width = 8, height = 4, units = "in")
# ggsave("conditionalLOR_DNA_ggalt.eps", plot = lor_dna_p, device = cairo_ps, 
#        dpi = 320, width = 7, height = 4, units = "in")


## plot of birds vs mammals comparison conditional odds ratios
lor_bvm_p <- ggplot(lor_bvm_dat) +
  geom_col(aes(x = brood_group, y = coeff), fill = "grey90", color = "grey70") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey20") +
  ## common LOR estimate
  geom_hline(yintercept = log(unname(mtest_bvm$estimate)), linetype = "longdash", size = 1.2) +
  ## common LOR error band for all brood groups
  geom_hline(yintercept = log(mtest_bvm$conf.int[1]), linetype = "solid",
             color = wes_palette("Darjeeling1")[2], size = 1.2) +
  geom_hline(yintercept = log(mtest_bvm$conf.int[2]), linetype = "solid",
             color = wes_palette("Darjeeling1")[2], size = 1.2) +
  geom_point(aes(x = brood_group, y = coeff, group = 1), 
             color = wes_palette("FantasticFox1")[5], size = 4) +
  geom_line(aes(x = brood_group, y = coeff, group = 1), 
            color = wes_palette("FantasticFox1")[5], size = 1) +
  ## conditional LOR error bars for each brood group
  # geom_errorbar(aes(x = brood_group, ymax = ucl, ymin = lcl), 
  #               color = wes_palette("FantasticFox1")[5],
  #               size = 1, width = 0.0) +
  # geom_dl(aes(x = 14.2, y = (log(unname(mtest_bvm$estimate)) + 0.5), label = "CLOR"), 
  #         method = list(dl.trans(x = x + .2), "last.points")) +
  labs(x = "Brood/litter size", y = "Log odds ratio") +
  scale_y_continuous(breaks = seq(-2, 4, 1)) +
  scale_x_discrete(labels = c(expression("" <= 2.5), "[2.5, 3.5)", "[3.5, 4.5)", 
                              "[4.5, 5.5)", "[5.5, 6.5)", expression("" > 6.5))) +
  # scale_x_continuous(limits = c(min(lor_bvm_dat$xAxis) - 0.5, max(lor_bvm_dat$xAxis) + 1),
  #                    breaks = (min(lor_bvm_dat$xAxis)):(max(lor_bvm_dat$xAxis) + 1), 
  #                    labels = c(lor_bvm_dat$brood_group,'')) +
  annotate("text", x = 7, y = log(unname(mtest_bvm$estimate)), hjust = -0.1,
           label = paste0("CLOR = ", format(round(log(unname(mtest_bvm$estimate)), 2), nsmall = 2))) +
  theme_classic() + theme(panel.border = element_rect(color='grey50', fill = NA),
                          text = element_text(size = 12),
                          plot.margin = unit(c(1, 6, 1, 1), "lines")) + 
  coord_cartesian(clip = "off")

ggsave("conditionalLOR_BvM_gg.eps", plot = lor_bvm_p, device = cairo_ps, 
       dpi = 320, width = 8, height = 5, units = "in")
# ggsave("conditionalLOR_BvM_ggalt.eps", plot = lor_bvm_p, device = cairo_ps, 
#        dpi = 320, width = 7, height = 5, units = "in")


## save both together in one figure
## alternative panel plot with labels
library(cowplot)
lor_gg <- plot_grid(lor_dna_p, lor_bvm_p, labels = c("(a)","(b)"), ncol = 1)
ggsave("conditionalLOR_labelled.eps", plot = lor_gg, device = cairo_ps, 
       dpi = 320, width = 8, height = 10, units = "in")


#### Alternative, more forest-plot like (not keen) ####
lor_dna_dat2 <- data.frame(brood_group3 = c(names(lor_dna$coefficients), "Common"),
                           yAxis = (length(names(lor_dna$coefficients)) + 1):1,
                           coeff = c(unname(lor_dna$coefficients), log(unname(mtest_dna$estimate))),
                           lcl = c(unname(lor_bvm_ci[,1]), log(mtest_dna$conf.int[1])),
                           ucl = c(unname(lor_bvm_ci[,2]), log(mtest_dna$conf.int[2])))

# Plot
p <- ggplot(lor_dna_dat2, aes(x = coeff, y = yAxis)) +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ucl, xmin = lcl), 
                 size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = wes_palette("FantasticFox1")[5]) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = lor_dna_dat2$yAxis, labels = lor_dna_dat2$brood_group3) +
  labs(x = "Log odds ratio", y = "Brood size") +
  scale_x_continuous(breaks = seq(-4, 8, 2)) # +
# annotate(geom = "text", y =1.1, x = 3.5, label ="Model p < 0.001\nPseudo R^2 = 0.10", size = 3.5, hjust = 0) + 
# ggtitle("Intention to remove box turtles from the road")

