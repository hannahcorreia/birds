### Cochran-Mantel-Haenszel test for socially monogamous vs not
library(tidyverse)
library(vcd)
# library(readxl)
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
## try integer with a "7+" group  
paternity_birds$brood_id <- ifelse((paternity_birds$avgbrood - trunc(paternity_birds$avgbrood)) < 0.5, 
                                   floor(paternity_birds$avgbrood), ceiling(paternity_birds$avgbrood))
paternity_birds$brood_group <- ifelse(paternity_birds$brood_id >= 7, ">6.5",
                                      ifelse(paternity_birds$brood_id <=2, "<2.5",
                                             as.character(paternity_birds$brood_id)))
paternity_birds$brood_group <- ordered(paternity_birds$brood_group, 
                                       levels = c("<2.5","3","4","5","6",">6.5"))

## What are the sizes of each of these groups?
paternity_birds %>% filter(!is.na(pmult)) %>% group_by(brood_group) %>% summarise(count = n())


## SocMon is currently a 3-category factor - change to 2-category (mono/xmono only)
paternity_birds$SocMon_bin <- ifelse(paternity_birds$SocMon %in% c("mono","xmono"), 
                                     paternity_birds$SocMon, NA)


## Three-dimensional table 
ptab <- paternity_birds %>%
  filter(!is.na(pmult)) %>%
  count(SocMon_bin, MP_cat, brood_group) %>%
  group_by(SocMon_bin) 

ptab4 <- xtabs(n ~ SocMon_bin + MP_cat + brood_group, ptab) 

## CMH test
mantelhaen.test(ptab4, exact = TRUE, correct = FALSE)

## Marginal odds ratios
AC <- margin.table(ptab4, c(2,1))

oddsratio(ptab4, 3, log = FALSE)

lor <- oddsratio(ptab4,3)
exp(confint(lor)) ## CI
summary(lor)
plot(lor, xlab="Brood size", ylab = "LOR(SocMono / MP)", main = "")





#### ggplot version of panel plot of conditional LOR for mono vs not comparison ####
library(ggplot2)
library(wesanderson)
lor_mono <- oddsratio(ptab4, 3)
lor_mono_ci <- confint(lor_mono)
mtest_mono <- mantelhaen.test(ptab4, exact = TRUE, correct = FALSE)


lor_mono_dat <- data.frame(brood_group = ordered(names(lor_mono$coefficients),
                                                 levels = c("<2.5","3","4","5","6",">6.5")),
                           xAxis = 1:length(names(lor_mono$coefficients)),
                           coeff = unname(lor_mono$coefficients),
                           lcl = unname(lor_mono_ci[,1]),
                           ucl = unname(lor_mono_ci[,2]))


## plot of socio mono comparison conditional odds ratios
lor_mono_p <- ggplot(lor_mono_dat) +
  geom_col(aes(x = as.factor(brood_group), y = coeff), fill = "grey90", color = "grey70") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey20") +
  ## common LOR error band for all brood groups as colored block (not pretty)
  # geom_ribbon(aes(ymin = log(mtest_mono$conf.int[1]), ymax = log(mtest_mono$conf.int[2]),
  #                 x = xAxis), fill = wes_palette("Darjeeling1")[2], alpha = 0.3) +
  ## common LOR estimate
  geom_hline(yintercept = log(unname(mtest_mono$estimate)), linetype = "longdash", size = 1.2) +
  ## common LOR error band for all brood groups
  geom_hline(yintercept = log(mtest_mono$conf.int[1]), linetype = "solid",
             color = wes_palette("Darjeeling1")[2], size = 1.2) +
  geom_hline(yintercept = log(mtest_mono$conf.int[2]), linetype = "solid",
             color = wes_palette("Darjeeling1")[2], size = 1.2) +
  geom_point(aes(x = as.factor(brood_group), y = coeff, group = 1), 
             color = wes_palette("FantasticFox1")[5], size = 4) +
  geom_line(aes(x = as.factor(brood_group), y = coeff, group = 1), 
            color = wes_palette("FantasticFox1")[5], size = 1) +
  ## conditional LOR error bars for each brood group
  # geom_errorbar(aes(x = as.factor(brood_group), ymax = ucl, ymin = lcl), 
  #               color = wes_palette("FantasticFox1")[5],
  #               size = 1, width = 0.0) +
  # geom_dl(aes(x = 14.2, y = (log(unname(mtest_mono$estimate)) + 0.5), label = "CLOR"), 
  #         method = list(dl.trans(x = x + .2), "last.points")) +
  labs(x = "Brood size", y = "Log odds ratio") +
  scale_y_continuous(breaks = seq(-2, 4, 1)) +
  scale_x_discrete(labels = c(expression("" <= 2.5), "[2.5, 3.5)", "[3.5, 4.5)", 
                              "[4.5, 5.5)", "[5.5, 6.5)", expression("" > 6.5))) +
  # scale_x_continuous(limits = c(min(lor_mono_dat$xAxis) - 0.5, max(lor_mono_dat$xAxis) + 1),
  #                    breaks = (min(lor_mono_dat$xAxis)):(max(lor_mono_dat$xAxis) + 1), 
  #                    labels = c(lor_mono_dat$brood_group,'')) +
  annotate("text", x = 7, y = log(unname(mtest_mono$estimate)), hjust = -0.1,
           label = paste0("CLOR = ", round(log(unname(mtest_mono$estimate)), 2))) +
  theme_classic() + theme(panel.border = element_rect(color='grey50', fill = NA),
                          text = element_text(size = 12),
                          plot.margin = unit(c(1, 6, 1, 1), "lines")) + 
  coord_cartesian(clip = "off")

ggsave("conditionalLOR_mono_gg.eps", plot = lor_mono_p, device = cairo_ps, 
       dpi = 320, width = 8, height = 4, units = "in")
# ggsave("conditionalLOR_mono_ggalt.eps", plot = lor_mono_p, device = cairo_ps, 
#        dpi = 320, width = 7, height = 4, units = "in")




