#### Clean and organize bird data for MCMC analyses ####
#### Last run 16 June 2021

library(readxl)
library(tidyverse)
setwd("~/Documents/Research/Masamu/paternity/birds")

## Current birds MP file
# bird_MP <- read_excel("BirdData16Jun21_final.xlsx")
bird_MP <- read_excel("BirdData4Mar22_final.xlsx")[,-1]
#View(bird_MP) 


# ## Quick merge of microsat indicator variable with current version of birds MP data
# 
# ## Microsat indicator file
# bird_fingerprint_MS <- read_excel("Bird Data - MP records - 20 May 2021.xlsx")
# #View(bird_fingerprint_MS)   
# 
# ## Keep only ID columns for merging and Fingerprint/Microsat
# bird_fingerprint_MS2 <- bird_fingerprint_MS[,c("CommonName","Scientific","SpeciesID",
#                                                "StudyID","PopulationID","Fingerprint(1)/Microsat(2)")]
# 
# bird_full <- merge(bird_MP, bird_fingerprint_MS2, 
#       by = c("CommonName","Scientific","SpeciesID","StudyID","PopulationID"))
# 
# colnames(bird_full)[colnames(bird_full)=="Fingerprint(1)/Microsat(2)"] <- "DNA_technique"
# 
# all.equal(bird_full$DNA, bird_full$DNA_technique)
# ## These are the same variables, so no need to merge after all!


bird_MP$species <- paste0(bird_MP$CommonName_cap, ", ", bird_MP$Scientific)

bird_MP2 <- bird_MP %>% relocate(species, .after = Scientific) %>%
  # 1 = fingerprint; 2 = microsatellite
  mutate(DNA = factor(DNA),
         # 1 = socially monogamous; 2 = NOT socially monogamous (remove "8" = no info)
         SocMon = ifelse(`Soc Monog`==1, 1, ifelse(`Soc Monog`==8, 99, 2)),
         SocMon = fct_recode(factor(SocMon), mono = "1", xmono = "2", no_info = "99"))

head(bird_MP2)

library(xlsx)
write.xlsx(bird_MP2, file = "paternity_birds.xlsx")

paternity_birds <- bird_MP2 %>%
  rename(nbrood = Nbr, avgbrood = MeanBS, minbrood = MinBS, maxbrood = MaxBS,
         pmult = MP, avgsire = MeanS, minsire = MinS, maxsire = MaxS) %>%
  select(species, nbrood, avgbrood, minbrood, maxbrood, pmult, avgsire, minsire, maxsire,
         DNA, SocMon, Reference)

write.csv(paternity_birds, file = "paternity_birds.csv", row.names = FALSE)


## Microsat-only data
paternity_birds_ms <- bird_MP2 %>%
  rename(nbrood = Nbr, avgbrood = MeanBS, minbrood = MinBS, maxbrood = MaxBS,
         pmult = MP, avgsire = MeanS, minsire = MinS, maxsire = MaxS) %>%
  filter(DNA == "microsatellite") %>%
  select(species, nbrood, avgbrood, minbrood, maxbrood, pmult, avgsire, minsire, maxsire,
         Reference)

write.csv(paternity_birds_ms, file = "paternity_birds_ms.csv", row.names = FALSE)


## Fingerprint-only data
paternity_birds_fp <- bird_MP2 %>%
  rename(nbrood = Nbr, avgbrood = MeanBS, minbrood = MinBS, maxbrood = MaxBS,
         pmult = MP, avgsire = MeanS, minsire = MinS, maxsire = MaxS) %>%
  filter(DNA == "fingerprint") %>%
  select(species, nbrood, avgbrood, minbrood, maxbrood, pmult, avgsire, minsire, maxsire,
         Reference)

write.csv(paternity_birds_fp, file = "paternity_birds_fp.csv", row.names = FALSE)



## Socially monogamous (microsat-only) data
paternity_birds_mono <- bird_MP2 %>%
  rename(nbrood = Nbr, avgbrood = MeanBS, minbrood = MinBS, maxbrood = MaxBS,
         pmult = MP, avgsire = MeanS, minsire = MinS, maxsire = MaxS) %>%
  filter(DNA == "microsatellite", SocMon == "mono") %>%
  select(species, nbrood, avgbrood, minbrood, maxbrood, pmult, avgsire, minsire, maxsire,
         Reference)

write.csv(paternity_birds_mono, file = "paternity_birds_mono.csv", row.names = FALSE)


## NOT socially monogamous (microsat-only) data
paternity_birds_xmono <- bird_MP2 %>%
  rename(nbrood = Nbr, avgbrood = MeanBS, minbrood = MinBS, maxbrood = MaxBS,
         pmult = MP, avgsire = MeanS, minsire = MinS, maxsire = MaxS) %>%
  filter(DNA == "microsatellite", SocMon == "xmono") %>%
  select(species, nbrood, avgbrood, minbrood, maxbrood, pmult, avgsire, minsire, maxsire,
         Reference)

write.csv(paternity_birds_xmono, file = "paternity_birds_xmono.csv", row.names = FALSE)

