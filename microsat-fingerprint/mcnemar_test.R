### McNemar test for DNA methods
library(tidyverse)
library(readxl)
paternity_birds <- read_excel("~/Documents/Research/Masamu/paternity/birds/paternity_birds.xlsx", 
                              col_types = c("skip", "text", "text", 
                                            "text", "text", "numeric", "text", 
                                            "text", "text", "text", "text", "text", 
                                            "numeric", "text", "text", "numeric", 
                                            "numeric", "numeric", "numeric", 
                                            "numeric", "numeric", "numeric", 
                                            "text", "text", "text", "text", "text", 
                                            "text", "numeric"))
# View(paternity_birds)

## Species represented more than once in the data
multDNA <- paternity_birds %>%
  group_by(species) %>%
  filter(n() > 1)

## Of those species represented more than once, microsat only
multMS <- paternity_birds %>%
  group_by(species) %>%
  filter(n() > 1 & DNA == "microsatellite")

## Of those species represented more than once, fingerprint only
multFP <- paternity_birds %>%
  group_by(species) %>%
  filter(n() > 1 & DNA == "fingerprint")

## Species found in both "FP" and "MS" datasets 
## (had at least one study of each DNA method)
bothDNA_species <- intersect(multMS$species, multFP$species)

bothDNA <- multDNA %>%
  group_by(species) %>%
  filter(species %in% bothDNA_species)

## Number of unique species with both DNA methods
unique(bothDNA$species)

bothDNA$MP_cat <- ifelse(bothDNA$MP==0, "zero", "positive")

## Table of species with count of studies by DNA method and MP (zero or positve)
DNAtabl <- bothDNA %>% group_by(species, DNA, MP_cat) %>% summarise(num_studies = n())

## Contingency table of species by MS/FP and zeroMP/positiveMP?
with(DNAtabl, table(DNA, MP_cat))

## Table of species DNA-MP interactions
DNApmult <- DNAtabl %>% unite(DNApmult, c("DNA", "MP_cat"))


## How many species had MP=zero in both FP and MS?
DNApmult %>% group_by(species) %>%
  filter(DNApmult=="fingerprint_zero" | DNApmult=="microsatellite_zero") %>% 
  filter(n() > 1) %>% summarise(n()) %>% count() #one, Lesser kestral

## How many species had MP=zero in FP and positiveMP in MS?
DNApmult %>% group_by(species) %>%
  filter(DNApmult=="fingerprint_zero" | DNApmult=="microsatellite_positive") %>% 
  filter(n() > 1) %>% summarise(n()) %>% count() #four species

## How many species had MP=positive in FP and MP=zero in MS?
DNApmult %>% group_by(species) %>%
  filter(DNApmult=="fingerprint_positive" | DNApmult=="microsatellite_zero") %>% 
  filter(n() > 1) %>% summarise(n()) %>% count() #one, Kentish plover

## How many species had MP=positive in both FP and MS?
DNApmult %>% group_by(species) %>%
  filter(DNApmult=="fingerprint_positive" | DNApmult=="microsatellite_positive") %>% 
  filter(n() > 1) %>% summarise(n()) %>% count() #eleven species

DNA_MP <- matrix(c(1,4,1,11), byrow = TRUE, nrow = 2, ncol = 2, 
                 dimnames = list("Fingerprint" = c("ZeroMP","PositiveMP"),
                                 "Microsatellite" = c("ZeroMP","PositiveMP")))

mcnemar.test(DNA_MP)


