##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estacion Biologica de Donana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 15th of December 2021

##############################################################
# Description of script and instructions
##############################################################

# This script is takes the correlations calculated in the script
# 02_BergmannsRule_correlationAnalysis.R and prepares the data for the phylogenetic
# meta-analysis for the paper:

# Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Uncovering the environmental drivers 
# of intraspecific body size variation in terrestrial vertebrates. 

##############################################################
# Packages needed
##############################################################

#load libraries
library(dplyr)

# clean environment
rm(list=ls())

# 1. Amphibians ----------------------------------------------------------------

# read in correlation results
results <- readRDS('Results/BergmannsRule_results_correlations_20211114.rds')

# subset results for amphibians
amphibians <- subset(results, class == 'amphibian')
# amphibians$Species_ph <- gsub(" ", "_", trimws(amphibians$speciesname))
amphibians <- amphibians[amphibians$env.var == 'prec' |
                           amphibians$env.var == 'npp'|
                            amphibians$env.var =='npp.sd',]

write.csv(amphibians,"Data/amphibians.csv", row.names = F)

# 2. Reptiles ------------------------------------------------------------------

# read in correlation results
results <- readRDS('Results/BergmannsRule_results_correlations_20211114.rds')

# subset results for reptiles
reptiles <- subset(results, class == 'reptile')
# reptiles$Species_ph <- gsub(" ", "_", trimws(reptiles$speciesname))
reptiles <- reptiles[reptiles$env.var == 'npp' | reptiles$env.var == 'npp.sd',]

write.csv(reptiles,"Data/reptiles.csv", row.names = F)

# 3. Mammals -------------------------------------------------------------------

# read in correlation results
results <- readRDS('Results/BergmannsRule_results_correlations_20211114.rds')

# subset results for mammals
mammals <- subset(results, class == 'mammal')

# read elton traits dataset and remove marine mammmals
elton_mam <- read.csv("Data/EltonTraits_Mammals_taxid.csv", header = T, stringsAsFactors = F)

mammals <- left_join(mammals, elton_mam[,c("Scientific", "ForStrat.Value")], by = c("speciesname" = "Scientific"))
mammals <-mammals[mammals$ForStrat.Value != "M",]
mammals <- mammals[mammals$env.var == 'tavg' | 
                     mammals$env.var =='tmax' | 
                     mammals$env.var == 'npp'| 
                     mammals$env.var =='npp.sd',]

mammals <-mammals[complete.cases(mammals$corr.coeff),]

write.csv(mammals,"Data/mammals.csv", row.names = F)


# 4. Birds ---------------------------------------------------------------------

# read in correlation results
results <- readRDS('Results/BergmannsRule_results_correlations_20211114.rds')

# subset results for birds
birds <- subset(results, class == 'bird')
#birds$Species_ph <- gsub(" ", "_", trimws(birds$speciesname))

# read elton traits dataset and remove marine mammmals
elton_bird <- read.csv("Data/BirdFuncDat.csv", header = T, stringsAsFactors = F)

birds <- left_join(birds, elton_bird[,c("Scientific", "PelagicSpecialist")], by = c("speciesname" = "Scientific"))
birds <-birds[c(birds$PelagicSpecialist == 0 | is.na(birds$PelagicSpecialist)),]
birds <- birds[birds$family != "Pelecanidae", ]
birds <- birds[birds$env.var == 'tavg' | 
                 birds$env.var =='tmax' | 
                 birds$env.var == 'npp'| 
                 birds$env.var =='npp.sd',]

birds <-birds[complete.cases(birds$corr.coeff),]

write.csv(birds,"Data/birds.csv", row.names = F)

# End of script ----