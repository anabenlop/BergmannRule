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

# Load data ---------------------------------------------------------------------
# read in correlation results
results <- readRDS('Results/BergmannsRule_results_correlations_20211224.rds')

# load environmental variation per species
dat_env <- read.csv("Data/Bergmann_envvariation.csv", header = T, stringsAsFactors = F)


# 1. Amphibians ----------------------------------------------------------------

# subset results for amphibians
amphibians <- subset(results, class == 'amphibian')

amphibians <- amphibians[amphibians$env.var == 'prec' |
                           amphibians$env.var == 'npp'|
                            amphibians$env.var =='npp.sd',]

# join env var data
amphibians <- left_join(amphibians, dat_env, by = c("speciesname" = "species"))

write.csv(amphibians,"Data/amphibians.csv", row.names = F)

# 2. Reptiles ------------------------------------------------------------------

# subset results for reptiles
reptiles <- subset(results, class == 'reptile')

reptiles <- reptiles[reptiles$env.var == 'npp' | reptiles$env.var == 'npp.sd',]

# join env var data
reptiles <- left_join(reptiles, dat_env, by = c("speciesname" = "species"))

write.csv(reptiles,"Data/reptiles.csv", row.names = F)

# 3. Mammals -------------------------------------------------------------------

# subset results for mammals
mammals <- subset(results, class == 'mammal')

# read elton traits dataset and remove marine mammmals
elton_mam <- read.csv("Data/EltonTraits_Mammals_taxid.csv", header = T, stringsAsFactors = F)

mammals <- left_join(mammals, elton_mam[,c("Scientific", "ForStrat.Value")], by = c("speciesname" = "Scientific"))
mammals <- mammals[!mammals$ForStrat.Value %in% c('M'),]
  
# fix errors in taxonomic classification 
mammals[mammals$speciesname == "Notiosorex crawfordi", "order"] <- "Eulipotyphla"
mammals[mammals$speciesname == "Notiosorex crawfordi", "family"] <- "Soricidae"

# fix old order name, Insectivora is now Eulipotyphla
mammals$order <- ifelse(mammals$order == "Insectivora" | 
                          mammals$order == "Erinaceomorpha" | 
                          mammals$order == "Soricomorpha" , "Eulipotyphla", mammals$order)

# keep env.var of interest
mammals <- mammals[mammals$env.var == 'tavg' | 
                     mammals$env.var =='tmax' | 
                     mammals$env.var == 'npp'| 
                     mammals$env.var =='npp.sd',]

# join env var data
mammals <- left_join(mammals, dat_env, by = c("speciesname" = "species"))

write.csv(mammals,"Data/mammals.csv", row.names = F) # 578


# 4. Birds ---------------------------------------------------------------------

# subset results for birds
birds <- subset(results, class == 'bird')

# read elton traits dataset and remove marine mammmals
elton_bird <- read.csv("Data/BirdFuncDat.csv", header = T, stringsAsFactors = F)

birds <- left_join(birds, elton_bird[,c("Scientific", "PelagicSpecialist")], by = c("speciesname" = "Scientific"))
birds <-birds[c(birds$PelagicSpecialist == 0 | is.na(birds$PelagicSpecialist)),]
birds <- birds[birds$family != "Pelecanidae", ]
birds <- birds[birds$env.var == 'tavg' | 
                 birds$env.var =='tmax' | 
                 birds$env.var == 'npp'| 
                 birds$env.var =='npp.sd',]

# join env var data
birds <- left_join(birds, dat_env, by = c("speciesname" = "species"))

write.csv(birds,"Data/birds.csv", row.names = F)

#5. Herps ----------------------------------------------------------------------
  
# subset results for herps
herps <- subset(results, class == 'amphibian' | class == 'reptile')
herps <- herps[herps$env.var == 'tavg',] # we only test the heat balance hyp for thermoconformers vs thermoregulators
    
# join env var data
herps <- left_join(herps, dat_env, by = c("speciesname" = "species"))

write.csv(herps,"Data/herps.csv", row.names = F)


# End of script ----