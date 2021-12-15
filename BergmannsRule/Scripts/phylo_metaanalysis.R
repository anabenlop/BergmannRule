##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 12th of December 2021

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script runs the phylogenetic meta-analysis for mammals for the paper: 


# Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Uncovering the environmental drivers 
# of intraspecific body size variation in terrestrial vertebrates. 


##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
# library(ggplot2)
# library(ggpubr)
library(tictoc)
# library(png)
# library(ggimage)
# library(grid)
# library(rsvg)
# library(grImport2)

#clean memory
rm(list=ls())

##############################################################
# Importing datasets                                      ####
##############################################################


# 3. Mammals ---------------------------------------------------------------------

#Load data
mammals_ph <- read.csv("Data/mammals_ph.csv", stringsAsFactors = F)
mammals_ph$Species_ph <- gsub(" ", "_", trimws(mammals_ph$speciesname))

# loading phylogenetic matrixes 
load("Data/mam_phylo_cor.Rdata") #mam_phylo_cor

# vector of environmental variables
env.vars <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

# define phylo vcov matrix and random effects
phylocor<-list(Species_ph  = mam_phylo_cor)
RE = list( ~1|speciesname, ~1|Species_ph)

# for loop running a meta-analysis for each environmental variable
tic("Run phylo meta-analysis in a loop")
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('mam.',env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(mammals_ph, env.var == env.vars[i]),
                random = RE, R = phylocor))
}

toc()

# save results
ma.ma <- list(mam.tavg, mam.tmin, mam.tmax, mam.prec, 
              mam.pet, mam.npp, mam.npp.sd)
names(ma.ma) <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

saveRDS(ma.ma,
        'Results/BergmannsRule_results_MA_mammals_phylo_nonphylo.rds')

# check results
# ma.ma_ph_nonph <- readRDS("Results/BergmannsRule_results_MA_mammals_phylo_nonphylo.rds")


# 4. Birds ---------------------------------------------------------------------

#Load data
birds_ph <- read.csv("Data/birds_ph.csv", stringsAsFactors = F)
# birds_ph$Species_ph <- gsub(" ", "_", trimws(birds_ph$speciesname))

# loading phylogenetic matrixes 
load("Data/bird_phylo_cor.Rdata") #bird_phylo_cor

# vector of environmental variables
env.vars <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

# define phylo vcov matrix and random effects
phylocor<-list(Species_ph  = bird_phylo_cor)
RE = list( ~1|speciesname, ~1|Species_ph)

# for loop running a meta-analysis for each environmental variable
tic("Run phylo meta-analysis in a loop")
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('bird.',env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(birds_ph, env.var == env.vars[i]),
                random = RE, R = phylocor))
}

toc()

# save results
bi.ma <- list(bird.tavg, bird.tmin, bird.tmax, bird.prec, 
              bird.pet, bird.npp, bird.npp.sd)
names(bi.ma) <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

saveRDS(bi.ma,
        'Results/BergmannsRule_results_MA_birds_phylo_nonphylo.rds')

# check results
# ma.ma_ph_nonph <- readRDS("Results/BergmannsRule_results_MA_birds_phylo_nonphylo.rds")
