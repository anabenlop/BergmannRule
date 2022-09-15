##############################################################
# Authors: 
# Erin Henry, Ana Benitez-Lopez (@anabenlop)
# Email: erinhenry55@gmail.com, abenitez81@gmail.com, ana.benitez@ugr.es
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# https://www.anabenitezlopez.com/

##############################################################
# Description of script and instructions
##############################################################

# This script runs the phylogenetic meta-analyses for each taxonomic group for testing several 
# existing hypotheses explaining within-species body size variation in terrestrial vertebrates, 
# including the heat balance, seasonality, resource availability and water conservation hypotheses 
# for ectotherms, and the heat conservation, heat dissipation, starvation resistance, and 
# resource availability hypotheses for endotherms.

##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
library(tictoc)

#clean environment
rm(list=ls())

##############################################################
# Importing datasets                                      ####
##############################################################
# 1. Amphibians ---------------------------------------------------------------------

#Load data
amphibians_ph <- read.csv("Data/amphdata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/Phylogeny/amph_phylo_cor.Rdata") #amph_phylo_cor

# vector of environmental variables
env.vars <- c('prec','npp','npp.sd')

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = amph_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# for loop running a meta-analysis for each environmental variable
tic("Run phylo meta-analysis in a loop")
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('amph.',env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(amphibians_ph, env.var == env.vars[i]),
                random = RE, R = phylocor))
}

toc()

# save results
amph.ma <- list(amph.prec, amph.npp, amph.npp.sd)
names(amph.ma) <- c('prec','npp','npp.sd')

saveRDS(amph.ma,
        'Results/BergmannsRule_results_MA_amphibians_phylo_nonphylo.rds')

# check results
# amph.ma_ph_nonph <- readRDS("Results/BergmannsRule_results_MA_amphibians_phylo_nonphylo.rds")

# 2. Reptiles ---------------------------------------------------------------------

#Load data
reptiles_ph <- read.csv("Data/reptdata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/Phylogeny/rept_phylo_cor.Rdata") #rept_phylo_cor

# vector of environmental variables
env.vars <- c('npp','npp.sd')

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = rept_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# for loop running a meta-analysis for each environmental variable
tic("Run phylo meta-analysis in a loop")
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('rept.',env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(reptiles_ph, env.var == env.vars[i]),
                random = RE, R = phylocor))
}

toc()

# save results
rept.ma <- list(rept.npp, rept.npp.sd)
names(rept.ma) <- c('npp','npp.sd')

saveRDS(rept.ma,
        'Results/BergmannsRule_results_MA_reptiles_phylo_nonphylo.rds')

# check results
# rept.ma_ph_nonph <- readRDS("Results/BergmannsRule_results_MA_reptiles_phylo_nonphylo.rds")

# 3. Mammals ---------------------------------------------------------------------

#Load data
mammals_ph <- read.csv("Data/mamdata_ph.csv", stringsAsFactors = F) # 567 species

# loading phylogenetic matrixes 
load("Data/Phylogeny/mam_phylo_cor.Rdata") #mam_phylo_cor

# vector of environmental variables
env.vars <- c('tavg','tmax','npp','npp.sd')

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = mam_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

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
ma.ma <- list(mam.tavg, mam.tmax, mam.npp, mam.npp.sd)
names(ma.ma) <- c('tavg','tmax','npp','npp.sd')

saveRDS(ma.ma,
        'Results/BergmannsRule_results_MA_mammals_phylo_nonphylo.rds')

# check results
# ma.ma_ph_nonph <- readRDS("Results/BergmannsRule_results_MA_mammals_phylo_nonphylo.rds")


# 4. Birds ---------------------------------------------------------------------

#Load data
birds_ph <- read.csv("Data/birddata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/Phylogeny/bird_phylo_cor.Rdata") #bird_phylo_cor

# vector of environmental variables
env.vars <- c('tavg','tmax','npp','npp.sd')

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = bird_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# for loop running a meta-analysis for each environmental variable
tic("Run phylo meta-analysis in a loop") #this may take long time to run
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
bi.ma <- list(bird.tavg, bird.tmax, bird.npp, bird.npp.sd)
names(bi.ma) <- c('tavg','tmax','npp','npp.sd')

saveRDS(bi.ma,
        'Results/BergmannsRule_results_MA_birds_phylo_nonphylo.rds')

# check results
# bi.ma_ph_nonph <- readRDS("Results/BergmannsRule_results_MA_birds_phylo_nonphylo.rds")

# End of script --- 
