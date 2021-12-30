###### Meta Regressions v5 ######

# Created on 29 September 2020
## Modified 14 November 2021 ##

# Packages and working directory -----------------------------------------------
library(dplyr)
#library(taxize)
library(metafor)
#library(beepr)

# setwd('D:/BergmannsRule_upload')

# 1. Heat Balance Hypothesis ---------------------------------------------------

#Load data
ecto <- read.csv("Data/herpdata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/herp_phylo_cor.Rdata") # need to create a ectotherm phylogenetic correlationmatrix

# new column: thermoregulator (reptiles and anura) or thermoconformer (caudata)
ecto$therm <- ifelse(ecto$order=="Caudata",
                     "Thermoconf", # if TRUE
                     "Thermoreg") # if FALSE

# check
table(ecto$order, ecto$therm)

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = herp_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# Meta-regression comparing correlation coefficients (SVL - Tavg) between 
# thermoregulators and thermoconformers.
mod.hb <- rma.mv(yi = z.cor.yi,
                 V = z.cor.vi,
                 data = ecto, 
                 mods = ~ therm - 1,
                 random = RE, R = phylocor)

summary(mod.hb)

# save results
saveRDS(mod.hb,"Results/BergmannsRule_results_MR_heatBalance.rds")


# 2. Migration Meta-Regression ------------------------------------------------

# correlation data
birds <- readRDS("Results/BergmannsRule_results_correlationsBirdsMR_20211115.rds")

# loading phylogenetic matrixes 
load("Data/bird_phylo_cor.Rdata") #bird_phylo_cor


# Tavg model
bird.tavg <- rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(birds,env.var=="tavg"),
                mods = ~ migration - 1,
                random = list(~1|order/family/speciesname))
saveRDS(bird.tavg,"Results/BergmannsRule_results_MR_mig_tavg.rds")
rm(bird.tavg)

# Tmin model
bird.tmin <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = subset(birds, env.var=='tmin'),
                    mods = ~ migration - 1,
                    random = list(~1|order/family/speciesname))
saveRDS(bird.tmin,'Results/BergmannsRule_results_MR_mig_tmin.rds')
rm(bird.tmin)

# Tmax model
bird.tmax <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = subset(birds, env.var=='tmax'),
                    mods = ~ migration - 1,
                    random = list(~1|order/family/speciesname))
saveRDS(bird.tmax,'Results/BergmannsRule_results_MR_mig_tmax.rds')
rm(bird.tmax)

# Precipitation model
bird.prec <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = subset(birds, env.var=='prec'),
                    mods = ~ migration - 1,
                    random = list(~1|order/family/speciesname))
saveRDS(bird.prec,'Results/BergmannsRule_results_MR_mig_prec.rds')
rm(bird.prec)

# PET model
bird.pet <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = subset(birds, env.var=='pet'),
                    mods = ~ migration - 1,
                    random = list(~1|order/family/speciesname))
saveRDS(bird.pet,'Results/BergmannsRule_results_MR_mig_pet.rds')
rm(bird.pet)

# NPP model
bird.npp <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = subset(birds, env.var=='npp'),
                    mods = ~ migration - 1,
                    random = list(~1|order/family/speciesname))
saveRDS(bird.npp,'Results/BergmannsRule_results_MR_mig_npp.rds')
rm(bird.npp)

# NPP sd model
bird.npp.sd <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = subset(birds, env.var=='npp.sd'),
                    mods = ~ migration - 1,
                    random = list(~1|order/family/speciesname))
saveRDS(bird.npp.sd,'Results/BergmannsRule_results_MR_mig_nppsd.rds')
rm(bird.npp.sd)



bird.tavg <- readRDS('Results/BergmannsRule_results_MR_mig_tavg.rds')
bird.tmin <- readRDS('Results/BergmannsRule_results_MR_mig_tmin.rds')
bird.tmax <- readRDS('Results/BergmannsRule_results_MR_mig_tmax.rds')
bird.prec <- readRDS('Results/BergmannsRule_results_MR_mig_prec.rds')
bird.pet <- readRDS('Results/BergmannsRule_results_MR_mig_pet.rds')
bird.npp <- readRDS('Results/BergmannsRule_results_MR_mig_npp.rds')
bird.npp.sd <- readRDS('Results/BergmannsRule_results_MR_mig_nppsd.rds')

# save results in list
bi.mr <- list(bird.tavg,bird.tmin,bird.tmax,bird.prec,bird.pet,bird.npp,bird.npp.sd)
names(bi.mr) <- c("tavg","tmin","tmax","prec","pet","npp","npp.sd")

saveRDS(bi.mr,'Results/BergmannsRule_results_MR_mig.rds')

# remove objects
rm(list=ls(pattern="bird."))
rm(env.vars,i)


