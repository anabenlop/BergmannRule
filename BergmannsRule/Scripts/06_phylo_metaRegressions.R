###### Meta Regressions v5 ######

# Created on 29 September 2020
## Modified 14 November 2021 ##

# Packages and working directory -----------------------------------------------
library(dplyr)
library(metafor)

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

# load bird data with migratory info
birds <- read.csv("Data/birds_ph_mig.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/bird_phylo_cor.Rdata") #bird_phylo_cor

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = bird_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# Tavg model
bird.tavg <- rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = birds,
                subset = env.var=="tavg",
                mods = ~ migratory-1,
                random = RE, R = phylocor)
summary(bird.tavg)

# save results
saveRDS(bird.tavg,"Results/BergmannsRule_results_MR_mig_tavg.rds")
rm(bird.tavg)

# Tavg model - no nomadic species
bird.tavg_n <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = birds,
                    subset = env.var=="tavg" & migratory != "nomadic",
                    mods = ~ migratory-1,
                    random = RE, R = phylocor)
summary(bird.tavg_n)

# save results
saveRDS(bird.tavg_n,"Results/BergmannsRule_results_MR_mig_tavg_n.rds")
rm(bird.tavg_n)

# Tmax model
bird.tmax <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = birds,
                    subset = env.var=="tmax",
                    mods = ~ migratory - 1,
                    random = RE,  R = phylocor)
summary(bird.tmax)

# save results
saveRDS(bird.tmax,'Results/BergmannsRule_results_MR_mig_tmax.rds')
rm(bird.tmax)

# Tmax model - no nomadic species
bird.tmax_n <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = birds,
                    subset = env.var=="tmax" & migratory != "nomadic",
                    mods = ~ migratory - 1,
                    random = RE,  R = phylocor)
summary(bird.tmax_n)

# save results
saveRDS(bird.tmax_n,'Results/BergmannsRule_results_MR_mig_tmax_n.rds')
rm(bird.tmax_n)

# NPP model
bird.npp <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = birds,
                    subset = env.var=='npp',
                    mods = ~ migratory - 1,
                    random = RE,  R = phylocor)
summary(bird.npp)

saveRDS(bird.npp,'Results/BergmannsRule_results_MR_mig_npp.rds')
rm(bird.npp)

# NPP model - no nomadic
bird.npp_n <- rma.mv(yi = z.cor.yi,
                   V = z.cor.vi,
                   data = birds,
                   subset = env.var=='npp' & migratory != "nomadic",
                   mods = ~ migratory - 1,
                   random = RE,  R = phylocor)
summary(bird.npp_n)

saveRDS(bird.npp_n,'Results/BergmannsRule_results_MR_mig_npp_n.rds')
rm(bird.npp_n)

# NPP sd model
bird.npp.sd <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = birds, 
                    subset = env.var=='npp.sd',
                    mods = ~ migratory - 1,
                    random = RE,  R = phylocor)
summary(bird.npp.sd)

saveRDS(bird.npp.sd,'Results/BergmannsRule_results_MR_mig_nppsd.rds')
rm(bird.npp.sd)

# NPP sd model -  no nomadic
bird.npp.sd_n <- rma.mv(yi = z.cor.yi,
                      V = z.cor.vi,
                      data = birds, 
                      subset = env.var=='npp.sd' & migratory != "nomadic",
                      mods = ~ migratory - 1,
                      random = RE,  R = phylocor)
summary(bird.npp.sd_n)

saveRDS(bird.npp.sd_n,'Results/BergmannsRule_results_MR_mig_nppsd_n.rds')
rm(bird.npp.sd_n)

# load fitted models
bird.tavg <- readRDS('Results/BergmannsRule_results_MR_mig_tavg.rds')
bird.tavg_n <- readRDS('Results/BergmannsRule_results_MR_mig_tavg_n.rds')
bird.tmax <- readRDS('Results/BergmannsRule_results_MR_mig_tmax.rds')
bird.tmax_n <- readRDS('Results/BergmannsRule_results_MR_mig_tmax_n.rds')
bird.npp <- readRDS('Results/BergmannsRule_results_MR_mig_npp.rds')
bird.npp_n <- readRDS('Results/BergmannsRule_results_MR_mig_npp_n.rds')
bird.npp.sd <- readRDS('Results/BergmannsRule_results_MR_mig_nppsd.rds')
bird.npp.sd_n <- readRDS('Results/BergmannsRule_results_MR_mig_nppsd_n.rds')

# save results in list
bi.mr <- list(bird.tavg_n,bird.tmax_n,bird.npp_n,bird.npp.sd_n)
names(bi.mr) <- c("tavg","tmax","npp","npp.sd")

saveRDS(bi.mr,'Results/BergmannsRule_results_MR_mig.rds')

# 3. Test mammals with body mass  ------------------------------------------------
# Freckleton et al. 2003 --> Bergmann’s rule does depend on body mass: larger species follow the rule more commonly than smaller ones

# load data
mammals_ph <- read.csv("Data/mamdata_ph.csv", stringsAsFactors = F)
# mammals_ph$Species_ph <- gsub(" ", "_", trimws(mammals_ph$speciesname))

# read elton traits dataset and save body mass
elton_mam <- read.csv("Data/EltonTraits_Mammals_taxid.csv", header = T, stringsAsFactors = F)

mammals_ph <- left_join(mammals_ph, elton_mam[,c("Scientific", "BM")], by = c("speciesname" = "Scientific"))




# loading phylogenetic matrixes 
load("Data/mam_phylo_cor.Rdata") #mam_phylo_cor

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = mam_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# Tavg model
mam.tavg.bm <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = birds,
                    subset = env.var=="tavg",
                    mods = ~ migratory-1,
                    random = RE, R = phylocor)
summary(bird.tavg)

# save results
saveRDS(bird.tavg,"Results/BergmannsRule_results_MR_mig_tavg.rds")
rm(bird.tavg)


# End of script

