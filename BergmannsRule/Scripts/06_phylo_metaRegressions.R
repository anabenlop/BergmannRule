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
mammals <- read.csv("Data/mamdata_ph.csv", stringsAsFactors = F)
# mammals_ph$Species_ph <- gsub(" ", "_", trimws(mammals_ph$speciesname))

# read elton traits dataset and save body mass
elton_mam <- read.csv("Data/EltonTraits_Mammals_taxid.csv", header = T, stringsAsFactors = F)

mammals <- left_join(mammals, elton_mam[,c("Scientific", "BodyMass.Value")], by = c("speciesname" = "Scientific"))

# join missing species
mis_bm <- read.csv("Data/mammals_bm.csv", stringsAsFactors = F) 
mammals <- left_join(mammals, mis_bm, by = "speciesname")
mammals$BodyMass.Value <-ifelse(is.na(mammals$BodyMass.Value), mammals$BodyMass, mammals$BodyMass.Value)

mammals$BM <- log10(mammals$BodyMass.Value) 


# loading phylogenetic matrixes 
load("Data/mam_phylo_cor.Rdata") #mam_phylo_cor

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = mam_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# Tavg model
mam.tavg.bm <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = mammals,
                    subset = env.var=="tavg",
                    mods = ~ BM,
                    random = RE, R = phylocor)
summary(mam.tavg.bm)

# save results
saveRDS(mam.tavg.bm,"Results/BergmannsRule_results_MR_mammals_BM_tavg.rds")
rm(mam.tavg.bm)

# Tmaxmodel
mam.tmax.bm <- rma.mv(yi = z.cor.yi,
                      V = z.cor.vi,
                      data = mammals,
                      subset = env.var=="tmax",
                      mods = ~ BM,
                      random = RE, R = phylocor)
summary(mam.tmax.bm) #marginally significant --> size-tmax corr becomes more negative for larger species

# save results
saveRDS(mam.tmax.bm,"Results/BergmannsRule_results_MR_mammals_BM_tmax.rds")
rm(mam.tmax.bm)


# NPP model
mam.npp.bm <- rma.mv(yi = z.cor.yi,
                      V = z.cor.vi,
                      data = mammals,
                      subset = env.var=="npp",
                      mods = ~ BM,
                      random = RE, R = phylocor)
summary(mam.npp.bm) 

# save results
saveRDS(mam.npp.bm,"Results/BergmannsRule_results_MR_mammals_BM_npp.rds")
rm(mam.npp.bm)



# NPPsd model
mam.nppsd.bm <- rma.mv(yi = z.cor.yi,
                      V = z.cor.vi,
                      data = mammals,
                      subset = env.var=="npp.sd",
                      mods = ~ BM,
                      random = RE, R = phylocor)
summary(mam.nppsd.bm) 


# save results
saveRDS(mam.nppsd.bm,"Results/BergmannsRule_results_MR_mammals_BM_nppsd.rds")
rm(mam.nppsd.bm)


# 4. Test effect of environmental variation  ------------------------------------------------

#### a) Mammals ####
# load mammal dataset 
mamdata <- read.csv("Data/mamdata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/mam_phylo_cor.Rdata") #bird_phylo_cor

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = mam_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# Tmean model with env variation
mam.tavg.env <- rma.mv(yi = z.cor.yi,
                    V = z.cor.vi,
                    data = mamdata,
                    subset = env.var=="tavg",
                    mods = ~ sd.tavg,
                    random = RE, R = phylocor)
summary(mam.tavg.env)

# save results
saveRDS(mam.tavg.env,"Results/BergmannsRule_results_MR_mam_tavg_env.rds")
rm(mam.tavg.env)


# Tmax model with env variation
mam.tmax.env <- rma.mv(yi = z.cor.yi,
                       V = z.cor.vi,
                       data = mamdata,
                       subset = env.var=="tmax",
                       mods = ~ sd.tmax,
                       random = RE, R = phylocor)
summary(mam.tmax.env)

# save results
saveRDS(mam.tmax.env,"Results/BergmannsRule_results_MR_mam_tmax_env.rds")
rm(mam.tmax.env)

# npp model with env variation
mam.npp.env <- rma.mv(yi = z.cor.yi,
                       V = z.cor.vi,
                       data = mamdata,
                       subset = env.var=="npp",
                       mods = ~ log10(sd.npp),
                       random = RE, R = phylocor)
summary(mam.npp.env)

# save results
saveRDS(mam.npp.env,"Results/BergmannsRule_results_MR_mam_npp_env.rds")
rm(mam.npp.env)

# npp.sd model with env variation
mam.npp.sd.env <- rma.mv(yi = z.cor.yi,
                      V = z.cor.vi,
                      data = mamdata,
                      subset = env.var=="npp.sd",
                      mods = ~ log10(sd.npp.sd),
                      random = RE, R = phylocor)
summary(mam.npp.sd.env)

# save results
saveRDS(mam.npp.sd.env,"Results/BergmannsRule_results_MR_mam_nppsd_env.rds")
rm(mam.npp.sd.env)

#### b) Birds ####
# load birds dataset 
birddata <- read.csv("Data/birddata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/bird_phylo_cor.Rdata") #bird_phylo_cor

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = bird_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# Tmean model with env variation
bird.tavg.env <- rma.mv(yi = z.cor.yi,
                       V = z.cor.vi,
                       data = birddata,
                       subset = env.var=="tavg",
                       mods = ~ sd.tavg,
                       random = RE, R = phylocor)
summary(bird.tavg.env)

# save results
saveRDS(bird.tavg.env,"Results/BergmannsRule_results_MR_bird_tavg_env.rds")
rm(bird.tavg.env)


# Tmax model with env variation
bird.tmax.env <- rma.mv(yi = z.cor.yi,
                       V = z.cor.vi,
                       data = birddata,
                       subset = env.var=="tmax",
                       mods = ~ sd.tmax,
                       random = RE, R = phylocor)
summary(bird.tmax.env) # negative signal for birds, the greater the environmental variation en max temp, the most likely is a species will have a stronger negative size-tmax relationship. 

# save results
saveRDS(bird.tmax.env,"Results/BergmannsRule_results_MR_bird_tmax_env.rds")
rm(bird.tmax.env)

# npp model with env variation
bird.npp.env <- rma.mv(yi = z.cor.yi,
                      V = z.cor.vi,
                      data = birddata,
                      subset = env.var=="npp",
                      mods = ~ log10(sd.npp),
                      random = RE, R = phylocor)
summary(bird.npp.env) #ns

# save results
saveRDS(bird.npp.env,"Results/BergmannsRule_results_MR_bird_npp_env.rds")
rm(bird.npp.env)

# npp.sd model with env variation
bird.npp.sd.env <- rma.mv(yi = z.cor.yi,
                         V = z.cor.vi,
                         data = birddata,
                         subset = env.var=="npp.sd",
                         mods = ~ log10(sd.npp.sd),
                         random = RE, R = phylocor)
summary(bird.npp.sd.env)

# save results
saveRDS(bird.npp.sd.env,"Results/BergmannsRule_results_MR_bird_nppsd_env.rds")
rm(bird.npp.sd.env)

#### c) Reptiles ####
# load reptiles dataset 
reptdata <- read.csv("Data/reptdata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/rept_phylo_cor.Rdata") #rept_phylo_cor

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = rept_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# npp model with env variation
rept.npp.env <- rma.mv(yi = z.cor.yi,
                       V = z.cor.vi,
                       data = reptdata,
                       subset = env.var=="npp",
                       mods = ~ log10(sd.npp),
                       random = RE, R = phylocor)
summary(rept.npp.env) # tendency to positive effect (large CI) with variation in npp across the gradient of body size records

# save results
saveRDS(rept.npp.env,"Results/BergmannsRule_results_MR_rept_npp_env.rds")
rm(rept.npp.env)

# npp.sd model with env variation
rept.npp.sd.env <- rma.mv(yi = z.cor.yi,
                          V = z.cor.vi,
                          data = reptdata,
                          subset = env.var=="npp.sd",
                          mods = ~ log10(sd.npp.sd),
                          random = RE, R = phylocor)
summary(rept.npp.sd.env) # no clear support for the seasonality hypothesis

# save results
saveRDS(rept.npp.sd.env,"Results/BergmannsRule_results_MR_rept_nppsd_env.rds")
rm(rept.npp.sd.env)

#### d) Amphibians ####
# load amphibians dataset 
amphdata <- read.csv("Data/amphdata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/amph_phylo_cor.Rdata") #amph_phylo_cor

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = amph_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# npp model with env variation
amph.npp.env <- rma.mv(yi = z.cor.yi,
                       V = z.cor.vi,
                       data = amphdata,
                       subset = env.var=="npp",
                       mods = ~ log10(sd.npp),
                       random = RE, R = phylocor)
summary(amph.npp.env) # negative effects, non significant

# save results
saveRDS(amph.npp.env,"Results/BergmannsRule_results_MR_amph_npp_env.rds")
rm(amph.npp.env)

# npp.sd model with env variation
amph.npp.sd.env <- rma.mv(yi = z.cor.yi,
                          V = z.cor.vi,
                          data = amphdata,
                          subset = env.var=="npp.sd",
                          mods = ~ log10(sd.npp.sd),
                          random = RE, R = phylocor)
summary(amph.npp.sd.env) # no clear support for the seasonality hypothesis

# prec model with env variation
amph.prec.env <- rma.mv(yi = z.cor.yi,
                          V = z.cor.vi,
                          data = amphdata,
                          subset = env.var=="prec",
                          mods = ~ sd.prec,
                          random = RE, R = phylocor)
summary(amph.prec.env) # tendency towards smaller size as species are exposed to more variation in precipitation

# save results
saveRDS(amph.prec.env,"Results/BergmannsRule_results_MR_amph_prec_env.rds")
rm(amph.prec.env)

# End of script ###

# End of script

