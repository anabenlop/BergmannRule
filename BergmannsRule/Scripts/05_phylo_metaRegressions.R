###### Meta Regressions v5 ######

# Created on 29 September 2020
## Modified 14 November 2021 ##

# Packages and working directory -----------------------------------------------
#library(plyr)
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

# 2a. Add  migration info for birds --------------------------------------------

# migration data from IUCN
migration <- read.csv("Data/bird_range_types.csv")

# species id data
ids <- read.csv("Data/birds_id.csv")

#Load data
birds_ph <- read.csv("Data/birddata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/bird_phylo_cor.Rdata") #bird_phylo_cor

# match ids in migration dataset to species names in ids
mig_b <- left_join(ids, migration, by = "tax_id")

# match mig_b with correlations dataset
birds_join <- left_join(birds_ph, mig_b, by = c("speciesname" = "Species"))



# change underscores to spaces for species names in ID data frame
# ids$Species <- gsub("_", " ",ids$Species)

# new column for species ID
# birds$id <- NA

# for loop: if species is in ID list, then add respective ID;  
# for(i in 1:nrow(birds)) {
#   if(birds[i,"speciesname"] %in% ids$Species) { 
#     birds[i,"id"] <- ids$tax_id[which(ids$Species==birds$speciesname[i])]
#     
#   } else { # else id = NA
#     birds[i,"id"] <- NA
#   }
# }

# match ids for NAs with ids from Species_ph
for(i in 1:nrow(birds)){
  if(is.na(birds[i,"id"])) {
    birds[i,"id"] <- birds[i,"tax_id"]
  }
}


# check NA ids
(birds.to.id <- subset(birds,is.na(id))[,c("speciesname","Species_RL",
                                           "Species_ph","tax_id")])
# save to find ids manually
#write.csv(birds.to.id,
#"C:/Users/Erin/Desktop/Bergmann's Rule/birds-to-id-v5.csv")

# read in data frame with ids for NAs
birds.to.id <- read.csv("Data/birds-to-id-v5.csv")

# check for duplicates with synonyms ***
x <- subset(birds, speciesname %in% birds.to.id$synonym)
nrow(x)
birds.to.id[which(birds.to.id$synonym=="Aphelocoma californica"),"speciesname"]
rm(x)

# match missing ids
for(i in 1:nrow(birds)) {
  
  # if species is in the list of missing ids, add id
  if(birds[i,"speciesname"] %in% birds.to.id$speciesname) { 
    birds[i,"id"] <- 
      birds.to.id[which(birds.to.id$speciesname==birds$speciesname[i]),"tax_id"]
  } 
}

# check (should be 0)
subset(birds,is.na(id))[,c("speciesname","Species_RL","Species_ph","tax_id")]

# Not 0, so add additional information from a second file
birds.to.id2 <- read.csv("Data/birds-to-id-v5-addition.csv")
for(i in 1:nrow(birds)) {
  
  # if species is in the list of missing ids, add id
  if(birds[i,"speciesname"] %in% birds.to.id2$speciesname) { 
    birds[i,"id"] <- 
      birds.to.id2[which(birds.to.id2$speciesname==birds$speciesname[i]),"tax_id"]
  } 
}

# check again (should be 0)
subset(birds, is.na(id))[,c("speciesname","Species_RL","Species_ph","tax_id")]


# add migration data using ids
birds$resident <- NA
for(i in 1:nrow(birds)){
  if(birds[i,"id"]%in%migration$tax_id) {
    birds[i,"resident"] <- 
      migration[which(migration$tax_id==birds[i,"id"]),"resident"]
  }
}

# new column: change true/false to resident/migratory
birds$migration <- ifelse(birds$resident,"Resident","Migratory")

# check for NAs
subset(birds,is.na(resident))[,c("speciesname","Species_RL","Species_ph","tax_id","resident","migration")]

# add migration info manually (using IUCN)
birds[which(birds$speciesname=="Calandrella cheleensis"),"resident"] <- FALSE
birds[which(birds$speciesname=="Calandrella cheleensis"),"migration"] <- "Migratory"
birds[which(birds$speciesname=="Gallicolumba beccarii"),"resident"] <- TRUE
birds[which(birds$speciesname=="Gallicolumba beccarii"),"migration"] <- "Resident"
birds[which(birds$speciesname=="Oceanodroma leucorhoa"),"resident"] <- FALSE
birds[which(birds$speciesname=="Oceanodroma leucorhoa"),"migration"] <- "Migratory"
birds[which(birds$speciesname=="Picoides dorsalis"),"resident"] <- TRUE
birds[which(birds$speciesname=="Picoides dorsalis"),"migration"] <- "Resident"

# check again (should be 0)
nrow(subset(birds,is.na(resident))[,c("speciesname","Species_RL","Species_ph","tax_id","resident","migration")])
nrow(subset(birds,is.na(resident)))/7 # divide by 7 to account for only one var


# Not 0, but only 36 species are NA. OK for now

# save data
saveRDS(birds,"Results/BergmannsRule_results_correlationsBirdsMR_20211115.rds")


# 2b. Migration Meta-Regression ------------------------------------------------

# correlation data
birds <- readRDS("Results/BergmannsRule_results_correlationsBirdsMR_20211115.rds")

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


