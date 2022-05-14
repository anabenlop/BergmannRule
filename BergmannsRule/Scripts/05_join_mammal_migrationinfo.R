##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# www.anabenitezlopez.com
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, Spain)
# Department of Zoology, University of Granada (UGR, Spain)
# Email: abenitez81@gmail.com; ana.benitez@ugr.es

# Script first created on the 13th of May 2022

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script loads the mammals correlations dataset, migratory information from Soriano-Redondo et al. (2020)
# https://www.nature.com/articles/s41467-020-19256-0
# and joins both datasets to run metaregressions for resident vs migratory mammals (script 06_phylo_metaRegressions)
# for the paper:

# Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Uncovering the environmental drivers 
# of intraspecific body size variation in terrestrial vertebrates. 

# clean environment
rm(list = ls())

# Packages and working directory -----------------------------------------------
library(dplyr)
library(tidyverse)
# library(taxize)

# Load mammals data and migration info --------------------------------------------
#Load data
mammals_ph <- read.csv("Data/mamdata_ph.csv", stringsAsFactors = F)

# migration data from Soriano-Redondo et al. (2020)
mig_status <- read.csv("Data/Soriano-Redondo_mig_mehaviour.csv", stringsAsFactors = FALSE) 
mig_status <- mig_status %>%
              mutate(animal  = str_replace(animal , "_", " "))

# Filter mammals
mig_mam <- mig_status %>%
          filter(Class == "Mammal") # 540 mammal species

# match species from mammals dataset with species in migration dataset 
mig_m <- left_join(mammals_ph, mig_mam[,c("animal","Mig_bi")], by = c("speciesname" = "animal"))
length(unique(mig_m[is.na(mig_m$Mig_bi),]$speciesname)) #450 species without mig status assigned

# migration data from Gnanadesikan et al 2018
mig_status2 <- read.csv("Data/Gnanadesikan_2018_mammals_mig.csv", stringsAsFactors = FALSE) 

###CONTINUE HERE####

# match species from mammals dataset with species in migration dataset 2
mig_m2 <- left_join(mig_m, mig_status2[,c("speciesname","Mig")], by = "speciesname")
mig_m2$Mig_status <- ifelse(is.na(mig_m2$Mig_bi), mig_m2$Mig, mig_m2$Mig_bi) 

length(unique(mig_m2[is.na(mig_m2$Mig_status),]$speciesname)) #312 species without mig status assigned

# migration data from Luca Santini and Leonardo Ancillotto (need a proper reference)
mig_status2 <- read.csv("Data/MIG_BEHAVIOUR_mam.csv", stringsAsFactors = FALSE) 

# match species from mammals dataset with species in migration dataset 2
mig_m2 <- left_join(mig_m, mig_status2[,c("speciesname","Mig")], by = "speciesname")
mig_m2$Mig_status <- ifelse(is.na(mig_m2$Mig_bi), mig_m2$Mig, mig_m2$Mig_bi) 

length(unique(mig_m2[is.na(mig_m2$Mig_status),]$speciesname)) #312 species without mig status assigned

missing <- unique(mig_m2[is.na(mig_m2$Mig_status),]$speciesname)
# missing <- data.frame(speciesname = missing)

# use taxize to find synonyms of species with missing migratory status
syn <- synonyms(missing, db = "itis")
syn_df <- synonyms_df(syn) # 62 retrieved
# 
# # save species for which we found a syn in itis
write.csv(syn_df, "Data/syn_itis_mam.csv", row.names = F)

# load here synonyms dataset to avoid searching again
syn_df <- read.csv("Data/syn_itis_mam.csv", stringsAsFactors = F)

# join syn with migratory dataset and keep those that match
temp <- left_join(syn_df[,c(".id", "syn_name", "acc_name")], mig_mam[,c("animal","Mig_bi")], by = c("syn_name" = "animal"))
temp2 <- left_join(temp, mig_mam[,c("animal","Mig_bi")], by = c("acc_name" = "animal"))
temp2$Mig_bi <- ifelse(is.na(temp2$Mig_bi.x), temp2$Mig_bi.y,temp2$Mig_bi.x)
temp2 <- temp2[!is.na(temp2$Mig_bi),]

# keep only those that were retrieved
syn_mig <- temp2[,c(".id", "syn_name", "acc_name","Mig_bi")]

# syn_mig <- syn_mig %>%  # this is done if more than one syn is retrieved
#   group_by(.id) %>% 
#   summarize(mig = first(mig))

# save species with mig info based on itis and Soriano-Redondo et al. (2020)
write.csv(syn_mig, "Data/syn_mig_mam.csv", row.names = F)

# repeat for another dataset (nbn) -- skip, too messy
# syn2 <- synonyms(missing, db = "nbn") # it takes quite some time
# syn_df2 <- synonyms_df(syn2) # 62 retrieved
# 
# # # save species for which we found a syn in itis
# # write.csv(syn_df2, "Data/syn_nbn.csv", row.names = F)
# 
# # load here synonyms dataset to avoid searching again
# syn_df2 <- read.csv("Data/syn_nbn.csv", stringsAsFactors = F)
# 
# # join syn with migratory dataset and keep those that match
# temp3 <- left_join(syn_df2[,c(".id", "nameString")], mig_status[,c("speciesname","Migratory_status", "Migratory_status_3")], by = c("nameString" = "speciesname"))
# temp3 <- temp3[!is.na(temp3$Migratory_status),]
# 
# # add mig column based on migratory status
# temp3$mig <- ifelse(temp3$Migratory_status == "resident", "resident", "migratory")
# syn_mig2 <- temp3[,c(".id", "nameString", "mig")]
# 
# syn_mig2 <- syn_mig2 %>% 
#   group_by(.id) %>% 
#   summarize(mig = first(mig))
# 
# # save species with mig info based on nbn and Eyres et al. 2017
# write.csv(syn_mig2, "Data/syn_mig2.csv", row.names = F)

# join original dataset with mig info, with missing species with mig status fixed after finding synonyms in itis and nbn
temp4 <- left_join(mig_m2,syn_mig, by = c("speciesname" = ".id"))
# temp5 <- left_join(temp4,syn_mig2, by = c("speciesname" = ".id"))

temp4$Mig_status <- ifelse(!is.na(temp4$Mig_status), temp4$Mig_status, temp4$Mig_bi.y) 

missing2 <- unique(temp4[is.na(temp4$Mig_status),]$speciesname) # 304 species need to be fixed manually
missing_df <- temp4[is.na(temp4$Mig_status),]
write.csv(missing_df, "Data/missing_mig_mam.csv", row.names = F)

# fix species without mig status manually (list of species in missing2)
temp5[temp5$speciesname  == "Alcedo cyanopectus", "migratory"] <- "resident" # https://mammalsoftheworld.org/bow/species/inbkin2/cur/movement

mammals_ph_mig <- temp5[,c("speciesname", "class", "order", "family", "freq" , "env.var", "corr.coeff",
                         "z.cor.yi", "z.cor.vi",
                         "sd.tavg", "sd.tmax",  "sd.npp", "sd.npp.sd", "rng.tavg", "rng.tmax", "rng.npp", "rng.npp.sd",
                         "SPID","migratory")]

write.csv(mammals_ph_mig, "Data/mammals_ph_mig.csv", row.names = F)

# End of script --- 