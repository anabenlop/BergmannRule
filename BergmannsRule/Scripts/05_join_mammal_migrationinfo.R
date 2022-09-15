##############################################################
# Authors: 
# Erin Henry, Ana Benitez-Lopez (@anabenlop)
# Email: erinhenry55@gmail.com, abenitez81@gmail.com, ana.benitez@ugr.es
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# https://www.anabenitezlopez.com/

##############################################################
# Description of script and instructions
##############################################################


# This script loads the mammals correlations dataset, migratory information from 
# Soriano-Redondo et al. (2020), https://www.nature.com/articles/s41467-020-19256-0; 
# Gnanadesikan et al. (2017), https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.3120
# Teitelbaum et al. (2015), https://onlinelibrary.wiley.com/doi/10.1111/ele.12435
# Bisson et al. (2009), https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007504
# and joins all datasets to run meta-regressions for resident vs migratory mammals (script 06_phylo_metaRegressions)
# Species synonyms to match both datasets are retrieved using the package taxize

# clean environment
rm(list = ls())

# Packages and working directory -----------------------------------------------
library(dplyr)
library(tidyverse)
library(taxize)

# Load mammals data and migration info --------------------------------------------
#Load data
mammals_ph <- read.csv("Data/mamdata_ph.csv", stringsAsFactors = F)

# migration data from Soriano-Redondo et al. (2020)
mig_status <- read.csv("Data/Soriano-Redondo_mig_behaviour.csv", stringsAsFactors = FALSE) 
mig_status <- mig_status %>%
              mutate(animal  = str_replace(animal , "_", " "))

# Filter mammals
mig_mam <- mig_status %>%
          filter(Class == "Mammal") # 540 mammal species

# match species from mammals dataset with species in migration dataset 1 (Soriano-Redondo et al. 2020)
mig_m <- left_join(mammals_ph, mig_mam[,c("animal","Mig_bi")], by = c("speciesname" = "animal"))
length(unique(mig_m[is.na(mig_m$Mig_bi),]$speciesname)) #450 species without mig status assigned

# migration data from Gnanadesikan et al 2018 (migratory dataset 2)
mig_status2 <- read.csv("Data/Gnanadesikan_2018_mammals_mig.csv", stringsAsFactors = FALSE) 
mig_status2$speciesname <- paste0(mig_status2$Genus," ", mig_status2$Species)

# match species from mammals dataset with species in migration dataset 2
mig_m2 <- left_join(mig_m, mig_status2[,c("speciesname","Movement")], by = "speciesname")
mig_m2$Mig_bi2 <- ifelse(mig_m2$Movement == "M",1, 
                            ifelse(mig_m2$Movement == "N", 0, NA)) 

mig_m2$Mig_status <- ifelse(is.na(mig_m2$Mig_bi), mig_m2$Mig_bi2, mig_m2$Mig_bi) 

length(unique(mig_m2[is.na(mig_m2$Mig_status),]$speciesname)) #393 species without mig status assigned

# migration data from Bisson et al. (2009)
mig_status3 <- read.csv("Data/Bisson_2009_Vespertilionidae_mig.csv", stringsAsFactors = FALSE) 

# match species from mammals dataset with species in migration dataset 2
mig_m3 <- left_join(mig_m2, mig_status3[,c("Species","Migration")], by = c("speciesname" = "Species"))
mig_m3$Mig_status <- ifelse(is.na(mig_m3$Mig_status), mig_m3$Migration, mig_m3$Mig_status) 

length(unique(mig_m3[is.na(mig_m3$Mig_status),]$speciesname)) #373 species without mig status assigned

# use taxize to find synonyms of species with missing migratory status
# Make batches so that the connection is not reset 

missing <- unique(mig_m3[is.na(mig_m3$Mig_status),]$speciesname)
# 
missing1 <- missing[1:50]
missing2 <- missing[51:100]
missing3 <- missing[101:150]
missing4 <- missing[151:200]
missing5 <- missing[201:250]
missing6 <- missing[251:300]
missing7 <- missing[301:373]

# 
miss_l <- list(missing1, missing2, missing3, missing4, missing5, missing6, missing7)
syn <- list()
 
for (i in 1:length(miss_l)) {
  syn[[i]] <- synonyms(miss_l[[i]], db = "itis")
  }

syn_df1 <- synonyms_df(syn[[1]])
syn_df2 <- synonyms_df(syn[[2]])
syn_df3 <- synonyms_df(syn[[3]])
syn_df4 <- synonyms_df(syn[[4]])
syn_df5 <- synonyms_df(syn[[5]])
syn_df6 <- synonyms_df(syn[[6]])
syn_df7 <- synonyms_df(syn[[7]])

# 
# #create extra columns in syn_df4 and syn_df5 to bind alltogether
# syn_df4$acc_name <- NA
# syn_df4$acc_author <- NA
# syn_df5$acc_name <- NA
# syn_df5$acc_author <- NA
# 
syn_df <- rbind(syn_df1, syn_df2, syn_df2, syn_df3, syn_df4, syn_df5, syn_df6, syn_df7)
                                                                                            # 
# # save species for which we found a syn in itis
write.csv(syn_df, "Data/syn_itis_mam2.csv", row.names = F)

# load here synonyms dataset to avoid searching again
syn_df <- read.csv("Data/syn_itis_mam2.csv", stringsAsFactors = F)

# join syn with migratory dataset 1 and keep those that match
temp <- left_join(syn_df[,c(".id", "syn_name", "acc_name")], mig_mam[,c("animal","Mig_bi")], by = c("syn_name" = "animal"))
temp2 <- left_join(temp, mig_mam[,c("animal","Mig_bi")], by = c("acc_name" = "animal"))
temp2$Mig_bi <- ifelse(is.na(temp2$Mig_bi.x), temp2$Mig_bi.y,temp2$Mig_bi.x)

# join syn + mig1 with migratory dataset 2 and keep those that match
temp3 <- left_join(temp2[,c(".id", "syn_name", "acc_name", "Mig_bi")], mig_status2[,c("speciesname","Movement")], by = c("syn_name" = "speciesname")) 
temp4 <- left_join(temp3, mig_status2[,c("speciesname","Movement")], by = c("acc_name" = "speciesname")) 
temp4$Mig_bi.x <- ifelse(temp4$Movement.x == "M", 1, 
                         ifelse(temp4$Movement.x == "N", 0, NA)) 
temp4$Mig_bi.y <- ifelse(temp4$Movement.y == "M", 1, 
                               ifelse(temp4$Movement.y == "N", 0, NA))

temp4$Mig_bi <- ifelse(is.na(temp4$Mig_bi), temp4$Mig_bi.x,
                       ifelse(is.na(temp4$Mig_bi.x), temp4$Mig_bi.y,temp4$Mig_bi)) # 13 species with mig info

# join syn + mig1 + mig 2 with migratory dataset 3 (Ancillotto and Santini) and keep those that match
temp5 <- left_join(temp4[,c(".id", "syn_name", "acc_name", "Mig_bi")], mig_status3[,c("speciesname","Mig")], by = c("syn_name" = "speciesname")) 
temp5 <- left_join(temp5, mig_status3[,c("speciesname","Mig")], by = c("acc_name" = "speciesname")) 

temp5$Mig_bi <- ifelse(is.na(temp5$Mig_bi), temp5$Mig.x,temp5$Mig_bi)
temp5$Mig_bi <- ifelse(is.na(temp5$Mig_bi), temp5$Mig.y,temp5$Mig_bi)

# keep only those that were retrieved
temp5 <- temp5[!is.na(temp5$Mig_bi),]
syn_mig <- temp5[,c(".id", "syn_name", "acc_name","Mig_bi")]

# syn_mig <- syn_mig %>%  # this is done if more than one syn is retrieved
#   group_by(.id) %>% 
#   summarize(mig = first(mig))

# save species with mig info based on itis and Soriano-Redondo et al. (2020)
write.csv(syn_mig, "Data/syn_mig_mam.csv", row.names = F)

# join original dataset with mig info, with missing species with mig status fixed after finding synonyms in itis
temp6 <- left_join(mig_m3,syn_mig, by = c("speciesname" = ".id"))

temp6$Mig_status <- ifelse(!is.na(temp6$Mig_status), temp6$Mig_status, temp6$Mig_bi.y) 

# remove duplicates
temp6 <- temp6[!duplicated(temp6),]

# there are still some species for which we do not have migratory info
miss_mam <- unique(temp6[is.na(temp6$Mig_status),]$speciesname) # 263 species need to be fixed manually
missing_df <- temp6[is.na(temp6$Mig_status),]
missing_df <- missing_df[,c("speciesname", "class", "order", "family", "freq" , "env.var", "corr.coeff",
                            "z.cor.yi", "z.cor.vi",
                            "sd.tavg", "sd.tmax",  "sd.npp", "sd.npp.sd", "rng.tavg", "rng.tmax", "rng.npp", "rng.npp.sd",
                            "SPID","Mig_status")]

write.csv(missing_df, "Data/missing_mig_mam.csv", row.names = F)

# fix species without mig status manually (list of species in miss_mam)
temp6[temp6$speciesname  == "Ursus americanus", "Mig_status"] <- 1 # Noyce, K. V., & Garshelis, D. L. (2011). Seasonal migrations of black bears (Ursus americanus): causes and consequences. Behavioral ecology and sociobiology, 65(4), 823-835.

# we will assume that all those missing are resident. we create a new var 
temp6$Mig_status2 <- ifelse(is.na(temp6$Mig_status), 0, temp6$Mig_status)

# temp6[temp6$speciesname  == "Abrothrix longipilis", "Mig_status"] <- "resident" #

mammals_ph_mig <- temp6[,c("speciesname", "class", "order", "family", "freq" , "env.var", "corr.coeff",
                         "z.cor.yi", "z.cor.vi",
                         "sd.tavg", "sd.tmax",  "sd.npp", "sd.npp.sd", "rng.tavg", "rng.tmax", "rng.npp", "rng.npp.sd",
                         "SPID","Mig_status", "Mig_status2")]

mammals_ph_mig$Mig_status <- as.character(mammals_ph_mig$Mig_status)
mammals_ph_mig$Mig_status2 <- as.character(mammals_ph_mig$Mig_status2)

mammals_ph_mig$Mig_status<- ifelse(mammals_ph_mig$Mig_status == "1", "Migratory", 
                                     ifelse(mammals_ph_mig$Mig_status == "0", "Resident", NA))
mammals_ph_mig$Mig_status2 <- ifelse(mammals_ph_mig$Mig_status2 == "1", "Migratory", "Resident")

write.csv(mammals_ph_mig, "Data/mammals_ph_mig.csv", row.names = F)

# End of script --- 