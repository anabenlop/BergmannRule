##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 30th of December 2021

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script loads the bird correlations dataset, migratory information from Eyres et al. 2017 
# (https://onlinelibrary.wiley.com/doi/pdf/10.1111/jav.01308; http://dataportal-senckenberg.de/database/metacat/bikf.10058.1/bikf)
# and joins both datasets to run metaregressions for resident vs migratory birds (script 06_phylo_metaRegressions)
# for the paper:


# Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Uncovering the environmental drivers 
# of intraspecific body size variation in terrestrial vertebrates. 

# clean environment
rm(list = ls())

# Packages and working directory -----------------------------------------------
library(dplyr)
library(metafor)

# 2a. Add  migration info for birds --------------------------------------------
#Load data
birds_ph <- read.csv("Data/birddata_ph.csv", stringsAsFactors = F)

# migration data from Eyres et al. 2017
mig_status<-read.csv("Data/SpeciesList3_1_migbehav_v2_0.csv", stringsAsFactors = FALSE) 
mig_status$speciesname <- paste0(mig_status$IOC3_1_Genus," ", mig_status$IOC3_1_Species)

# match species from birds dataset with species in migration dataset 
mig_b <- left_join(birds_ph, mig_status[,c("speciesname","Migratory_status", "Migratory_status_3")], by = "speciesname")
length(unique(mig_b[is.na(mig_b$Migratory_status),]$speciesname)) #111 species without mig status assigned
length(unique(mig_b[is.na(mig_b$Migratory_status_3),]$speciesname))#111 species without mig status 3 assigned

missing <- unique(mig_b[is.na(mig_b$Migratory_status_3),]$speciesname)

# save species without migratory status 
write.csv(missing, "Data/missing_mig.csv", row.names = F)

# fix species without mig status
birddata_temp[birddata_temp$Binomial == "Acanthis flammea", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") #based on christinae, which includes latouchii
birddata_temp[birddata_temp$Binomial == "Aethopyga latouchii", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") #based on christinae, which includes latouchii
birddata_temp[birddata_temp$Binomial == "Alcedo azurea", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # recorded as Ceyx azureus 
birddata_temp[birddata_temp$Binomial == "Alcippe cinereiceps", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # All members in genus Alcippe are resident
birddata_temp[birddata_temp$Binomial == "Psittacara leucophthalmus", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # recorded as Aratinga leucophtalma
birddata_temp[birddata_temp$Binomial == "Eupsittula astec", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Includes Aratinga nana/astec
birddata_temp[birddata_temp$Binomial == "Zanda funerea", c("Migratory_status", "Migratory_status_3")] <- c("resident","partial resident") # Recorded as Calyptorhynchus funereus
birddata_temp[birddata_temp$Binomial == "Synoicus ypsilophorus", c("Migratory_status", "Migratory_status_3")] <- c("resident","partial resident") # Recorded as Coturnix ypsilophora
birddata_temp[birddata_temp$Binomial == "Dryobates scalaris", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Picoides scalaris
birddata_temp[birddata_temp$Binomial == "Dyaphorophyia castanea", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Platysteira castanea
birddata_temp[birddata_temp$Binomial == "Dyaphorophyia chalybea", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Platysteira chalybea
birddata_temp[birddata_temp$Binomial == "Garrulax elliotii", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Trochalopteron elliotii
birddata_temp[birddata_temp$Binomial == "Haemorhous mexicanus", c("Migratory_status", "Migratory_status_3")] <- c("directional migratory", "partial directional migrant") # Recorded as Carpodacus mexicanus
birddata_temp[birddata_temp$Binomial == "Hemicircus sordidus", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as	Hemicircus concretus/sordidus
birddata_temp[birddata_temp$Binomial == "Nesoptilotis leucotis", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Lichenostomus leucotis
birddata_temp[birddata_temp$Binomial == "Macronous kelleyi", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Macronus kelleyi
birddata_temp[birddata_temp$Binomial == "Macronous ptilosus", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Macronus ptilosus
birddata_temp[birddata_temp$Binomial == "Meiglyptes grammithorax", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Meiglyptes tristis
birddata_temp[birddata_temp$Binomial == "Pardaliparus venustulus", c("Migratory_status", "Migratory_status_3")] <- c("dispersive migratory", "partial dispersive migrant") # Recorded as Periparus venustulus
birddata_temp[birddata_temp$Binomial == "Zoothera interpres", c("Migratory_status", "Migratory_status_3")] <- c("resident","partial resident") # Recorded as Geockila interpres
birddata_temp[birddata_temp$Binomial == "Turdus libonyanus", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Turdus libonyana
birddata_temp[birddata_temp$Binomial == "Trogon ambiguus", c("Migratory_status", "Migratory_status_3")] <- c("directional migratory", "partial directional migrant") # Included in Trogon elegans
birddata_temp[birddata_temp$Binomial == "Passer italiae", c("Migratory_status", "Migratory_status_3")] <- c("directional migratory", "partial directional migrant") # Included in Passer hispaniolensis
birddata_temp[birddata_temp$Binomial == "Sephanoides sephaniodes", c("Migratory_status", "Migratory_status_3")] <- c("directional migratory", "partial directional migrant") # Recorded as Sephanoides sephanoides
birddata_temp[birddata_temp$Binomial == "Spinus psaltria", c("Migratory_status", "Migratory_status_3")] <- c("directional migratory", "partial directional migrant") # Recorded as Carduelis psaltria
birddata_temp[birddata_temp$Binomial == "Trichastoma tickelli", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Pellorneum tickelli
birddata_temp[birddata_temp$Binomial == "Trichoglossus meyeri", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Included in Trichoglossus flavoviridis


