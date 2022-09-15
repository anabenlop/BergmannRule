##############################################################
# Authors: 
# Erin Henry, Ana Benitez-Lopez (@anabenlop)
# Email: erinhenry55@gmail.com, abenitez81@gmail.com, ana.benitez@ugr.es
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# https://www.anabenitezlopez.com/

##############################################################
# Description of script and instructions
##############################################################

# This script loads the bird correlations dataset, migratory information from Eyres et al. 2017 
# (https://onlinelibrary.wiley.com/doi/pdf/10.1111/jav.01308; http://dataportal-senckenberg.de/database/metacat/bikf.10058.1/bikf)
# and joins both datasets to run metaregressions for resident vs migratory birds (script 06_phylo_metaRegressions)


# clean environment
rm(list = ls())

# Packages and working directory -----------------------------------------------
library(dplyr)
library(taxize)

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
# missing <- data.frame(speciesname = missing)

# use taxize to find synonyms of species with missing migratory status
# syn <- synonyms(missing, db = "itis")
# syn_df <- synonyms_df(syn) # 62 retrieved
# 
# # save species for which we found a syn in itis
# write.csv(syn_df, "Data/syn_itis.csv", row.names = F)

# load here synonyms dataset to avoid searching again
syn_df <- read.csv("Data/syn_itis.csv", stringsAsFactors = F)

# join syn with migratory dataset and keep those that match
temp <- left_join(syn_df[,c(".id", "syn_name", "acc_name")], mig_status[,c("speciesname","Migratory_status", "Migratory_status_3")], by = c("syn_name" = "speciesname"))
temp2 <- left_join(temp, mig_status[,c("speciesname","Migratory_status", "Migratory_status_3")], by = c("acc_name" = "speciesname"))
temp2$Migratory_status_3 <- ifelse(is.na(temp2$Migratory_status_3.x), temp2$Migratory_status_3.y,temp2$Migratory_status_3.x)
temp2$Migratory_status <- ifelse(is.na(temp2$Migratory_status.x), temp2$Migratory_status.y,temp2$Migratory_status.x)
temp2 <- temp2[!is.na(temp2$Migratory_status),]

# remove species with unknown mig status (extinct)
temp2 <- temp2[temp2$Migratory_status != "unknown",]
temp2$mig <- ifelse(temp2$Migratory_status == "resident", "resident", "migratory")
syn_mig <- temp2[,c(".id", "syn_name", "acc_name","mig")]

syn_mig <- syn_mig %>% 
            group_by(.id) %>% 
            summarize(mig = first(mig))

# save species with mig info based on itis and Eyres et al. 2017
write.csv(syn_mig, "Data/syn_mig.csv", row.names = F)

# repeat for another dataset (nbn)
# syn2 <- synonyms(missing, db = "nbn") # it takes quite some time
# syn_df2 <- synonyms_df(syn2) # 62 retrieved

# # save species for which we found a syn in itis
# write.csv(syn_df2, "Data/syn_nbn.csv", row.names = F)

# load here synonyms dataset to avoid searching again
syn_df2 <- read.csv("Data/syn_nbn.csv", stringsAsFactors = F)

# join syn with migratory dataset and keep those that match
temp3 <- left_join(syn_df2[,c(".id", "nameString")], mig_status[,c("speciesname","Migratory_status", "Migratory_status_3")], by = c("nameString" = "speciesname"))
temp3 <- temp3[!is.na(temp3$Migratory_status),]

# add mig column based on migratory status
temp3$mig <- ifelse(temp3$Migratory_status == "resident", "resident", "migratory")
syn_mig2 <- temp3[,c(".id", "nameString", "mig")]

syn_mig2 <- syn_mig2 %>% 
  group_by(.id) %>% 
  summarize(mig = first(mig))

# save species with mig info based on nbn and Eyres et al. 2017
write.csv(syn_mig2, "Data/syn_mig2.csv", row.names = F)

# join original dataset with mig info, with missing species with mig status fixed after finding synonyms in itis and nbn
temp4 <- left_join(mig_b,syn_mig, by = c("speciesname" = ".id"))
temp5 <- left_join(temp4,syn_mig2, by = c("speciesname" = ".id"))

temp5$Migratory_status <- ifelse(temp5$Migratory_status == "resident", "resident",
                                 ifelse(temp5$Migratory_status == "nomadic","nomadic", "migratory"))

temp5$migratory <- ifelse(!is.na(temp5$Migratory_status), temp5$Migratory_status, 
                         ifelse(!is.na(temp5$mig.x),temp5$mig.x, temp5$mig.y)) 

missing2 <- unique(temp5[is.na(temp5$migratory),]$speciesname) # 45 species need to be fixed manually

# fix species without mig status manually (list of species in missing2)
temp5[temp5$speciesname  == "Alcedo cyanopectus", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/inbkin2/cur/movement
temp5[temp5$speciesname  == "Aramides cajanea", "migratory"] <- "resident" #https://birdsoftheworld.org/bow/species/gycwor1/cur/introduction#mignat
temp5[temp5$speciesname  == "Alaudala cheleensis", "migratory"] <- "migratory" # Eyres et al. 2017, recorded as Calandrella cheleensis
temp5[temp5$speciesname  == "Chlorospingus flavopectus", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/cobtan1/cur/introduction#mig
temp5[temp5$speciesname  == "Kittacincla luzoniensis", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/whbsha1/cur/introduction#mig, recorded as Copsychus luzoniensis
temp5[temp5$speciesname  == "Melloria quoyi", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/blabut1/cur/introduction#mig, recorded as Cracticus quoyi
temp5[temp5$speciesname  == "Origma murina", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/rumwar1/cur/introduction#mig, recorded as Crateroscelis murina
temp5[temp5$speciesname  == "Origma robusta", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/momwar1/cur/introduction#mig, recorded as Crateroscelis robusta
temp5[temp5$speciesname  == "Cyanoloxia cyanoides", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/bubgro1/cur/distribution
temp5[temp5$speciesname  == "Cyanoderma ruficeps", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/rucbab1/cur/introduction#mig
temp5[temp5$speciesname  == "Schoeniclus elegans", "migratory"] <- "migratory" #  Eyres et al. 2017, recorded as Emberiza elegans
temp5[temp5$speciesname  == "Schoeniclus spodocephala", "migratory"] <- "migratory" #  Eyres et al. 2017, recorded as Emberiza spodocephala
temp5[temp5$speciesname  == "Fringillaria tahapisi", "migratory"] <- "migratory" #  Eyres et al. 2017, recorded as Emberiza tahapisi
temp5[temp5$speciesname  == "Alopecoenas beccarii", "migratory"] <- "resident" #  https://birdsoftheworld.org/bow/species/brgdov1/cur/introduction#mignat (lumped taxon)
temp5[temp5$speciesname  == "Setopagis parvulus", "migratory"] <- "migratory" #  Eyres et al. 2017, recorded as Caprimulgus parvulus
temp5[temp5$speciesname  == "Kurochkinegramma hypogrammicum", "migratory"] <- "resident" #  Eyres et al. 2017, recorded as Hypogramma hypogrammicum
temp5[temp5$speciesname  == "Bolemoreus frenatus", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/brihon1/cur/introduction#mig
temp5[temp5$speciesname  == "Mixornis gularis", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/sttbab1/cur/introduction#mig
temp5[temp5$speciesname  == "Melaenornis chocolatinus", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Dioptrornis chocolatinus
temp5[temp5$speciesname  == "Kieneria crissalis", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Melozone crissalis
temp5[temp5$speciesname  == "Kieneria fusca", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Melozone fusca
temp5[temp5$speciesname  == "Sporophila angolensis", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/cbsfin/cur/introduction#mig
temp5[temp5$speciesname  == "Otus flammeolus", "migratory"] <- "migratory" # Eyres et al. 2017, recorded as Megascops flammeolus
temp5[temp5$speciesname  == "Pardaliparus elegans", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Periparus elegans
temp5[temp5$speciesname  == "Machlolophus spilonotus", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Parus spilonotus
temp5[temp5$speciesname  == "Nannopterum auritus", "migratory"] <- "migratory" # Eyres et al. 2017, recorded as Phalacrocorax auritus
temp5[temp5$speciesname  == "Phyllergates cucullatus", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Phyllergates heterolaemus
temp5[temp5$speciesname  == "Seicercus borealis", "migratory"] <- "migratory" # Eyres et al. 2017, recorded as Phylloscopus borealis
temp5[temp5$speciesname  == "Abrornis inornata", "migratory"] <- "migratory" # Eyres et al. 2017, recorded as Phylloscopus inornatus
temp5[temp5$speciesname  == "Abrornis maculipennis", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Phylloscopus maculipennis
temp5[temp5$speciesname  == "Seicercus tenellipes", "migratory"] <- "migratory" # Eyres et al. 2017, recorded as Phylloscopus tenellipes
temp5[temp5$speciesname  == "Pipra erythrocephala", "migratory"] <- "resident" # https://birdsoftheworld.org/bow/species/gohman1/cur/distribution, recorded as Ceratopipra erythrocephala
temp5[temp5$speciesname  == "Pipraeidea bonariensis", "migratory"] <- "migratory" # Eyres et al. 2017, recorded as Thraupis bonariensis
temp5[temp5$speciesname  == "Poecilodryas albispecularis", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Heteromyias albispecularis
temp5[temp5$speciesname  == "Rhipidura spilodera", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Rhipidura verreauxi
temp5[temp5$speciesname  == "Myrmelastes leucostigma", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Schistocichla leucostigma
temp5[temp5$speciesname  == "Aethomyias arfakianus", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Sericornis arfakianus
temp5[temp5$speciesname  == "Neosericornis citreogularis", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Sericornis citreogularis
temp5[temp5$speciesname  == "Aethomyias spilodera", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Sericornis spilodera
temp5[temp5$speciesname  == "Curruca subcoerulea", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Sylvia subcaerulea
temp5[temp5$speciesname  == "Luscinia cyanura", "migratory"] <- "migratory" # Eyres et al. 2017, recorded as Tarsiger cyanurus
temp5[temp5$speciesname  == "Tangara episcopus", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Thraupis episcopus
temp5[temp5$speciesname  == "Tangara palmarum", "migratory"] <- "resident" # Eyres et al. 2017, recorded as Thraupis palmarum
temp5[temp5$speciesname  == "Tangara sayaca", "migratory"] <- "migratory" # Eyres et al. 2017, recorded as Thraupis sayaca
temp5[temp5$speciesname  == "Troglodytes musculus", "migratory"] <- "migratory" # Eyres et al. 2017, recorded as Troglodytes aedon

birds_ph_mig <- temp5[,c("speciesname", "class", "order", "family", "freq" , "env.var", "corr.coeff",
                         "z.cor.yi", "z.cor.vi",
                         "sd.tavg", "sd.tmax",  "sd.npp", "sd.npp.sd", "rng.tavg", "rng.tmax", "rng.npp", "rng.npp.sd",
                         "SPID","migratory")]

write.csv(birds_ph_mig, "Data/birds_ph_mig.csv", row.names = F)

# End of script --- 