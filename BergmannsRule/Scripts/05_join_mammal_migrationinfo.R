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
mig_m3$Mig_bi3 <- ifelse(mig_m3$Migration == 0, 0, 1)
                         
mig_m3$Mig_status <- ifelse(is.na(mig_m3$Mig_status), mig_m3$Mig_bi3, mig_m3$Mig_status) 

length(unique(mig_m3[is.na(mig_m3$Mig_status),]$speciesname)) #373 species without mig status assigned

# use taxize to find synonyms of species with missing migratory status
# Make batches so that the connection is not reset 

# missing <- unique(mig_m3[is.na(mig_m3$Mig_status),]$speciesname)
# # 
# missing1 <- missing[1:50]
# missing2 <- missing[51:100]
# missing3 <- missing[101:150]
# missing4 <- missing[151:200]
# missing5 <- missing[201:250]
# missing6 <- missing[251:300]
# missing7 <- missing[301:373]
# 
# # 
# miss_l <- list(missing1, missing2, missing3, missing4, missing5, missing6, missing7)
# syn <- list()
#  
# for (i in 1:length(miss_l)) {
#   syn[[i]] <- synonyms(miss_l[[i]], db = "itis")
#   }
# 
# syn_df1 <- synonyms_df(syn[[1]])
# syn_df2 <- synonyms_df(syn[[2]])
# syn_df3 <- synonyms_df(syn[[3]])
# syn_df4 <- synonyms_df(syn[[4]])
# syn_df5 <- synonyms_df(syn[[5]])
# syn_df6 <- synonyms_df(syn[[6]])
# syn_df7 <- synonyms_df(syn[[7]])
# 
# # 
# #create extra columns in syn_df4 and syn_df5 to bind all together
# syn_df2$acc_name <- NA
# syn_df2$acc_author <- NA
# syn_df7$acc_name <- NA
# syn_df7$acc_author <- NA
# # 
# syn_df <- rbind(syn_df1, syn_df2, syn_df2, syn_df3, syn_df4, syn_df5, syn_df6, syn_df7)
                                                                                            # 
# # save species for which we found a syn in itis
# write.csv(syn_df, "Data/syn_itis_mam.csv", row.names = F)

# load here synonyms dataset to avoid searching again
syn_df <- read.csv("Data/syn_itis_mam.csv", stringsAsFactors = F)

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

# join syn + mig1 + mig 2 with migratory dataset 3 (Bisson et al. 2009) and keep those that match
temp5 <- left_join(temp4[,c(".id", "syn_name", "acc_name", "Mig_bi")], mig_status3[,c("Species","Migration")], by = c("syn_name" = "Species")) 
temp5 <- left_join(temp5, mig_status3[,c("Species","Migration")], by = c("acc_name" = "Species")) 
temp5$Migration.x <- ifelse(temp5$Migration.x == 0, 0, 1)
temp5$Migration.y <- ifelse(temp5$Migration.y == 0, 0, 1)

temp5$Mig_bi <- ifelse(is.na(temp5$Mig_bi), temp5$Migration.x,
                       ifelse(is.na(temp5$Mig_bi.x), temp5$Migration.y, temp5$Mig_bi))
                              

# keep only those that were retrieved
temp5 <- temp5[!is.na(temp5$Mig_bi),]
syn_mig <- temp5[,c(".id", "syn_name", "acc_name","Mig_bi")]
syn_mig$Mig_bi <- ifelse(syn_mig$Mig_bi == 0, 0, 1)

# Collapse synonyms
syn_mig <- syn_mig %>%  # this is done if more than one syn is retrieved
  group_by(.id) %>%
  summarize(mig = first(Mig_bi))

# save species with mig info based on itis and Soriano-Redondo et al. (2020), Gnanadesikan et al. (2017), and Bisson et al. (2009)
write.csv(syn_mig, "Data/syn_mig_mam.csv", row.names = F)

# join original dataset with mig info, with missing species with mig status fixed after finding synonyms in itis
temp6 <- left_join(mig_m3,syn_mig, by = c("speciesname" = ".id"))

temp6$Mig_status <- ifelse(!is.na(temp6$Mig_status), temp6$Mig_status, temp6$mig) 

# remove duplicates
temp6 <- temp6[!duplicated(temp6),]

# there are still some species for which we do not have migratory info
miss_mam <- unique(temp6[is.na(temp6$Mig_status),]$speciesname) # 367 species need to be fixed manually
missing_df <- temp6[is.na(temp6$Mig_status),]
missing_df <- missing_df[,c("speciesname", "class", "order", "family", "freq" , "env.var", "corr.coeff",
                            "z.cor.yi", "z.cor.vi",
                            "sd.tavg", "sd.tmax",  "sd.npp", "sd.npp.sd", "rng.tavg", "rng.tmax", "rng.npp", "rng.npp.sd",
                            "SPID","Mig_status")]

write.csv(missing_df, "Data/missing_mig_mam.csv", row.names = F)

# fix species without mig status manually (list of species in miss_mam)
temp6[temp6$speciesname  == "Ursus americanus", "Mig_status"] <- 1 # Noyce, K. V., & Garshelis, D. L. (2011). Seasonal migrations of black bears (Ursus americanus): causes and consequences. Behavioral ecology and sociobiology, 65(4), 823-835.
temp6[temp6$speciesname  == "Pteropus alecto", "Mig_status"] <- 1 # Popa-Lisseanu, A. G., & Voigt, C. C. (2009). Bats on the move. Journal of Mammalogy, 90(6), 1283-1289.
temp6[temp6$speciesname  == "Parastrellus hesperus", "Mig_status"] <- 0 # These bats generally are regarded as nonmigratory, but seasonal status is poorly known in the northwestern portion of the range in Washington and Oregon (Verts, B. J., and L. N. Carraway. 1998. Land mammals of Oregon. University of California Press, Berkeley)
temp6[temp6$speciesname  == "Chaerephon pumilus", "Mig_status"] <- 0 # https://animaldiversity.org/accounts/Chaerephon_pumilus/
temp6[temp6$speciesname  == "Eumops perotis", "Mig_status"] <- 0 # https://animaldiversity.org/accounts/Eumops_perotis/
temp6[temp6$speciesname  == "Molossus rufus", "Mig_status"] <- 0 # https://sta.uwi.edu/fst/lifesciences/sites/default/files/lifesciences/documents/ogatt/Molossus_rufus%20-%20Black%20Mastiff%20Bat.pdf
temp6[temp6$speciesname  == "Mormoops megalophylla", "Mig_status"] <- 1 # https://animaldiversity.org/accounts/Mormoops_megalophylla/
temp6[temp6$speciesname  == "Noctilio albiventris", "Mig_status"] <- 0 # https://animaldiversity.org/accounts/Noctilio_albiventris/
temp6[temp6$speciesname  == "Noctilio leporinus", "Mig_status"] <- 0 # Villavicencio, C. H. N. (2019). Unusual Record For Greater Bulldog Bat Noctilio leporinus Linnaeus, 1758 (Chiroptera: Noctilionidae) in the Southern Andes of Ecuador. ACI Avances en Ciencias e Ingenierías, 11(3), 34-39.
temp6[temp6$speciesname  == "Desmodus rotundus", "Mig_status"] <- 0 # Fraser, K. C., McKinnon, E. A., & Diamond, A. W. (2010). Migration, diet, or molt? Interpreting stable‐hydrogen isotope values in Neotropical bats. Biotropica, 42(4), 512-517.
temp6[temp6$speciesname  == "Glossophaga longirostris", "Mig_status"] <- 0 # Sosa, M., & Soriano, P. J. (1993). Solapamiento de dieta entre Leptonycteris curasoae y Glossophaga longirostris (Mammalia: Chiroptera). Revista de Biología Tropical, 529-532.
temp6[temp6$speciesname  == "Glossophaga soricina", "Mig_status"] <- 0 # Fleming, T. H., Nuñez, R. A., & Sternberg, L. D. S. L. (1993). Seasonal changes in the diets of migrant and non-migrant nectarivorous bats as revealed by carbon stable isotope analysis. Oecologia, 94(1), 72-75.
temp6[temp6$speciesname  == "Hylonycteris underwoodi", "Mig_status"] <- 0 # https://animaldiversity.org/accounts/Hylonycteris_underwoodi/
temp6[temp6$speciesname  == "Lophostoma silvicolum", "Mig_status"] <- 0 # Moussy, C., Hosken, D. J., Mathews, F., Smith, G. C., Aegerter, J. N., & Bearhop, S. (2013). Migration and dispersal patterns of bats and their influence on genetic structure. Mammal Review, 43(3), 183-195.
temp6[temp6$speciesname  == "Macrotus waterhousii", "Mig_status"] <- 0 # https://animalia.bio/waterhouses-leaf-nosed-bat
temp6[temp6$speciesname  == "Micronycteris megalotis", "Mig_status"] <- 0 # Holanda, G. M., Oliveira, E. H. C. D., & Ribeiro, N. A. B. (2012). Dispersão geográfica da família Phyllostomidae (Chiroptera) baseado nas sequências do cifocromo b. Revista Pan-Amazônica de Saúde, 3(3), 21-31.
temp6[temp6$speciesname  == "Phyllostomus hastatus", "Mig_status"] <- 0 # Williams, T. C., Williams, J. M., & Griffin, D. R. (1966). The homing ability of the neotropical bat Phyllostomus hastatus, with evidence for visual orientation. Animal Behaviour, 14(4), 468-473.
temp6[temp6$speciesname  == "Platyrrhinus helleri", "Mig_status"] <- 0 # McGuire, L. P., & Boyle, W. A. (2013). Altitudinal migration in bats: evidence, patterns, and drivers. Biological Reviews, 88(4), 767-786.
temp6[temp6$speciesname  == "Platyrrhinus lineatus", "Mig_status"] <- 1 # Esbérard, C. E., Godoy, M. S., Renovato, L., & Carvalho, W. D. (2017). Novel long-distance movements by Neotropical bats (Mammalia: Chiroptera: Phyllostomidae) evidenced by recaptures in southeastern Brazil. Studies on neotropical fauna and environment, 52(1), 75-80.
temp6[temp6$speciesname  == "Pygoderma bilabiatum", "Mig_status"] <- 0 # Esbérard, C. E., Lima, I. P. D., Nobre, P. H., Althoff, S. L., Jordão-Nogueira, T., Dias, D., ... & Sobrinho, A. S. (2011). Evidence of vertical migration in the Ipanema bat Pygoderma bilabiatum (Chiroptera: Phyllostomidae: Stenodermatinae). Zoologia (curitiba), 28, 717-724.
temp6[temp6$speciesname  == "Rhinophylla pumilio", "Mig_status"] <- 0 # https://animaldiversity.org/accounts/Rhinophylla_pumilio/
temp6[temp6$speciesname  == "Trachops cirrhosus", "Mig_status"] <- 0 # https://animaldiversity.org/accounts/Trachops_cirrhosus/
temp6[temp6$speciesname  == "Uroderma bilobatum", "Mig_status"] <- 0 # Meyer, C. F. (2007). Effects of rainforest fragmentation on neotropical bats-land-bridge islands as a model system (Doctoral dissertation, Universität Ulm).
temp6[temp6$speciesname  == "Eonycteris spelaea", "Mig_status"] <- 0 # Waldien, D.L., Adleson, S. & Wilson, Z. 2020. Eonycteris spelaea. The IUCN Red List ofThreatened Species 2020: e.T7787A22128326. https://dx.doi.org/10.2305/IUCN.UK.2020-3.RLTS.T7787A22128326.en
temp6[temp6$speciesname  == "Nyctimene albiventer", "Mig_status"] <- 0 # Pattiselanno, F. Chiropteran Mammals: Representing the Lowland Rain Forest at the Northern Region of Papua, Indonesia.
temp6[temp6$speciesname  == "Ptenochirus jagori", "Mig_status"] <- 0 # Moussy, C., Hosken, D. J., Mathews, F., Smith, G. C., Aegerter, J. N., & Bearhop, S. (2013). Migration and dispersal patterns of bats and their influence on genetic structure. Mammal Review, 43(3), 183-195.
temp6[temp6$speciesname  == "Pteropus alecto", "Mig_status"] <- 1 # Popa-Lisseanu, A. G., & Voigt, C. C. (2009). Bats on the move. Journal of Mammalogy, 90(6), 1283-1289.
temp6[temp6$speciesname  == "Rousettus amplexicaudatus", "Mig_status"] <- 1 # https://animaldiversity.org/accounts/Rousettus_amplexicaudatus/
temp6[temp6$speciesname  == "Syconycteris australis", "Mig_status"] <- 1 # https://animaldiversity.org/accounts/Syconycteris_australis/
temp6[temp6$speciesname  == "Aeorestes cinereus", "Mig_status"] <- 1 # https://www.depts.ttu.edu/nsrl/mammals-of-texas-online-edition/Accounts_Chiroptera/Aeorestes_cinereus.php#:~:text=These%20bats%20are%20migratory%20and,states%2C%20generally%20in%20montane%20areas

# we will assume that all those missing are resident. we create a new var 
temp6$Mig_status2 <- ifelse(is.na(temp6$Mig_status), 0, temp6$Mig_status)

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