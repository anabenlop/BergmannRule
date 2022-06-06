##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estacion Biologica de Donana (EBD-CSIC, ESP)
# Department of Zoology, University of Granada (UGR, Spain)
# Email: ana.benitez@ugr.es; abenitez81@gmail.com

# Script first created on the 2nd of June 2022

##############################################################
# Description of script and instructions
##############################################################

# This script is used to extract the family-level phylogeny for mammals for the paper:

# Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Uncovering the environmental drivers 
# of intraspecific body size variation in terrestrial vertebrates. 

##############################################################
# Packages needed
##############################################################

#load libraries
library(dplyr)
library(ape)
library(treebase)
library(rotl)
library(diagram)
library(stringr)

# Clear memory
rm(list=ls())

##############################################################
# Importing datasets
##############################################################

# load mammals taxonomy from https://www.mammaldiversity.org/ 
mam_tax <- read.csv('Data/MDD_v1.9_6596species.csv', stringsAsFactors = F)
mam_tax$family <- str_to_sentence(mam_tax$family) #convert to upper case
mam_tax$order <- str_to_sentence(mam_tax$order) #convert to upper case

# remove extinct species
mam_tax <- mam_tax[mam_tax$extinct != 1,]
mam_tax <- mam_tax[!mam_tax$iucnStatus %in% c("EW", "EX"),]

# remove marine species
mam_tax <- mam_tax[!mam_tax$order %in% c('Sirenia'),]
mam_tax <- mam_tax[!mam_tax$family %in% c('Octodontidae','Otariidae','Balaenidae','Balaenopteridae','Ziphiidae','Neobalaenidae','Delphinidae','Phocidae','Monodontidae','Dugongidae','Iniidae','Physeteridae','Phocoenidae','Odobenidae','Trichechidae', 'Platanistidae'),]

fam <- mam_tax %>% 
  group_by(family) %>%
  summarize(freq = n()) # 6345 species and 145 families

# generating list of species
mam_fam <- sort(unique(as.character(fam$family))) #145 families

##############################################################
# Formatting family data
##############################################################

# obtaining dataframe listing the Open Tree identifiers potentially 
# matching our list of families

taxa <- tnrs_match_names(mam_fam) 
taxa[taxa$approximate_match==TRUE,] # Prionodontidae, Heterocephalidae and Zenkerellidae are not well matched

# check all possible options that match
taxa[taxa$number_matches != 1,]
ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa, ott_id = ott_id_tocheck[i]))
}

# fix some of these families to older placements or subfamily names 
mam_tax[mam_tax$family =="Heterocephalidae", "family"] <- "Bathyergidae"
mam_tax[mam_tax$family =="Prionodontidae",  "family"] <- "Prionodontinae"
mam_tax[mam_tax$family =="Zenkerellidae",  "family"] <- "Anomaluridae"
mam_tax[mam_tax$family =="Zapodidae",  "family"] <- "Dipodidae"
mam_tax[mam_tax$family =="Sminthidae",  "family"] <- "Dipodidae"
mam_tax[mam_tax$family =="Potamogalidae",  "family"] <- "Tenrecidae"

fam[fam$family =="Heterocephalidae", "family"] <- "Bathyergidae"
fam[fam$family =="Prionodontidae",  "family"] <- "Prionodontinae"
fam[fam$family =="Zenkerellidae",  "family"] <- "Anomaluridae"
fam[fam$family =="Zapodidae",  "family"] <- "Dipodidae"
fam[fam$family =="Sminthidae",  "family"] <- "Dipodidae"
fam[fam$family =="Potamogalidae",  "family"] <- "Tenrecidae"

# mam_fam[mam_fam =="Heterocephalidae"] <- "Bathyergidae"
# mam_fam[mam_fam =="Prionodontidae"] <- "Prionodontinae"
# mam_fam[mam_fam =="Zenkerellidae"] <- "Anomaluridae"

# rerun again
mam_fam <- sort(unique(as.character(fam$family))) 

taxa.c <- tnrs_match_names(names = mam_fam)
taxa.c[taxa.c$approximate_match==TRUE,] # all good now

# exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa <- taxa.c
taxa[taxa$number_matches != 1,]
ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa, ott_id = ott_id_tocheck[i]))
}

# it's ok, they are synonyms

# check synonyms and change name accordingly
fix_taxa <- taxa[taxa$is_synonym == TRUE,] 

fix_taxa <-fix_taxa[,c("search_string", "unique_name")]
fix_taxa$family <- str_to_sentence(fix_taxa$search_string) #convert to upper case to join with original dataset

mam_tax <- left_join(mam_tax,fix_taxa, by =c("family" = "family"))
mam_tax$family <- ifelse(!is.na(mam_tax$unique_name), mam_tax$unique_name, mam_tax$family)

fam <- mam_tax %>% 
  group_by(family) %>%
  summarize(freq = n()) # 6345 species and 145 families

# generating list of families
mam_fam <- sort(unique(as.character(fam$family))) #145 families

# rerun 3
taxa.c2 <- tnrs_match_names(names = mam_fam)

taxa.c2[taxa.c2$approximate_match==TRUE,] # no species returned
taxa.c2[taxa.c2$is_synonym==TRUE,] # no species returned

##############################################################
# Retrieving phylogenetic relationships
##############################################################

# retrieving phylogenetic relationships among taxa in the form 
# of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa.c2[["ott_id"]], label_format = "name") # Zapodidae and Sminthidae not found
plot(tree, cex=.5, label.offset =.1, no.margin = TRUE)

##############################################################
# Dealing with polytomies
##############################################################

# we run the following function to check for the existence of polytomies.
# If polytomies exist, the output will be `FALSE`, and vice versa.

is.binary(tree) # there are some polytomies

# to take care of these polytomies, we are going to use a 
# randomization approach
set.seed(111) #making it replicable, at least for this version of R (i.e. v.4.0.2)
tree_random <- multi2di(tree,random=TRUE)
is.binary(tree_random)

##############################################################
# Final checks
##############################################################

# exploring whether our tree covers all the species we wanted 
# it to include, and making sure that the species names in our 
# database match those in the tree. We use the following code.

tree_random$tip.label <- gsub("_"," ", tree_random$tip.label)
intersect(as.character(tree_random$tip.label), as.character(species))

species[!species %in% as.character(tree_random$tip.label)] #listed in our database but not in the tree
tree_random$tip.label[!as.character(tree_random$tip.label) %in% species] # listed in the tree but not in our database

##I have a problem with a species labeled "mrcaott3599545ott4131616", and another labeled 
# "mrcaott80776ott602508" in the tree

## Also, "Anas cyanoptera"  "Loxia leucoptera" "Regulus regulus"  are in my data but not in the tree

# try to see which species is that

test<-tnrs_match_names(names = c("Anas cyanoptera", "Loxia leucoptera", "Regulus regulus" ))

tree_test <- tol_induced_subtree(ott_ids = c("4131616",  "602508", "82411"), label_format = "name")
tree_test <- tol_induced_subtree(ott_ids = c("4131616",  "602508", "206533"), label_format = "name")

# we fix them here
tree_random$tip.label[tree_random$tip.label =="mrcaott3599545ott4131616"] <-"Regulus regulus"
tree_random$tip.label[tree_random$tip.label =="mrcaott80776ott602508"] <-"Loxia leucoptera"
tree_random$tip.label[tree_random$tip.label =="mrcaott82415ott206533"] <-"Anas cyanoptera"

tiff("Results/bird_phylogenetic_tree_pruned.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# we can now save the tree
save(tree_random, file = "Data/bird_tree_random.Rdata")

##############################################################
# Computing branch lengths
##############################################################

# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we need to make sure that tree labels and database
# use the same nomenclature
setdiff(birddata$speciesname, as.character(tree_random$tip.label)) # Anas cyanoptera
setdiff(as.character(tree_random$tip.label),birddata$speciesname)

# exclude species in the tree that are not in dataset
drops <- tree_random$tip.label[!tree_random$tip.label %in% birddata$speciesname]
bird.tree_random.fixed <- drop.tip(tree_random, drops)

# exclude species in the dataset that are not in the tree
drop_sp <- birddata$speciesname[!birddata$speciesname %in% tree_random$tip.label]
birddata <- birddata[birddata$speciesname != drop_sp,] # 1545 sp

# save the new tree
write.tree(bird.tree_random.fixed, file = "Data/bird.tree_random.fixed.tre")

# compute branch lengths of tree
phylo_branch <- compute.brlen(bird.tree_random.fixed, method = "Grafen", power = 1)

# check if tree is ultrametric
is.ultrametric(phylo_branch) # TRUE


##############################################################
# Phylogenetic matrix
##############################################################

# matrix to be included in the models
bird_phylo_cor <- vcv(phylo_branch, cor = T)

# remove rows not in correlation matrix
birddata_ph <- birddata[which(birddata$speciesname %in% rownames(bird_phylo_cor)),] 

##create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID<-data.frame(speciesname = unique(birddata_ph$speciesname), SPID = paste0("SP",1:length(unique(birddata_ph$speciesname))))
SpID$speciesname<-as.character(SpID$speciesname)
birddata_ph<-inner_join(birddata_ph,SpID, by = "speciesname")

# finally, save matrix for future analyses
save(bird_phylo_cor, file = "Data/bird_phylo_cor.Rdata")


# exporting fixed dataset for analyses
write.csv(birddata_ph, 
          "Data/birddata_ph.csv", row.names = FALSE)

# End of script ####