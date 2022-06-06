##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estacion Biologica de Donana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 15th of December 2021

##############################################################
# Description of script and instructions
##############################################################

# This script is used to build the phylogeny for birds and estimate the 
# phylogenetic relatedness among the species included in:

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

# load database
birddata<-read.csv("Data/birds.csv", header = TRUE, stringsAsFactors = FALSE) # 

# generating list of species
species <- sort(unique(as.character(birddata$speciesname))) #1561 species

##############################################################
# Formatting species data
##############################################################

# obtaining dataframe listing the Open Tree identifiers potentially 
# matching our list of species.

taxa <- tnrs_match_names(species) 

# 2 species not matched, we changed the names according to avibase
species[species=="Alcippe schaefferi"] <- "Alcippe davidi"
species[species=="Rhopospina fruticeti"] <- "Phrygilus fruticeti"

birddata[birddata$speciesname=="Alcippe schaefferi","speciesname"] <- "Alcippe davidi"
birddata[birddata$speciesname=="Rhopospina fruticeti","speciesname"] <- "Phrygilus fruticeti"

# rerun
taxa.c <- tnrs_match_names(names = species)
taxa.c[taxa.c$approximate_match==TRUE,] 

# according to the `approximate_match` column, there might be 
# 4 typos in the species list 

fix_taxa <- taxa.c[taxa.c$approximate_match==TRUE,] # 4 species returned

fix_taxa <-fix_taxa[,c("search_string", "unique_name")]
fix_taxa$species <- str_to_sentence(fix_taxa$search_string) #convert to upper case to join with original dataset

birddata <- left_join(birddata,fix_taxa, by =c("speciesname" = "species"))
birddata$speciesname <-ifelse(!is.na(birddata$unique_name), birddata$unique_name, birddata$speciesname)
birddata <- birddata[,-c(26:27)] # remove join columns

# rerun again
species <- sort(unique(as.character(birddata$speciesname))) #1546 species

taxa.c2 <- tnrs_match_names(names = species)
taxa.c2[taxa.c2$approximate_match==TRUE,] # all good now

# exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa <- taxa.c2
taxa[taxa$number_matches != 1,]
ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa, ott_id = ott_id_tocheck[i]))
}

# it's ok, they are synonyms

# check synonyms and change name accordingly
fix_taxa <- taxa[taxa$is_synonym == TRUE,] 

fix_taxa <-fix_taxa[,c("search_string", "unique_name")]
fix_taxa$species <- str_to_sentence(fix_taxa$search_string) #convert to upper case to join with original dataset

birddata <- left_join(birddata,fix_taxa, by =c("speciesname" = "species"))
birddata$speciesname <-ifelse(!is.na(birddata$unique_name), birddata$unique_name, birddata$speciesname)
birddata <- birddata[,-c(26:27)] # remove join columns

species <- sort(unique(as.character(birddata$speciesname))) #1546 species, we lose some species which were recorded as separate species, now they are duplicates

# rerun 3
taxa.c3 <- tnrs_match_names(names = species)

taxa.c3[taxa.c3$approximate_match==TRUE,] # no species returned
taxa.c3[taxa.c3$is_synonym==TRUE,] # no species returned

##############################################################
# Retrieving phylogenetic relationships
##############################################################

# retrieving phylogenetic relationships among taxa in the form 
# of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa.c3[["ott_id"]], label_format = "name")
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