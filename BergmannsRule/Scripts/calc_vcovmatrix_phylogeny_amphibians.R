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

# This script loads the pruned amphibian tree and calculates the phylogenetic vcov matrix to run 
# phylogenetic meta-analysis for amphibians for the paper:

# Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Uncovering the environmental drivers 
# of intraspecific body size variation in terrestrial vertebrates. 


##############################################################
# Packages needed
##############################################################

#load libraries
library(dplyr)
library(phytools)
#library(ape)
#library(treebase)
# library(rotl)
# library(diagram)

# Clear memory
rm(list=ls())

##############################################################
# Importing datasets
##############################################################

# Load data
# data
amph.tree <- read.tree("Data/amph_shl_dates.tre")
amphibians <- readRDS("Results/BergmannsRule_results_correlations_20211114.rds")
amphibians <- subset(amphibians,class=='amphibian')
amphibians$Species_ph <- gsub(" ", "_", trimws(amphibians$speciesname))

# generating list of species
species <- sort(unique(as.character(amphibians$speciesname))) # 36 species

# check that tree labels and database
# use the same nomenclature
setdiff(amphibians$Species_ph, as.character(amph.tree$tip.label)) # "Artibeus jamaicensis" "Geomys bursarius"     "Myotis lucifugus"     "Sturnira luisi" still missing
setdiff(as.character(amph.tree$tip.label),amphibians$Species_ph)


#exclude species in the tree that are not in your dataset
drops<-amph.tree$tip.label[!amph.tree$tip.label %in% amphibians$Species_ph]
amph.tree<-drop.tip(amph.tree, drops)

# compute branch lengths of tree 
# phylo_branch <- compute.brlen(amph.tree, method = "Grafen", power = 1)

# check if tree is ultrametric
is.ultrametric(amph.tree) # TRUE

##############################################################
# Phylogenetic matrix
##############################################################

# matrix to be included in the models
amph_phylo_cor <- vcv(amph.tree, cor = T)

# remove rows not in correlation matrix
amphibians_ph <- amphibians[which(amphibians$Species_ph %in% rownames(amph_phylo_cor)),] # we lose 30 species!
species_ph <- sort(unique(as.character(amphibians_ph$Species_ph))) # 24 species


##create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID<-data.frame(speciesname = unique(amphibians_ph$speciesname), SPID = paste0("SP",1:length(unique(amphibians_ph$speciesname))))
SpID$speciesname<-as.character(SpID$speciesname)
amphibians_ph<-inner_join(amphibians_ph,SpID, by = "speciesname")

# finally, save matrix for future analyses
save(amph_phylo_cor, file = "Data/amph_phylo_cor.Rdata")

# exporting fixed dataset for analyses
write.csv(amphibians_ph, 
          "Data/amphibians_ph.csv", row.names = FALSE)


# End of script ####