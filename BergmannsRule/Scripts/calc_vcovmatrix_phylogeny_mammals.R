##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estacion Biologica de Donana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 12th of December 2021

##############################################################
# Description of script and instructions
##############################################################

# This script loads the pruned mammal tree and calculates the phylogenetic vcov matrix to run 
# phylogenetic meta-analysis for the size-NPP correlations for the paper:

# Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Uncovering the environmental drivers 
# of intraspecific body size variation in terrestrial vertebrates. 


##############################################################
# Packages needed
##############################################################

#load libraries
library(dplyr)
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
 mammals <- read.csv("Data/mammals.csv", stringsAsFactors = F)
 
# generating list of species
 # species <- sort(unique(as.character(mammals$speciesname))) # 596 species

# load pruned tree
 mam.tree <- read.tree("Data/mam.tree.tre")

# compute branch lengths of tree 
# phylo_branch <- compute.brlen(mam.tree, method = "Grafen", power = 1)

# check if tree is ultrametric
is.ultrametric(mam.tree) # TRUE

##############################################################
# Phylogenetic matrix
##############################################################

# matrix to be included in the models
mam_phylo_cor <- vcv(mam.tree, cor = T)

# remove rows not in correlation matrix
mammals_ph <- mammals[which(mammals$Species_ph %in% rownames(mam_phylo_cor)),] # we lose 30 species!

##create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID<-data.frame(speciesname = unique(mammals_ph$speciesname), SPID = paste0("SP",1:length(unique(mammals_ph$speciesname))))
SpID$speciesname<-as.character(SpID$speciesname)
mammals_ph<-inner_join(mammals_ph,SpID, by = "speciesname")

# finally, save matrix for future analyses
save(mam_phylo_cor, file = "Data/mam_phylo_cor.Rdata")

# exporting fixed dataset for analyses
write.csv(mammals_ph, 
           "Data/mammals_ph.csv", row.names = FALSE)


# End of script ####