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

# This script is used to build the phylogeny for reptiles and estimate the 
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

# load database and elton traits data to remove marine species
reptdata<-read.csv("Data/reptiles.csv", header = TRUE, stringsAsFactors = FALSE) # 

# generating list of species
species <- sort(unique(as.character(reptdata$speciesname))) #81 species

##############################################################
# Formatting species data
##############################################################

# obtaining dataframe listing the Open Tree identifiers potentially 
# matching our list of species.

taxa <- tnrs_match_names(species) 

# according to the `approximate_match` column, there might be 
# 0 typos in the species list 
# nrow(taxa[taxa$approximate_match==TRUE,])
taxa[taxa$approximate_match==TRUE,] # no species returned

# fixing those typos (example in case they were unresolved matches)
#species[species=="XX"] <- "XX"

#reptdata$speciesname<-as.character(reptdata$speciesname)
#reptdata[reptdata$speciesname=="XX","speciesname"] <- "XX"

# rerun
# taxa.c <- tnrs_match_names(names = species)
# 
# taxa.c[taxa.c$approximate_match==TRUE,] # no species returned

# exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa.c <- taxa
taxa.c[taxa.c$number_matches != 1,]
ott_id_tocheck <- taxa.c[taxa.c$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
   print(inspect(taxa.c, ott_id = ott_id_tocheck[i]))
 }

# No multiple matches

# check synonyms and change name accordingly
fix_taxa <- taxa.c[taxa.c$is_synonym==TRUE,] # 10 species returned

fix_taxa <-fix_taxa[,c("search_string", "unique_name")]
fix_taxa$species <- str_to_sentence(fix_taxa$search_string) #convert to upper case to join with original dataset

reptdata <- left_join(reptdata,fix_taxa, by =c("speciesname" = "species"))
reptdata$speciesname <-ifelse(!is.na(reptdata$unique_name), reptdata$unique_name, reptdata$speciesname)
reptdata <- reptdata[,-c(10:11)] # remove join columns

species <- sort(unique(as.character(reptdata$speciesname))) #81 species

# rerun 2
taxa.c2 <- tnrs_match_names(names = species)

taxa.c2[taxa.c2$approximate_match==TRUE,] # no species returned
taxa.c2[taxa.c2$is_synonym==TRUE,] # no species returned


##############################################################
# Retrieving phylogenetic relationships
##############################################################

# retrieving phylogenetic relationships among taxa in the form 
# of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa.c2[["ott_id"]], label_format = "name")
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


tiff("Results/rept_phylogenetic_tree_pruned.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# we can now save the tree
save(tree_random, file = "Data/rept_tree_random.Rdata")


##############################################################
# Computing branch lengths
##############################################################

# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we need to make sure that tree labels and database
# use the same nomenclature
setdiff(reptdata$speciesname, as.character(tree_random$tip.label)) # 0 species
setdiff(as.character(tree_random$tip.label),reptdata$speciesname)

# exclude species in the tree that are not in dataset
drops <- tree_random$tip.label[!tree_random$tip.label %in% reptdata$speciesname]
rept.tree_random.fixed <- drop.tip(tree_random, drops)

# save the new tree
write.tree(rept.tree_random.fixed, file = "Data/rept.tree_random.fixed.tre")

# compute branch lengths of tree
phylo_branch <- compute.brlen(rept.tree_random.fixed, method = "Grafen", power = 1)

# check if tree is ultrametric
is.ultrametric(phylo_branch) # TRUE


##############################################################
# Phylogenetic matrix
##############################################################

# matrix to be included in the models
rept_phylo_cor <- vcv(phylo_branch, cor = T)

# remove rows not in correlation matrix
reptdata_ph <- reptdata[which(reptdata$speciesname %in% rownames(rept_phylo_cor)),] # we do not lose any species

##create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID<-data.frame(speciesname = unique(reptdata_ph$speciesname), SPID = paste0("SP",1:length(unique(reptdata_ph$speciesname))))
SpID$speciesname<-as.character(SpID$speciesname)
reptdata_ph<-inner_join(reptdata_ph,SpID, by = "speciesname")

# finally, save matrix for future analyses
save(rept_phylo_cor, file = "Data/rept_phylo_cor.Rdata")

# exporting fixed dataset for analyses
write.csv(reptdata_ph, 
          "Data/reptdata_ph.csv", row.names = FALSE)

# End of script ####