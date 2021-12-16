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

# This script is used to build the phylogeny for mammals and estimate the 
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
mamdata<-read.csv("Data/mammals.csv", header = TRUE, stringsAsFactors = FALSE) # 596

# remove unknown species(genus level)
mamdata <- mamdata[mamdata$speciesname != "Carollia carollia", ]
mamdata <- mamdata[mamdata$speciesname != "Glossophaga glossophaga", ]
mamdata <- mamdata[mamdata$speciesname != "Oryzomys oryzomys", ]
mamdata <- mamdata[mamdata$speciesname != "Reithrodontomys reithrodontomys", ]
mamdata <- mamdata[mamdata$speciesname != "Sturnira sturnira", ]
mamdata <- mamdata[mamdata$speciesname != "Lophostoma aequatorialis", ]
mamdata <- mamdata[mamdata$speciesname != "Peromyscus peromyscus", ]

# generating list of species
species <- sort(unique(as.character(mamdata$speciesname))) #556 species

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
#species[species=="Crocidura attenuatta"] <- "Crocidura attenuata"

#mamdata$speciesname<-as.character(mamdata$speciesname)
#mamdata[mamdata$speciesname=="Crocidura attenuatta","speciesname"] <- "Crocidura attenuata"

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

# check synonyms and change name accordingly
fix_taxa <- taxa.c[taxa.c$is_synonym==TRUE,] # 9 species returned

fix_taxa <-fix_taxa[,c("search_string", "unique_name")]
fix_taxa$species <- str_to_sentence(fix_taxa$search_string) #convert to upper case to join with original dataset

mamdata <- left_join(mamdata,fix_taxa, by =c("speciesname" = "species"))
mamdata$speciesname <-ifelse(!is.na(mamdata$unique_name), mamdata$unique_name, mamdata$speciesname)
mamdata <- mamdata[,-c(10:12)] # remove join columns

species <- sort(unique(as.character(mamdata$speciesname))) #556 species

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

##I have a problem with a species labeled "mrcaott319315ott366155", and another labeled 
# "mrcaott62482ott62486" in the tree

## Also, "Artibeus jamaicensis"     "Artibeus planirostris"    "Geomys bursarius" 
# "Miniopterus schreibersii" "Myotis lucifugus" "Sturnira luisi" are in my data but not in the tree

# try to see which species is that

test<-tnrs_match_names(names = c("Artibeus jamaicensis","Artibeus planirostris","Geomys bursarius",
                                 "Miniopterus schreibersii", "Myotis lucifugus", "Sturnira luisi"))

tree_test <- tol_induced_subtree(ott_ids = c("366155",  "62486"), label_format = "name")

# Miniopterus schreibersii is mrcaott319315ott366155 and Artibeus planirostris is 62486
# we fix them here
tree_random$tip.label[tree_random$tip.label =="mrcaott319315ott366155"] <-"Miniopterus schreibersii"
tree_random$tip.label[tree_random$tip.label =="mrcaott62482ott62486"] <-"Artibeus planirostris"

tiff("Results/mam_phylogenetic_tree_pruned.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# we can now save the tree
save(tree_random, file = "Data/mam_tree_random.Rdata")

##############################################################
# Computing branch lengths
##############################################################

# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we need to make sure that tree labels and database
# use the same nomenclature
setdiff(mamdata$speciesname, as.character(tree_random$tip.label)) # "Artibeus jamaicensis" "Geomys bursarius"     "Myotis lucifugus"     "Sturnira luisi" still missing
setdiff(as.character(tree_random$tip.label),mamdata$speciesname)

# exclude species in the tree that are not in dataset
  drops <- tree_random$tip.label[!tree_random$tip.label %in% mamdata$speciesname]
  mam.tree_random.fixed <- drop.tip(tree_random, drops)

# save the new tree
write.tree(mam.tree_random.fixed, file = "Data/mam.tree_random.fixed.tre")

# compute branch lengths of tree
phylo_branch <- compute.brlen(mam.tree_random.fixed, method = "Grafen", power = 1)

# check if tree is ultrametric
is.ultrametric(phylo_branch) # TRUE


##############################################################
# Phylogenetic matrix
##############################################################

# matrix to be included in the models
mam_phylo_cor <- vcv(phylo_branch, cor = T)

# remove rows not in correlation matrix
mamdata_ph <- mamdata[which(mamdata$speciesname %in% rownames(mam_phylo_cor)),] # 554

##create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID<-data.frame(speciesname = unique(mamdata_ph$speciesname), SPID = paste0("SP",1:length(unique(mamdata_ph$speciesname))))
SpID$speciesname<-as.character(SpID$speciesname)
mamdata_ph<-inner_join(mamdata_ph,SpID, by = "speciesname")

# finally, save matrix for future analyses
save(mam_phylo_cor, file = "Data/mam_phylo_cor.Rdata")


# exporting fixed dataset for analyses
write.csv(mamdata_ph, 
           "Data/mamdata_ph.csv", row.names = FALSE)

# End of script ####