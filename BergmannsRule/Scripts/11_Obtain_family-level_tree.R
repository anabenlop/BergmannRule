# Pruning trees to one member per genus, 
# Liam Revell solution
# http://blog.phytools.org/2014/11/pruning-trees-to-one-member-per-genus.html

# load libraries
library(tidyverse)
library(viridis)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(scales)
library(dplyr)
library(ape)
library(ggtree)
library(phytools)
library(TreeTools)
library(treeio)
library(wesanderson)

rm(list = ls())

## mammals ####
# load mammals taxonomy from PHYLACINE
mam_tax <- read.csv('Data/PHYLACINE_Trait_data.csv', stringsAsFactors = F)

# remove extinct species
mam_tax <- mam_tax[!mam_tax$IUCN.Status.1.2 %in% c("EW", "EX", "EP"),]

# remove marine species
mam_tax <- mam_tax[mam_tax$Marine != 1,]
mam_tax <- mam_tax[!mam_tax$Order.1.2 %in% c('Sirenia'),]
mam_tax <- mam_tax[!mam_tax$Family.1.2 %in% c('Octodontidae','Otariidae','Balaenidae','Balaenopteridae','Ziphiidae','Neobalaenidae','Delphinidae','Phocidae','Monodontidae','Dugongidae','Iniidae','Physeteridae','Phocoenidae','Odobenidae','Trichechidae', 'Platanistidae'),]

# fix CingulataFam cases
mam_tax$Family.1.2 <- ifelse(mam_tax$Family.1.2 == "CingulataFam", "Chlamyphoridae", mam_tax$Family.1.2) 
mam_tax$Family.1.2 <- ifelse(mam_tax$Genus.1.2 == "Dasypus", "Dasypodidae", mam_tax$Family.1.2) 

fam <- mam_tax %>% 
  group_by(Family.1.2) %>%
  summarize(Species = first(Binomial.1.2),
            freq = n()) # 5327 species

# load species-level tree from PHYLACINE
m_tree_1000 <- ape::read.nexus(file = "Data/Phylogeny/mam_complete_phylogeny.nex") # not included in repo as it is too large
# class(m_tree_1000[[1]]) #multiphylo
m_tree <- m_tree_1000[[1]] # pick first tree, we are not concerned about resolution as we are interested in families. using other trees does not change the final tree
rm(m_tree_1000)

# prune tree by excluding all species but one per family 
setdiff(fam$Species, as.character(m_tree$tip.label)) # check differences
# setdiff(as.character(m_tree$tip.label),fam$Species)

drops <- m_tree$tip.label[!m_tree$tip.label %in% fam$Species]
m_tree2 <- drop.tip(m_tree, drops)

# change tip.labels to family names
tips <- data.frame(Species = m_tree2$tip.label)
tips_df <- left_join(tips, fam, by = c("Species"))
m_tree2$tip.label <- c(tips_df$Family.1.2)

write.nexus(m_tree2, file = "Data/Phylogeny/mam_tree_family.nex")

p <- ggtree(m_tree2, layout = "circular") 
p %<+% fam + geom_tree(size = 1) + geom_tiplab() 

# ggsave("Figures/phylo_mam_test.png", 
#        width= 250, height= 250, units = 'mm', dpi=300)
# 

## reptiles ####
# load squamates taxonomy from Tonini et al. 2017
rept_tax <- read.csv('Data/squam_shl_new_Classification.csv', stringsAsFactors = F)

fam <- rept_tax %>% 
  group_by(Family) %>%
  summarize(Species = first(Species), 
            freq = n()) # 9755 species, and 73 families

fam$Species <- sub(" ", "_", fam$Species)

# Load species-level tree - squamates (Tonini et al. 2017)
r_tree <- ape::read.tree(file = "Data/Phylogeny/squam_shl_new_Consensus_9755.tre")

# prune tree by excluding all species but one per family 
setdiff(fam$Species, as.character(r_tree$tip.label)) # check differences
# setdiff(as.character(m_tree$tip.label),fam$Species)

drops <- r_tree$tip.label[!r_tree$tip.label %in% fam$Species]
r_tree2 <- drop.tip(r_tree, drops)

# change tip.labels to family names
tips <- data.frame(Species = r_tree2$tip.label)
tips_df <- left_join(tips, fam, by = c("Species"))
r_tree2$tip.label <- c(tips_df$Family)

write.nexus(r_tree2, file = "Data/Phylogeny/rept_tree_family.nex")

p <- ggtree(r_tree2, layout = "circular") 
p %<+% fam + geom_tree(size = 1) + geom_tiplab() 

## amphibians ####
# load data - Jetz and Pyron 2019
amph_tax <- read.csv('Data/amph_shl_new_Classification.csv', stringsAsFactors = F)
amph_tax <- amph_tax[-c(7239:7240),] # remove rows with no info or info on the Outgroup

fam <- amph_tax %>% 
  group_by(Family) %>%
  summarize(Species = first(Scientific.Name), 
            freq = n()) # 75 families 7238 species

fam$Species <- sub(" ", "_", fam$Species)

# Load species-level tree  (Jetz and Pyron 2019)
a_tree <- ape::read.tree(file = "Data/Phylogeny/amph_shl_new_Consensus_7238.txt")

# prune tree by excluding all species but one per family 
setdiff(fam$Species, as.character(a_tree$tip.label)) # check differences
# setdiff(as.character(m_tree$tip.label),fam$Species)

drops <- a_tree$tip.label[!a_tree$tip.label %in% fam$Species]
a_tree2 <- drop.tip(a_tree, drops)

# change tip.labels to family names
tips <- data.frame(Species = a_tree2$tip.label)
tips_df <- left_join(tips, fam, by = c("Species"))
a_tree2$tip.label <- c(tips_df$Family)

write.nexus(a_tree2, file = "Data/Phylogeny/amph_tree_family.nex")

p <- ggtree(a_tree2, layout = "circular") 
p %<+% fam + geom_tree(size = 1) + geom_tiplab() 

# End of script ----

