##############################################################
# Authors: 
# Ana Benítez
# Email: abenitez81@gmail.com

# Script to check family coverage in the dataset 
# Made on 31 May 2022

# clean environment
rm(list = ls())

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

# Phylogenetic trees for the four taxa ####
# mammals ####
# load mammals data from PHYLACINE
mam_tax <- read.csv('Data/PHYLACINE_Trait_data.csv', stringsAsFactors = F)

fam <- mam_tax %>% 
  group_by(Family.1.2) %>%
  summarize(freq = n()) # 5831 species

# load tree - Tobias et al. 2022 - AVONET - BirdLife
# m_tree <- ape::read.nexus(file = "Data/Phylogeny/mam_complete_phylogeny.nex") # 1000 trees
# m_tree
# plot(m_tree)

m_tree_sm <- ape::read.nexus(file = "Data/Phylogeny/mam_small_phylogeny.nex")
class(m_tree_sm[1]) #multiphylo
plot(m_tree_sm[])



# load mammal dataset 
mamdata <- read.csv('Data/mamdata_ph.csv', stringsAsFactors = F)
mamdata <- mamdata[mamdata$env.var == 'tavg', ]



ET2<-ET[!duplicated(ET$MSWFamilyLatin),] # 132 families
mamdata2<-mamdata[!duplicated(mamdata$family),] # 49 families 49/132*100 = 37.1 %



## birds ####
# load AVONET Datase BirdLife taxonomy
avonet<-read.csv('Data/AVONET1_BirdLife.csv', stringsAsFactors = F)

# habfam <- avonet %>% 
#           group_by(Habitat, Family1) %>%
#           summarize(freq = n())

# marine <- habfam %>% 
#           filter(Habitat == "Marine")

fam <- avonet %>% 
      group_by(Family1) %>%
      summarize(freq = n()) # 11009 species

# load tree - Tobias et al. 2022 - AVONET - BirdLife
b_tree <- ape::read.nexus(file = "Data/Birdlife_Family_Phylogeny.nex")
b_tree
plot(b_tree)

# load data bergmann
birddata <- read.csv('Data/birddata_ph.csv', stringsAsFactors = F)

# summarize by env.var and family 
birddata2 <- birddata %>% 
  group_by(env.var, family) %>%
  summarise(mean = mean(corr.coeff),
            nsp = n())

# use data for tavg only - discrete plot
bird_tavg <- birddata2 %>% 
  filter(env.var == "tavg") # 1545 species

# create dataframe with all described families 
bfam <- data.frame(family = b_tree$tip.label)

# remove families with marine species from AVONET and tree
bfam <- bfam[!bfam$family %in% c("Alcidae","Chionidae", "Diomedeidae", "Fregatidae", "Hydrobatidae", "Spheniscidae", 'Pelecanidae', "Procellariidae", "Phaethontidae"),]
bfam <- data.frame(family = bfam)
b_tree2 <- drop.tip(b_tree, tip = c("Alcidae","Chionidae", "Diomedeidae", "Fregatidae", "Hydrobatidae", "Spheniscidae", 'Pelecanidae', "Procellariidae", "Phaethontidae"))

# join dataset with families included in Bergmann's rule paper
bfam <- left_join(bfam, bird_tavg, by = "family")
bfam$inc <- ifelse(is.na(bfam$mean), "no", "yes")
bfam$bin <- ifelse(bfam$inc == "yes", 1, 0)
bfam$bin <- as.integer(bfam$bin)
bfam$hyp <- ifelse(bfam$mean < 0, "yes","no")

# palette <- c("yes" = "dark purple", "no" = "yellow")

p <- ggtree(b_tree2, layout = "circular") 
p %<+% bfam + geom_tree(aes(color = inc), size = 1) + geom_tiplab(aes(color=inc)) +
  scale_color_manual(values= wes_palette("Cavalcanti1", n = 2)) 

ggsave("Figures/phylo_birds.png", 
       width= 250, height= 250, units = 'mm', dpi=300)


# exclude families in the tree that are not in dataset
bfam2 <- bfam[bfam$bin == 1,]
drops <- b_tree2$tip.label[!b_tree2$tip.label %in% bfam2$family]
b_tree3 <- drop.tip(b_tree2, drops)

# create tree
p2 <-  ggtree(b_tree3, layout = "circular") 

# nsp per family
p2 %<+% bfam2 + geom_tree(aes(color = log10(nsp)), size = 1) + geom_tiplab(aes(color = log10(nsp))) +
        scale_color_viridis() 

# mean corr coeff tavg per family
p2 %<+% bfam2 + geom_tree(aes(color = mean), size = 1) + geom_tiplab(aes(color = mean)) +
  scale_color_viridis() 

# Adherence to heat conserv hypothesis birds (bergmann's rule)
p2 %<+% bfam2 + geom_tree(aes(color = hyp), size = 1) + geom_tiplab(aes(color = hyp)) +
         scale_color_manual(values= wes_palette("Cavalcanti1", n = 2)) 

ggsave("Figures/phylo_birds_tavg.png", 
       width= 250, height= 250, units = 'mm', dpi=300)


table(bfam2$hyp) # 104 families follow bergmann's rule

# Test phylogenetic clustering of included avian families. Are most included families related?

# The D-statistic is the sum of state changes along the edges of a phylogeny. 
# The D statistic is equal to 1 if the observed binary trait has a phylogenetically random distribution across 
# the tips of the phylogeny and to 0 if the observed trait is as clumped as if it had evolved by Brownian motion under 
# our threshold model. Values of D can fall outside this range. 
# Continuous vertical lines represent the mean of random D-statistic values obtained from two null models: 
# Brownian motion (red) and phylogenetic randomness (blue). 
# Dashed gray lines represent observed D-statistics.

bfam_D <- comparative.data(b_tree2, bfam[,c("family","bin")], family)
PhyloD <- phylo.d(bfam_D, binvar=bin)
plot(PhyloD) # some phylogenetic clumping, but nothing too severe 

# Test phylogenetic clustering of families adhering to the heat conservation hypothesis. Are the majority of families that follow the heat conserv hyp clustered?
bfam_D2 <- comparative.data(b_tree3, bfam2[,c("family","hyp")], family) 
PhyloD2 <- phylo.d(bfam_D2, binvar = hyp)
plot(PhyloD2) # phylogenetically random pattern, almost overdispersed across families.

## reptiles ####
rept_tax <- read.csv('Data/reptile_checklist_2022_03.csv', stringsAsFactors = F)

reptdata <- read.csv('Data/reptdata_ph.csv', stringsAsFactors = F)
reptdata <- reptdata[reptdata$env.var == 'tavg', ]

ET2<-ET[!duplicated(ET$BLFamilyLatin),] # 184 families
reptdata2<-reptdata[!duplicated(reptdata$family),] # 140 families (140/184*100 = 76%)

r_tree <- ape::read.tree(file = "Data/rep_phyl_Pyron.txt")
r_tree 

plot(r_tree)




# plotBranchbyTrait(b_tree, bfam$mean, # x2 variable used to colour

#                   type = "fan", edge.width=0.5, cex=0.4, show.tip.label=FALSE,
#                   mode='tips', 
#                   palette=palette, legend=TRUE)#viridis




## amphibians ####




