##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Email: ana.benitez@ugr.es

# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# www.anabenitezlopez.com
# Department of Environmental Science, Radboud University (the Netherlands)
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, Spain)
# Department of Zoology, University of Granada (UGR, Spain)

# First created on 31 May 2022

##############################################################
# Description of script and instructions                  ####
##############################################################

# Script to check family coverage in the dataset 

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
library(treeio)
library(wesanderson)
library(caper)

# Phylogenetic trees for the four taxa ####
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

mfam <- mam_tax %>% 
  group_by(Family.1.2) %>%
  summarize(freq = n()) # 5327 species

names(mfam) <- c("family", "freq")

# Load family-level tree
m_tree <- ape::read.nexus(file = "Data/Phylogeny/mam_tree_family.nex")

# load mammal dataset 
mamdata <- read.csv('Data/mamdata_ph.csv', stringsAsFactors = F)

# summarize by env.var and family 
mamdata2 <- mamdata %>% 
  group_by(env.var, family) %>%
  summarise(mean = mean(corr.coeff),
            nsp = n())

# use data for tavg only - discrete plot
mam_tavg <- mamdata2 %>% 
  filter(env.var == "tavg") 

# join families dataset with families included in Bergmann's rule paper
mfam <- left_join(mfam, mam_tavg, by = "family")
mfam$inc <- ifelse(is.na(mfam$mean), "no", "yes")
mfam$bin <- ifelse(mfam$inc == "yes", 1, 0)
mfam$bin <- as.integer(mfam$bin)
mfam$hyp <- ifelse(mfam$mean < 0, "yes","no")

# palette <- c("yes" = "dark purple", "no" = "yellow")

p <- ggtree(m_tree, layout = "circular") 
pm <- p %<+% mfam + geom_tree(aes(color = inc), size = 1) + geom_tiplab(aes(color=inc)) +
      scale_color_manual(values = wes_palette("Cavalcanti1", n = 2), breaks = c("no", "yes")) +
      theme(plot.margin = unit(c(3,3,3,3), "cm"))
pm

ggsave(filename = "Figures/phylo_mammals_included.png", 
       width= 250, height= 250, units = 'mm', dpi=300)

# exclude families in the tree that are not in dataset
mfam2 <- mfam[mfam$bin == 1,]
drops <- m_tree$tip.label[!m_tree$tip.label %in% mfam2$family]
m_tree2 <- drop.tip(m_tree, drops)

# create tree
p2 <-  ggtree(m_tree2, layout = "circular") 

# nsp per family
p2 %<+% mfam2 + geom_tree(aes(color = log10(nsp)), size = 1) + geom_tiplab(aes(color = log10(nsp))) +
  scale_color_viridis() 

# mean corr coeff tavg per family
p2 %<+% mfam2 + geom_tree(aes(color = mean), size = 1) + geom_tiplab(aes(color = mean)) +
  scale_color_viridis() 

# Adherence to heat conserv hypothesis birds (bergmann's rule)
p2 %<+% mfam2 + geom_tree(aes(color = hyp), size = 1) + geom_tiplab(aes(color = hyp)) +
  scale_color_manual(values= wes_palette("Cavalcanti1", n = 2)) + theme(legend.position = "none")

ggsave("Figures/phylo_mammals_tavg.png", 
       width= 250, height= 250, units = 'mm', dpi=300)


table(mfam2$hyp) # 26 out of 46 families follow bergmann's rule (56%)

# Test phylogenetic clustering of included mammalian families. Are most included families related?

# The D-statistic is the sum of state changes along the edges of a phylogeny. 
# The D statistic is equal to 1 if the observed binary trait has a phylogenetically random distribution across 
# the tips of the phylogeny and to 0 if the observed trait is as clumped as if it had evolved by Brownian motion under 
# our threshold model. Values of D can fall outside this range. 
# Continuous vertical lines represent the mean of random D-statistic values obtained from two null models: 
# Brownian motion (red) and phylogenetic randomness (blue). 
# Dashed gray lines represent observed D-statistics.
library(caper)

# order mfam
tips <- data.frame(family = m_tree$tip.label)
mfam <- left_join(tips, mfam, by = c("family"))

mfam_D <- comparative.data(m_tree, mfam[,c("family","bin")], family)
PhyloD <- phylo.d(mfam_D, binvar=bin) # Estimated D :  0.4084372, P no random = 0, P Browniam = 0.12
plot(PhyloD) # some phylogenetic clumping, but nothing too severe 

# order mfam
tips <- data.frame(family = m_tree2$tip.label)
mfam2 <- left_join(tips, mfam2, by = c("family"))

# Test phylogenetic clustering of families adhering to the heat conservation hypothesis. Are the majority of families that follow the heat conserv hyp clustered?
mfam_D2 <- comparative.data(m_tree2, mfam2[,c("family","hyp")], family) 
PhyloD2 <- phylo.d(mfam_D2, binvar = hyp) # Estimated D :  0.3867931, P no random: 0.078, P Brownian: 0.271
plot(PhyloD2) # phylogenetically random pattern, almost overdispersed across families.

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
b_tree <- ape::read.nexus(file = "Data/Phylogeny/Birdlife_Family_Phylogeny.nex")
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
pb <- p %<+% bfam + geom_tree(aes(color = inc), size = 1) + geom_tiplab(aes(color=inc)) +
  scale_color_manual(values= wes_palette("Cavalcanti1", n = 2)) + 
  theme(legend.position = "none", plot.margin = unit(c(3,3,3,3), "cm"))
pb

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
PhyloD <- phylo.d(bfam_D, binvar=bin) # Estimated D :  0.62, P no random: 0.006, P Brownian: 0.014
plot(PhyloD) # some phylogenetic clumping, but nothing too severe 

# Test phylogenetic clustering of families adhering to the heat conservation hypothesis. Are the majority of families that follow the heat conserv hyp clustered?
bfam_D2 <- comparative.data(b_tree3, bfam2[,c("family","hyp")], family) 
PhyloD2 <- phylo.d(bfam_D2, binvar = hyp)
plot(PhyloD2) # phylogenetically random pattern, almost overdispersed across families.

## reptiles ####
# load squamates taxonomy from Tonini et al. 2017
rept_tax <- read.csv('Data/squam_shl_new_Classification.csv', stringsAsFactors = F)

rfam <- rept_tax %>% 
  group_by(Family) %>%
  summarize(freq = n()) # 9755 species

# Load family-level tree - squamates
r_tree <- ape::read.nexus(file = "Data/Phylogeny/rept_tree_family.nex")

# load reptiles dataset - all squamates
reptdata <- read.csv('Data/reptdata_ph.csv', stringsAsFactors = F)

# summarize by env.var and family 
reptdata2 <- reptdata %>% 
  group_by(env.var, family) %>%
  summarise(mean = mean(corr.coeff),
            nsp = n())

# use data for npp only - discrete plot
rept_npp <- reptdata2 %>% 
  filter(env.var == "npp") 

colnames(rept_npp)[2] <- "Family"

# join families dataset with families included in Bergmann's rule paper
rfam <- left_join(rfam, rept_npp, by = "Family")
rfam$inc <- ifelse(is.na(rfam$mean), "no", "yes")
rfam$bin <- ifelse(rfam$inc == "yes", 1, 0)
rfam$bin <- as.integer(rfam$bin)
rfam$hyp <- ifelse(rfam$mean < 0, "yes","no")

# palette <- c("yes" = "dark purple", "no" = "yellow")

p <- ggtree(r_tree, layout = "circular") 
pr <- p %<+% rfam + geom_tree(aes(color = inc), size = 1) + geom_tiplab(aes(color=inc)) +
  scale_color_manual(values = wes_palette("Cavalcanti1", n = 2), breaks = c("no", "yes")) +
  theme(legend.position = "none", plot.margin = unit(c(3,3,3,3), "cm"))
pr

ggsave(filename = "Figures/phylo_reptiles_included.png", 
       width= 250, height= 250, units = 'mm', dpi=300)

# exclude families in the tree that are not in dataset
rfam2 <- rfam[rfam$bin == 1,]
drops <- r_tree$tip.label[!r_tree$tip.label %in% rfam2$Family]
r_tree2 <- drop.tip(r_tree, drops)

# create tree
p2 <-  ggtree(r_tree2, layout = "circular") 

# nsp per family
p2 %<+% rfam2 + geom_tree(aes(color = log10(nsp)), size = 1) + geom_tiplab(aes(color = log10(nsp))) +
  scale_color_viridis() 

# mean corr coeff tavg per family
p2 %<+% rfam2 + geom_tree(aes(color = mean), size = 1) + geom_tiplab(aes(color = mean)) +
  scale_color_viridis() 

# Adherence to resource availability hypothesis squamates 
p2 %<+% rfam2 + geom_tree(aes(color = hyp), size = 1) + geom_tiplab(aes(color = hyp)) +
  scale_color_manual(values= wes_palette("Cavalcanti1", n = 2)) + theme(legend.position = "none")

ggsave("Figures/phylo_reptiles_npp.png", 
       width= 250, height= 250, units = 'mm', dpi=300)

table(rfam2$hyp) # 9 out of 19 families follow resource avail. hyp (47%)

# Test phylogenetic clustering of included mammalian families. Are most included families related?

# The D-statistic is the sum of state changes along the edges of a phylogeny. 
# The D statistic is equal to 1 if the observed binary trait has a phylogenetically random distribution across 
# the tips of the phylogeny and to 0 if the observed trait is as clumped as if it had evolved by Brownian motion under 
# our threshold model. Values of D can fall outside this range. 
# Continuous vertical lines represent the mean of random D-statistic values obtained from two null models: 
# Brownian motion (red) and phylogenetic randomness (blue). 
# Dashed gray lines represent observed D-statistics.

# order rfam
tips <- data.frame(Family = r_tree$tip.label)
rfam <- left_join(tips, rfam, by = c("Family"))

rfam_D <- comparative.data(r_tree, rfam[,c("Family","bin")], Family)
PhyloD <- phylo.d(rfam_D, binvar=bin) # Estimated D :  0.8684831, P no random = 0.324, P Browniam = 0.05
plot(PhyloD) # some tendency for overdispersion

# order rfam
tips <- data.frame(Family = r_tree2$tip.label)
rfam2 <- left_join(tips, rfam2, by = c("Family"))

# Test phylogenetic clustering of families adhering to the resource availability hypothesis. Are the majority of families that follow the heat conserv hyp clustered?
rfam_D2 <- comparative.data(r_tree2, rfam2[,c("Family","hyp")], Family) 
PhyloD2 <- phylo.d(rfam_D2, binvar = hyp) # Estimated D :  0.377288, P no random: 0.184, P Brownian: 0.446
plot(PhyloD2) # phylogenetically clustering pattern.


## amphibians ####
# load data - Jetz and Pyron 2019
amph_tax <- read.csv('Data/amph_shl_new_Classification.csv', stringsAsFactors = F)
amph_tax <- amph_tax[-c(7239:7240),] # remove rows with no info or info on the Outgroup

afam <- amph_tax %>% 
  group_by(Family) %>%
  summarize(freq = n()) # 7238 species and 73 families

# Load family-level tree - squamates
a_tree <- ape::read.nexus(file = "Data/Phylogeny/amph_tree_family.nex")

# load amphibians dataset 
amphdata <- read.csv('Data/amphdata_ph.csv', stringsAsFactors = F)

# summarize by env.var and family 
amphdata2 <- amphdata %>% 
  group_by(env.var, family) %>%
  summarise(mean = mean(corr.coeff),
            nsp = n())

# use data for npp only - discrete plot
amph_npp <- amphdata2 %>% 
  filter(env.var == "npp") 

colnames(amph_npp)[2] <- "Family"

# join families dataset with families included in Bergmann's rule paper
afam <- left_join(afam, amph_npp, by = "Family")
afam$inc <- ifelse(is.na(afam$mean), "no", "yes")
afam$bin <- ifelse(afam$inc == "yes", 1, 0)
afam$bin <- as.integer(afam$bin)
afam$hyp <- ifelse(afam$mean < 0, "yes","no")

# palette <- c("yes" = "dark purple", "no" = "yellow")

p <- ggtree(a_tree, layout = "circular") 
pa <- p %<+% afam + geom_tree(aes(color = inc), size = 1) + geom_tiplab(aes(color=inc)) +
  scale_color_manual(values = wes_palette("Cavalcanti1", n = 2), breaks = c("no", "yes")) +
  theme(legend.position = "none", plot.margin = unit(c(3,3,3,3), "cm"))
pa

ggsave(filename = "Figures/phylo_amphibians_included.png", 
       width= 250, height= 250, units = 'mm', dpi=300)

# exclude families in the tree that are not in dataset
afam2 <- afam[afam$bin == 1,]
drops <- a_tree$tip.label[!a_tree$tip.label %in% afam2$Family]
a_tree2 <- drop.tip(a_tree, drops)

# create tree
p2 <-  ggtree(a_tree2, layout = "circular") 

# nsp per family
p2 %<+% afam2 + geom_tree(aes(color = log10(nsp)), size = 1) + geom_tiplab(aes(color = log10(nsp))) +
  scale_color_viridis() 

# mean corr coeff npp per family
p2 %<+% afam2 + geom_tree(aes(color = mean), size = 1) + geom_tiplab(aes(color = mean)) +
  scale_color_viridis() 

# Adherence to resource availability hypothesis amphibians
p2 %<+% afam2 + geom_tree(aes(color = hyp), size = 1) + geom_tiplab(aes(color = hyp)) +
  scale_color_manual(values= wes_palette("Cavalcanti1", n = 2)) + theme(legend.position = "none")

ggsave("Figures/phylo_amphibians_npp.png", 
       width= 250, height= 250, units = 'mm', dpi=300)

table(afam2$hyp) # 4 out of 11 families follow resource avail. hyp (36%)

# Test phylogenetic clustering of included mammalian families. Are most included families related?

# The D-statistic is the sum of state changes along the edges of a phylogeny. 
# The D statistic is equal to 1 if the observed binary trait has a phylogenetically random distribution across 
# the tips of the phylogeny and to 0 if the observed trait is as clumped as if it had evolved by Brownian motion under 
# our threshold model. Values of D can fall outside this range. 
# Continuous vertical lines represent the mean of random D-statistic values obtained from two null models: 
# Brownian motion (red) and phylogenetic randomness (blue). 
# Dashed gray lines represent observed D-statistics.

# order afam
tips <- data.frame(Family = a_tree$tip.label)
afam <- left_join(tips, afam, by = c("Family"))

afam_D <- comparative.data(a_tree, afam[,c("Family","bin")], Family)
PhyloD <- phylo.d(afam_D, binvar=bin) # Estimated D :  0.64918, P no random = 0.189, P Brownian = 0.197
plot(PhyloD) # families not phylogenetically clustered, no evidence from non random distribution either

# order afam
tips <- data.frame(Family = a_tree2$tip.label)
afam2 <- left_join(tips, afam2, by = c("Family"))

# Test phylogenetic clustering of families adhering to the resource availability hypothesis. Are the majority of families that follow the heat conserv hyp clustered?
afam_D2 <- comparative.data(a_tree2, afam2[,c("Family","hyp")], Family) 
PhyloD2 <- phylo.d(afam_D2, binvar = hyp) # Estimated D :  1.493811, P no random: 0.57, P Brownian: 0.23
plot(PhyloD2) # phylogenetically overdispersed pattern?? too small of a sample...

## Plot all phylogenies together into Fig. S1 ----
### Combine 4 paneles (pa,pr,pb,pm) into Fig. S1)
# modify plot. margins and legend size
pa2 <- pa + theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 18))
pr2 <- pr + theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 18))
pb2 <- pb + theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 18))
pm2 <- pm + theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"),
                  legend.title = element_blank(),
                        legend.text = element_text(size = 18))

ggarrange(pa2, pr2, pb2, pm2, ncol=2,nrow=2,
          labels= "auto",#hjust=-5,vjust=2,
          common.legend = T,
          legend = "bottom") 

### Save figure
ggsave(filename = 'Figures/Figure_S1.png', 
       width = 350, height = 350, units = 'mm', dpi=300)
