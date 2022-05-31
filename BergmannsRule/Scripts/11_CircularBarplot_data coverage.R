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
ET<-read.csv('Data/EltonTraits_Mammals_taxid.csv', stringsAsFactors = F)
ET<-ET[!ET$Order %in% c('Sirenia'),]
ET<-ET[!ET$MSWFamilyLatin %in% c('Octodontidae','Otariidae','Balaenidae','Balaenopteridae','Ziphiidae','Neobalaenidae','Delphinidae','Phocidae','Monodontidae','Dugongidae','Iniidae','Physeteridae','Phocoenidae','Odobenidae','Trichechidae', 'Platanistidae'),]

mamdata <- read.csv('Data/mamdata_ph.csv', stringsAsFactors = F)
mamdata <- mamdata[mamdata$env.var == 'tavg', ]

# ET$Order<-as.character(ET$Order)
# #ET$Order[ET$Order %in% c('Proboscidea', 'Perissodactyla','Cetartiodactyla')]<-'Ungulata'
# DT$Order<-as.character(DT$Order)
# #DT$Order[DT$Order %in% c('Proboscidea', 'Perissodactyla','Cetartiodactyla')]<-'Ungulata'
# D$Order<-as.character(D$Order)
# #D$Order[D$Order %in% c('Proboscidea', 'Perissodactyla','Cetartiodactyla')]<-'Ungulata'

ET2<-ET[!duplicated(ET$MSWFamilyLatin),] # 132 families
mamdata2<-mamdata[!duplicated(mamdata$family),] # 49 families 49/132*100 = 37.1 %

TabET<-table(ET$Order) #n species
TabETf<-table(ET2$Order) #n families
TabSP<-table(mamdata$order) #n species
TabSPf<-table(mamdata2$order) #n families

ET<-data.frame(Order=names(TabET), Nfm=as.numeric(TabETf), Nsps=as.numeric(TabET))
mam_dt<-data.frame(Order=names(TabSP), NfmInSamp=as.numeric(TabSPf), NspInSamp= as.numeric(TabSP))

data<-left_join(ET, mam_dt, by = 'Order')
data$NspInSamp[is.na(data$NspInSamp)]<-0
data$NfmInSamp[is.na(data$NfmInSamp)]<-0

data$Percentage<-round(100*data$NspInSamp/data$Nsps)
data$PercentageMissing<-100-data$Percentage
data$PercentageFm<-round(100*data$NfmInSamp/data$Nfm)
data$PercentageFm<-ifelse(data$PercentageFm>100, 100, data$PercentageFm)
data$PercentageMissingFm<-100-data$PercentageFm

data2 <- data

#PERCENTAGE OF FAMILIES
data2$labelFam<-paste0(data2$Order)

data2 <- as.data.frame(pivot_longer(data, cols=c(PercentageFm, PercentageMissingFm), names_to = "Var", values_to = "Value"))

data2$id<-as.numeric(as.factor(data2$Order))
data2$angle= 90 - 360 * (data2$id-0.5) / length(unique(data2$id))     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
data2$Var<-factor(data2$Var, levels=c('PercentageMissingFm','PercentageFm'))
data2$hjust<-ifelse(data2$angle < -90, 1, 0)
data2$angle<-ifelse(data2$angle < -90, data2$angle+180, data2$angle)
data2$LabNum<-ifelse(data2$Var =='PercentageFm', data2$Value, '')
data2$NfmInSamp<-ifelse(data2$NfmInSamp>data2$Nfm, data2$Nfm, data2$NfmInSamp)
data2$Lab2<-ifelse(data2$Var =='PercentageFm', paste0(' ', data2$LabNum, '% (',data2$NfmInSamp,'/',data2$Nfm,') '), '')
#data2[data2$Var=='Nsps', 'Value']<-100
#data2<-data2[data2$Var!='NspInSamp',]

# Create dataset
p1<-ggplot(data2) +      
  # Add the stacked bar
  geom_bar(aes(x=Order, y=Value, fill=Var), stat="identity", alpha=0.5) +
 # scale_fill_viridis(discrete=TRUE, direction=-1) +  
  scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(1,0,0,1), "cm")      
  ) +
  coord_polar(clip='off') +
  # Add labels on top of each bar
  geom_text(aes(x=id, y=40, label=Order, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4, angle= data2$angle, inherit.aes = FALSE) +
  geom_text(aes(x=id, y=100, label=Lab2, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4.2, angle= data2$angle, inherit.aes = FALSE )
p1

#PERCENTAGE OF SPECIES
data2 <- as.data.frame(pivot_longer(data, cols=c(Percentage, PercentageMissing), names_to = "Var", values_to = "Value"))

data2$id<-as.numeric(as.factor(data2$Order))
data2$angle= 90 - 360 * (data2$id-0.5) / length(unique(data2$id))     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
data2$Var<-factor(data2$Var, levels=c('PercentageMissing','Percentage'))
data2$hjust<-ifelse(data2$angle < -90, 1, 0)
data2$angle<-ifelse(data2$angle < -90, data2$angle+180, data2$angle)

data2$LabNum<-ifelse(data2$Var =='Percentage', data2$Value, '')
data2$Lab2<-ifelse(data2$Var =='Percentage', paste0(' ',data2$LabNum, '% (',data2$NspInSamp,'/',data2$Nsps,') '), '')

# Create dataset
p2<-ggplot(data2) +      
  # Add the stacked bar
  geom_bar(aes(x=Order, y=Value, fill=Var), stat="identity", alpha=0.5) +
  #scale_fill_viridis(discrete=TRUE, direction=-1) +  
  scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(1,0,0,1), "cm")      
  ) +
  coord_polar(clip='off') +
  # Add labels on top of each bar
  geom_text(aes(x=id, y=40, label=Order, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4, angle= data2$angle, inherit.aes = FALSE) +
  geom_text(aes(x=id, y=100, label=Lab2, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4.2, angle= data2$angle, inherit.aes = FALSE )
p2

#SAVE ALL
pdf('Figures/FamilyCoverage.pdf', width=18, height=12)#width=13, height=8)
ggarrange(p2, p1, labels=c('(a)', '(b)'), font.label=list(size=20, face='plain'),
          nrow=2, ncol=2)
dev.off()

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

table(bfam2$bergmann) # 104 families follow bergmann's rule


# phylogenetic clustering

bfam$bin <- ifelse(bfam$inc == "yes", 1, 0)
bfam$bin <- as.integer(bfam$bin)

bfam_D <- comparative.data(b_tree2, bfam, family) # does not work, it only shows families included in the analysis.
PhyloD <- phylo.d(bfam_D, binvar=bin)

# reptiles ####
ET<-read.csv('Data/reptile_checklist_2022_03.csv', stringsAsFactors = F)

reptdata <- read.csv('Data/reptdata_ph.csv', stringsAsFactors = F)
reptdata <- reptdata[reptdata$env.var == 'tavg', ]

ET2<-ET[!duplicated(ET$BLFamilyLatin),] # 184 families
reptdata2<-reptdata[!duplicated(reptdata$family),] # 140 families (140/184*100 = 76%)

TabET<-table(ET$IOCOrder) #n species
TabETf<-table(ET2$IOCOrder) #n families
TabSP<-table(reptdata$order) #n species
TabSPf<-table(reptdata2$order) #n families

ET<-data.frame(Order=names(TabET), Nfm=as.numeric(TabETf), Nsps=as.numeric(TabET))
rept_dt<-data.frame(Order=names(TabSP), NfmInSamp=as.numeric(TabSPf), NspInSamp= as.numeric(TabSP))

data<-left_join(ET, rept_dt, by = 'Order')
data$NspInSamp[is.na(data$NspInSamp)]<-0
data$NfmInSamp[is.na(data$NfmInSamp)]<-0

data$Percentage<-round(100*data$NspInSamp/data$Nsps)
data$PercentageMissing<-100-data$Percentage
data$PercentageFm<-round(100*data$NfmInSamp/data$Nfm)
data$PercentageFm<-ifelse(data$PercentageFm>100, 100, data$PercentageFm)
data$PercentageMissingFm<-100-data$PercentageFm

data2 <- data

#PERCENTAGE OF FAMILIES
data2$labelFam<-paste0(data2$Order)

data2 <- as.data.frame(pivot_longer(data, cols=c(PercentageFm, PercentageMissingFm), names_to = "Var", values_to = "Value"))

data2$id<-as.numeric(as.factor(data2$Order))
data2$angle= 90 - 360 * (data2$id-0.5) / length(unique(data2$id))     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
data2$Var<-factor(data2$Var, levels=c('PercentageMissingFm','PercentageFm'))
data2$hjust<-ifelse(data2$angle < -90, 1, 0)
data2$angle<-ifelse(data2$angle < -90, data2$angle+180, data2$angle)
data2$LabNum<-ifelse(data2$Var =='PercentageFm', data2$Value, '')
data2$NfmInSamp<-ifelse(data2$NfmInSamp>data2$Nfm, data2$Nfm, data2$NfmInSamp)
data2$Lab2<-ifelse(data2$Var =='PercentageFm', paste0(' ', data2$LabNum, '% (',data2$NfmInSamp,'/',data2$Nfm,') '), '')
#data2[data2$Var=='Nsps', 'Value']<-100
#data2<-data2[data2$Var!='NspInSamp',]

# Create dataset
p1<-ggplot(data2) +      
  # Add the stacked bar
  geom_bar(aes(x=Order, y=Value, fill=Var), stat="identity", alpha=0.5) +
  # scale_fill_viridis(discrete=TRUE, direction=-1) +  
  scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(1,0,0,1), "cm")      
  ) +
  coord_polar(clip='off') +
  # Add labels on top of each bar
  geom_text(aes(x=id, y=40, label=Order, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4, angle= data2$angle, inherit.aes = FALSE) +
  geom_text(aes(x=id, y=100, label=Lab2, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4.2, angle= data2$angle, inherit.aes = FALSE )
p1

#PERCENTAGE OF SPECIES
data2 <- as.data.frame(pivot_longer(data, cols=c(Percentage, PercentageMissing), names_to = "Var", values_to = "Value"))

data2$id<-as.numeric(as.factor(data2$Order))
data2$angle= 90 - 360 * (data2$id-0.5) / length(unique(data2$id))     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
data2$Var<-factor(data2$Var, levels=c('PercentageMissing','Percentage'))
data2$hjust<-ifelse(data2$angle < -90, 1, 0)
data2$angle<-ifelse(data2$angle < -90, data2$angle+180, data2$angle)

data2$LabNum<-ifelse(data2$Var =='Percentage', data2$Value, '')
data2$Lab2<-ifelse(data2$Var =='Percentage', paste0(' ',data2$LabNum, '% (',data2$NspInSamp,'/',data2$Nsps,') '), '')

# Create dataset
p2<-ggplot(data2) +      
  # Add the stacked bar
  geom_bar(aes(x=Order, y=Value, fill=Var), stat="identity", alpha=0.5) +
  #scale_fill_viridis(discrete=TRUE, direction=-1) +  
  scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(1,0,0,1), "cm")      
  ) +
  coord_polar(clip='off') +
  # Add labels on top of each bar
  geom_text(aes(x=id, y=40, label=Order, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4, angle= data2$angle, inherit.aes = FALSE) +
  geom_text(aes(x=id, y=100, label=Lab2, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4.2, angle= data2$angle, inherit.aes = FALSE )
p2

#SAVE ALL
pdf('Figures/FamilyCoverageRept.pdf', width=18, height=12)#width=13, height=8)
ggarrange(p2, p1, labels=c('(a)', '(b)'), font.label=list(size=20, face='plain'),
          nrow=2, ncol=2)
dev.off()





# plotBranchbyTrait(b_tree, bfam$mean, 
#                   type = "fan", edge.width=0.5, cex=0.4, show.tip.label=FALSE,
#                   mode='tips', 
#                   palette=palette, legend=TRUE)#viridis




# reptiles ####
r_tree <- ape::read.tree(file = "Data/rep_phyl_Pyron.txt")
r_tree 

plot(r_tree)

# x2 variable used to colour


