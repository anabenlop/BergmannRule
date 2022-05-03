# library
library(tidyverse)
library(viridis)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(scales)
library(dplyr)

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
mamdata2<-mamdata[!duplicated(mamdata$family),] # 49 families

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



#KERNELS

# ET<-read.csv('~/Documents/TraitDatabases/EltonTraits/EltonTraits_Mammals_taxid.csv')
# ET<-ET[!ET$Order %in% c('Sirenia'),]
# ET<-ET[!ET$MSWFamilyLatin %in% c('Octodontidae','Otariidae','Balaenidae','Balaenopteridae','Ziphiidae','Neobalaenidae','Delphinidae','Phocidae','Monodontidae','Dugongidae','Iniidae','Physeteridae','Phocoenidae','Odobenidae','Trichechidae', 'Platanistidae'),]
# ET_BM<-ET$BodyMass.Value
# DT_BM<-read.csv('~/Documents/Progetti/Density/MAMMALS/TRAITS_Mammals.csv', sep=';')$BM_g
# 
# BM<-data.frame(Var=c(rep('All mammals', length(ET_BM)), rep('TetraDENSITY', length(DT_BM))), Value=c(ET_BM, DT_BM))
# 
# p3<-ggplot(BM, aes(x=Value, fill= Var, colour=Var)) +
#   geom_density(alpha = 0.1) + 
#   xlab('log10 Body mass (g)') + ylab('') +
#   scale_x_continuous(trans='log10', labels=trans_format("log10", math_format(10^.x))) +
#   theme_classic() +
#   theme(panel.background = element_rect(fill = "transparent"), 
#         plot.background = element_rect(fill = "transparent", color = NA),
#         plot.margin = unit(c(2.5,1.5,2.5,1.5), "cm"), legend.title = element_blank(), axis.text = element_text(size = 13), axis.title.x = element_text(size = 18), legend.text = element_text(size=14))
# p3
# 
# 
# #plot distribution
# library(rworldmap)
# 
# wmap <- getMap(resolution = "low")
# 
# p4<-ggplot() + 
#   geom_polygon(data=wmap, aes(long, lat, group=group)) +
#   geom_point(data=D, aes(x=X, y=Y), size=1, shape=19, color="#44AA99") +
#   xlab('') + ylab('') +
#   coord_cartesian(ylim=c(-60,100)) +
#   theme_void()
# p4




