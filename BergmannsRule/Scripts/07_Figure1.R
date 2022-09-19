##############################################################
# Authors: 
# Erin Henry, Ana Benitez-Lopez (@anabenlop)
# Email: erinhenry55@gmail.com, abenitez81@gmail.com, ana.benitez@ugr.es
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# https://www.anabenitezlopez.com/

##############################################################
# Description of script and instructions
##############################################################

# Packages and working directory -----------------------------------------------
library(metafor)
library(ggplot2)
library(ggpubr)
library(rphylopic)
library(raster)  # intersect()
library(rworldmap)
library(rworldxtra)

# clean environment
rm(list=ls())

# Figure 1: Map of Occurrences -------------------------------------------------

# get data before divided into cells
# old.dat <- readRDS("Data/BergmannsRule_data_final_20211031.rds")

# get data divided into cells
dat <- readRDS("Data/BergmannsRule_data_forCorrelations_20211114.rds")

# get data used in the analyses 
amphdata <- read.csv("Data/amphdata_ph.csv")
reptdata <- read.csv("Data/reptdata_ph.csv")
birddata <- read.csv("Data/birddata_ph.csv")
mamdata <- read.csv("Data/mamdata_ph.csv")

# merge all data for Figure 1
finaldata <- rbind(amphdata, reptdata, birddata, mamdata)

# subset data to only include species in correlation results
dat <- subset(dat,speciesname %in% finaldata$speciesname)

# make world map
world <- getMap(resolution = "high")
#world <- spTransform(world, CRS("+proj=moll +ellps=WGS84"))
map <- ggplot() + geom_polygon(data = world,
                               aes(x = long, y = lat, group = group),
                               fill = "gray70", 
                               colour = NA) +  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(1,1,1,1),"points"))
map

### Birds

# subset data
birds <- subset(dat,class=="bird")

# get phylopic
img <- image_data("19f3f55c-f942-464b-a61f-27794b9000f7", size = "512")[[1]]

# create map for birds
birdmap <- map + 
  geom_point(data=birds,
             aes(x=x, y=y),
             color="#4477AA",
             size=.15) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  add_phylopic(img,1,-150,-55,ysize=30) # orig 60

# plot
birdmap

### Mammals

# subset data
mammals <- subset(dat,class=="mammal")

# get phylopic
img <- image_data("bb553480-e37f-4236-8c69-ce9fa8116b39", size = "512")[[1]]

# map of mammals
mammap <- map + 
  geom_point(data=mammals,
             aes(x=x, y=y),
             color="#4477AA",
             size=.15) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  add_phylopic(img,1,-150,-50,ysize=45)

# plot
mammap

### Reptiles

# subset data
reptiles <- subset(dat,class=="reptile")

# get phylopic
img <- image_data("0e2a08ed-13a1-4b9e-a047-ef4045e7d88f", size = "512")[[1]]

# reptile map
repmap <- map + 
  geom_point(data=reptiles,
             aes(x=x, y=y),
             color="#4477AA",
             size=.15) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  add_phylopic(img,1,-150,-50,ysize=45)

# plot
repmap

### Amphibians

# subset data
amphibians <- subset(dat,class=="amphibian")

# get phylopic
img <- image_data("4679516b-405b-444f-974d-9775876716e2", size = "512")[[1]]

# amphibian map
amphmap <- map + 
  geom_point(data=amphibians,
             aes(x=x, y=y), #,color=order),
             color="#4477AA",
             size=.15) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  add_phylopic(img,1,-150,-50,ysize=35) # +
#scale_color_manual(name="order",values=c("#4477AA", "#EE6677")) +
#theme(legend.position = "none") # blue = anura, red = caudata
amphmap

### Add all maps to one figure
ggarrange(amphmap,repmap,birdmap,mammap,labels="auto") #,font.label = list(size=8))

# save
ggsave(filename='Figures/Figure_1.png', 
       width=180, height=120, units = 'mm', dpi=600)

# End of script ---
