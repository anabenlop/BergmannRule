# Script to create figures in the Bergmann's rule paper
# Made 18 October 2020
# Adapted on 17 November 2021

# Packages and working directory -----------------------------------------------
library(metafor)
library(ggplot2)
library(ggpubr)
library(rphylopic)
library(raster)  # intersect()
library(rworldmap)
library(rworldxtra)

setwd("D:/BergmannsRule_upload")

# Figure 1: Map of Occurrences -------------------------------------------------

# get data before divided into cells
old.dat <- readRDS("Data/BergmannsRule_data_final_20211031.rds")

# get data divided into cells
data <- readRDS("Data/BergmannsRule_data_forCorrelations_20211114.rds")

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
birds <- subset(data,class=="bird")

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
mammals <- subset(data,class=="mammal")

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
reptiles <- subset(data,class=="reptile")

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
amphibians <- subset(data,class=="amphibian")

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


# Figure 2: Meta-analysis -----------------------------------------------------

### birds ###

# get data
bird.mods <- readRDS("Results/BergmannsRule_results_MA_birds_20211115.rds")
names(bird.mods)

# create table of values for plot
model <- data.frame(beta = unlist(lapply(bird.mods,"[","beta")),
                    env.var = c("MT","MinT","MaxT","MP","PET","NPP","NPPsd"),
                    ci.lb = unlist(lapply(bird.mods,"[","ci.lb")),
                    ci.ub = unlist(lapply(bird.mods,"[","ci.ub")),
                    n = unlist(sapply(bird.mods,"[","k")))
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)
model$hyp <- c("Heat conservation","Heat conservation","Heat dissipation","NA",
               "Heat dissipation","Resource availability",
               "Starvation resistance")

# remove redundant variables
model <- subset(model, env.var!="MP" & env.var!="MinT" & env.var!="PET")

# reorder env.vars
model$env.var <- factor(model$env.var, levels = c("NPPsd","NPP","MaxT","MT"))  

xmin <- -.2 #round(min(model$ci.lb-.01),digits=2)
xmax <- .2 #round(max(model$ci.ub+.01),digits=2)

# get silhouette
bird <- name_search(text = "thraupidae", options = "namebankID")[[1]]
name_images(uuid = bird$uid[1])
img <- image_data("19f3f55c-f942-464b-a61f-27794b9000f7", size = "512")[[1]] # use first image

# make plot
labels <- seq(round(xmin,digits=1), round(xmax,digits=1), by=.1)
pb <- ggplot() +
  geom_vline(xintercept = 0, color = "gray80",size=.5) +
  geom_point(mapping=aes(x=beta,y=env.var,color=hyp),
             show.legend=T,data=model,size=2) + 
  geom_errorbarh(data=model,show.legend=T,size=.75, 
                 mapping=aes(x=beta,y=env.var,xmin=ci.lb,xmax=ci.ub,color=hyp),
                 height=0) +
  theme_classic2() + labs(title="Birds",x="Spearman's r",y=NULL) +
  theme(axis.title=element_text(size=9),
        plot.title = element_text(size=12),
        axis.text=element_text(size=6),
        plot.margin=unit(c(.5,.5,.5,.5),"lines"),
        panel.border=element_rect(fill=NA),
        legend.title=element_text(size=9),
        legend.text=element_text(size=6),
        legend.spacing.y = unit(.1,'cm'),
        legend.spacing.x = unit(.1,'cm'),
        legend.box.margin=margin(-10,-10,-10,-10)) + 
  scale_color_manual(name="Hypothesis",
                     values=c("#4477AA", "#EE6677","#228833","#AA3377")) + 
  xlim(-.2,.2) + 
  add_phylopic(img,1,.17,4.15,ysize=.7) # orig size = 1
pb



### mammals ###

# get data
mam.mods <- readRDS("Results/BergmannsRule_results_MA_mammals_20211115.rds")
names(mam.mods)

# create table of values for plot
model <- data.frame(beta = unlist(lapply(mam.mods,"[","beta")),
                    env.var = c("MT","MinT","MaxT","MP","PET","NPP","NPPsd"),
                    ci.lb = unlist(lapply(mam.mods,"[","ci.lb")),
                    ci.ub = unlist(lapply(mam.mods,"[","ci.ub")),
                    n = unlist(sapply(mam.mods,"[","k")))
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)
model$hyp <- c("Heat conservation","Heat conservation","Heat dissipation",
               "NA","Heat dissipation","Resource availability",
               "Starvation resistance")

# remove redundant variables
model <- subset(model,env.var!="MP" & env.var!="MinT" & env.var!="PET")

# reorder env.vars
model$env.var <- factor(model$env.var, levels = c("NPPsd","NPP","MaxT","MT"))  

xmin = -.2 #round(min(model$ci.lb-.01),digits=2)
xmax = .2 #round(max(model$ci.ub+.01),digits=2)

# get silhouette
mam <- name_search(text = "pliocervus turolensis", options = "namebankID")[[1]]
name_images(uuid = mam$uid[1])
img <- image_data("bb553480-e37f-4236-8c69-ce9fa8116b39", size = "512")[[1]] # use first image

# make plot
labels <- seq(round(xmin,digits=1), round(xmax,digits=1), by=.1)
pm <- ggplot() +
  geom_vline(xintercept = 0, color = "gray80",size=.5) +
  geom_point(mapping=aes(x=beta,y=env.var,color=hyp),
             show.legend=T,data=model,size=2) + 
  geom_errorbarh(data=model,show.legend=T,size=.75, 
                 mapping=aes(x=beta,y=env.var,xmin=ci.lb,xmax=ci.ub,color=hyp),
                 height=0) +
  theme_classic2() +
  labs(title="Mammals",x="Spearman's r",y=NULL) +
  theme(axis.title=element_text(size=9),
        plot.title = element_text(size=12),
        axis.text=element_text(size=6),
        plot.margin=unit(c(.5,.5,.5,.5),"lines"),
        panel.border=element_rect(fill=NA),
        legend.title=element_text(size=9),
        legend.text=element_text(size=6),
        legend.spacing.y = unit(.1,'cm'),
        legend.spacing.x = unit(.1,'cm'),
        legend.box.margin=margin(-10,-10,-10,-10)) + 
  scale_color_manual(name="Hypothesis",
                     values=c("#4477AA", "#EE6677","#228833","#AA3377")) + 
  #xlim(-.2,.2) + 
  #scale_x_continuous(breaks=c(-.1,0, .1)) +
  add_phylopic(img,1,.65,4,ysize=1) 
pm


### reptiles ###

# get data
rep.mods <- readRDS("Results/BergmannsRule_results_MA_reptiles_20211115.rds")
names(rep.mods)

# create table of values for plot
model <- data.frame(beta = unlist(lapply(rep.mods,"[","beta")),
                    env.var = c("MT","MinT","MaxT","MP","PET","NPP","NPPsd"),
                    ci.lb = unlist(lapply(rep.mods,"[","ci.lb")),
                    ci.ub = unlist(lapply(rep.mods,"[","ci.ub")),
                    n = unlist(sapply(rep.mods,"[","k")))
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)
model$hyp <- c("NA","NA","NA","NA","NA","Resource availability","Seasonality")

# keep variables for resource availability and seasonality
model <- subset(model,env.var=="NPP"|env.var=="NPPsd")


# reorder env.vars
model$env.var <- factor(model$env.var, levels = c("NPPsd","NPP"))  

xmin = -.2 #round(min(model$ci.lb-.01),digits=2)
xmax = .2 #round(max(model$ci.ub+.01),digits=2)

# get silhouette
rep <- name_search(text = "egernia saxatilis", options = "namebankID")[[1]]
name_images(uuid = rep$uid[1])
img <- image_data("0e2a08ed-13a1-4b9e-a047-ef4045e7d88f", size = "512")[[1]] # use first image

# make plot
labels <- seq(round(xmin,digits=1), round(xmax,digits=1), by=.1)
pr <- ggplot() +
  geom_vline(xintercept = 0, color = "gray80",size=.5) +
  geom_point(mapping=aes(x=beta,y=env.var,color=hyp),
             show.legend=T,data=model,size=2) + 
  geom_errorbarh(data=model,show.legend=T,size=.75, 
                 mapping=aes(x=beta,y=env.var,xmin=ci.lb,xmax=ci.ub,color=hyp),
                 height=0) +
  theme_classic2() +
  labs(title="Reptiles",x="Spearman's r",y=NULL) +
  theme(axis.title=element_text(size=9),
        plot.title = element_text(size=12),
        axis.text=element_text(size=6),
        plot.margin=unit(c(.5,.5,.5,.5),"lines"),
        panel.border=element_rect(fill=NA),
        legend.title=element_text(size=9),
        legend.text=element_text(size=6),
        legend.spacing.y = unit(.1,'cm'),
        legend.spacing.x = unit(.1,'cm'),
        legend.box.margin=margin(-10,-10,-10,-10)) + 
  scale_color_manual(name="Hypothesis",values=c("#228833","#AA3377")) + 
  xlim(-.2,.2) +
  add_phylopic(img,1,.166,2.25,ysize=.5)
pr



### amphibians ###

# get data
amph.mods <-  readRDS("Results/BergmannsRule_results_MA_amphibians_20211115.rds")
names(amph.mods)

# create table of values for plot
model <- data.frame(beta = unlist(lapply(amph.mods,"[","beta")),
                    env.var = c("MT","MinT","MaxT","MP","PET","NPP","NPPsd"),
                    ci.lb = unlist(lapply(amph.mods,"[","ci.lb")),
                    ci.ub = unlist(lapply(amph.mods,"[","ci.ub")),
                    n = unlist(sapply(amph.mods,"[","k")))
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)
model$hyp <- c("NA","NA","NA","Water conservation","Water conservation","Resource availability","Seasonality")

# keep variables for resource availability, seasonality, and water conservation hypotheses
model <- subset(model,env.var!="MT"&env.var!="MaxT"&env.var!="MinT"&env.var!="PET")

# reorder env.vars
model$env.var <- factor(model$env.var, levels = c("NPPsd","NPP","MP"))  

xmin = -.2 #round(min(model$ci.lb-.01),digits=2)
xmax = .2 #round(max(model$ci.ub+.01),digits=2)

# get silhouette
amph <- name_search(text = "dendrobates azureus", options = "namebankID")[[1]]
name_images(uuid = amph$uid[1])
img <- image_data("4679516b-405b-444f-974d-9775876716e2", size = "512")[[1]] # use first image

# make plot
labels <- seq(round(xmin,digits=1), round(xmax,digits=1), by=.1)
pa <- ggplot() +
  geom_vline(xintercept = 0, color = "gray80",size=.5) +
  geom_point(mapping=aes(x=beta,y=env.var,color=hyp),show.legend=T,data=model,size=2) + 
  geom_errorbarh(data=model,show.legend=T,size=.75, 
                 mapping=aes(x=beta,y=env.var,xmin=ci.lb,xmax=ci.ub,color=hyp),
                 height=0) +
  theme_classic2() +
  labs(title="Amphibians",x="Spearman's r",y=NULL) +
  theme(axis.title=element_text(size=9),
        plot.title = element_text(size=12),
        axis.text=element_text(size=6),
        plot.margin=unit(c(.5,.5,.5,.5),"lines"),
        panel.border=element_rect(fill=NA),
        legend.title=element_text(size=9),
        legend.text=element_text(size=6),
        legend.spacing.y = unit(.1,'cm'),
        legend.spacing.x = unit(.1,'cm'),
        legend.box.margin=margin(-10,-10,-10,-10)) + 
  scale_color_manual(name="Hypothesis",
                     values=c("#228833","#AA3377","#66CCEE")) + 
  xlim(-.3,.3) +
  add_phylopic(img,1,.27,3.2,ysize=.6)
pa


### put all together ###
ggarrange(pa,pr,pb,pm,ncol=2,nrow=2,
          labels="auto",hjust=-2,vjust=1.53,
          common.legend = F) # save as 1100 x 700

# save
ggsave(filename='Figures/Figure_2.png', 
       width=180, height=120, units = 'mm', dpi=600)


# Figure 3: Meta-regressions (Heat Balance and Migration) ----------------------

### Migration meta-regression panel
# get migration models
mig.mods <- readRDS("Results/Bergmannsrule_results_MR_mig.rds")

# Get sample sizes for figure legend
birds <- readRDS("Results/BergmannsRule_results_correlationsBirdsMR_20211115.rds")
nrow(subset(birds,migration=="Resident"))/7 # N resident species
nrow(subset(birds,migration=="Migratory"))/7 # N migratory species

names(mig.mods)

# dataframe of models
model <- data.frame(beta = unlist(lapply(mig.mods,"[","beta")),
                    env.var = rep(c("MT","MinT","MaxT","MP","PET","NPP","NPPsd"),each=2),
                    ci.lb = unlist(lapply(mig.mods,"[","ci.lb")),
                    ci.ub = unlist(lapply(mig.mods,"[","ci.ub")),
                    n = unlist(sapply(mig.mods,"[","k")),
                    migration = rep(c("Migratory","Resident"),times=7),
                    ID = 1:14)
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)

model <- subset(model,env.var!="MP" & env.var!="MinT" & env.var!="PET")
model$env.var <- factor(model$env.var,levels=c("MT","MaxT","NPP","NPPsd"))

# make plot
p <- ggplot(model,aes(x=beta,y=migration),show.legend=T) +
  geom_vline(xintercept = 0, color = "gray80",size=.5) +
  geom_point(mapping=aes(x=beta,y=migration,color=env.var),show.legend=T,
             data=model,size=2) + 
  geom_errorbarh(data=model,show.legend=F,size=.75, 
                 mapping=aes(x=beta,y=migration,color=env.var,xmin=ci.lb,
                             xmax=ci.ub,height=0)) +
  theme_classic2() +
  labs(x="Spearman's r",y=NULL) +
  scale_color_manual(name="Variable",
                     values=c("#4477AA","#EE6677","#228833","#AA3377")) +
  theme(axis.title=element_text(size=9),
        plot.title = element_text(size=12),
        axis.text=element_text(size=8),
        plot.margin=unit(c(1,5,1,1),"lines"),
        panel.border=element_rect(fill=NA),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9))

p <- p + facet_wrap(~env.var,nrow=4) + 
  theme(strip.text.x = element_blank(), # remove facet labels
        strip.background = element_blank())

p



### Heat balance meta-regression (Panel b)
hb.mod <- readRDS("Results/BergmannsRule_results_MR_heatBalance.rds")

# Get sample sizes for figure legend
hb <- readRDS("Results/BergmannsRule_results_correlations_20211114.rds")
nrow(subset(hb,class=="reptile" |order=="Anura"))/7 # N thermoregulators
nrow(subset(hb,order=="Caudata"))/7 # N thermoconformers
rm(hb)

# dataframe of models
model <- data.frame(beta = hb.mod$beta,
                    env.var <- c("MT","MT"),
                    ci.lb = hb.mod$ci.lb,
                    ci.ub = hb.mod$ci.ub,
                    therm = c("TC","TR"))
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)

model <- subset(model,env.var!="MP" & env.var!="MinT" & env.var!="PET")
model$therm <- factor(model$therm,levels=c("TR","TC"))

# make plot
p2 <- ggplot(model,aes(x=beta,y=therm),show.legend=F) +
  geom_vline(xintercept = 0, color = "gray80",size=.55) +
  geom_point(mapping=aes(x=beta,y=therm,color=env.var),show.legend=F,
             data=model,size=2) + 
  geom_errorbarh(data=model,show.legend=F,size=.75, 
                 mapping=aes(x=beta,y=therm,color=env.var,xmin=ci.lb,
                             xmax=ci.ub,height=0)) +
  theme_classic2() +
  labs(x="Spearman's r",y=NULL) +
  scale_color_manual(name="Variable",
                     values=c("#4477AA")) +
  theme(axis.title=element_text(size=9),
        plot.title = element_text(size=12),
        axis.text=element_text(size=8),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9),
        plot.margin=unit(c(1,5,1,1),"lines"),
        panel.border=element_rect(fill=NA),
        legend.spacing.x = unit(.1,'cm'),
        )
p2

### Combine (Into Fig 3.)
ggarrange(p2,p,ncol=2,nrow=1,
          labels="auto",#hjust=-5,vjust=2,
          common.legend = T,
          legend = "bottom") # save as 1000 x 400

### Save figure
ggsave(filename='Figures/Figure_3.png', 
       width=180, height=80, units = 'mm', dpi=600)


# Table S3.2: Meta-analysis results --------------------------------------------

# get meta analysis data
b.mods <- readRDS("Results/BergmannsRule_results_MA_birds_20211115.rds")
m.mods <- readRDS("Results/BergmannsRule_results_MA_mammals_20211115.rds")
r.mods <- readRDS("Results/BergmannsRule_results_MA_reptiles_20211115.rds")
a.mods <- readRDS("Results/BergmannsRule_results_MA_amphibians_20211115.rds")

# view results
x <- a.mods$prec # fill in model

round(transf.ztor(x$beta), digits=3) # estimate
round(transf.ztor(x$ci.lb), digits=3) # ci lower
round(transf.ztor(x$ci.ub), digits=3) # ci upper
x$pval # p value
round(x$QE, digits=1) # QE
x$QEp # p val for QE

### percentages
dat <- readRDS('Results/BergmannsRule_results_correlations_20211114.rds')
am <- subset(dat,class=='amphibian')
re <- subset(dat,class=='reptile')
bi <- subset(dat,class=='bird')
ma <- subset(dat,class=='mammal')

x <- subset(am,env.var=='prec') # set data

round(nrow(subset(x, z.cor.yi > 0)) / nrow(x) * 100,digits=2) #% pos
round(nrow(subset(x, z.cor.yi < 0)) / nrow(x) * 100,digits=2) #% neg


# Table S3.3: Migration meta-regression ---------------------------------------

b.mods <- readRDS("Results/BergmannsRule_results_MR_mig.rds")

# view results
x <- b.mods$tavg # fill in model

round(transf.ztor(x$beta), digits=3) # estimate
round(transf.ztor(x$ci.lb), digits=3) # ci lower
round(transf.ztor(x$ci.ub), digits=3) # ci upper
x$pval # p value
round(x$QM,digits=2) # QM
x$QMp # p val for QM
round(x$QE, digits=1) # QE
x$QEp # p val for QE

# Appendix S1: Data sources ----------------------------------------------------

# get data used for correlations
cor.dat <- readRDS('Results/BergmannsRule_results_correlations_20211114.rds')

# get data before pooling into cells
dat <- readRDS('Data/BergmannsRule_data_final_20211031.rds')

# subset data to only include species in correlation results
dat <- subset(dat,speciesname %in% cor.dat$speciesname)

EH <- subset(dat,database=='Henry Lit Search')
x <- plyr::count(dat,vars='citation')


# Results section --------------------------------------------------------------

# get data used for correlations
cor.dat <- readRDS('Results/BergmannsRule_results_correlations_20211114.rds')
cor.dat <- subset(cor.dat,env.var=='tavg')

# species numbers
nrow(cor.dat)
nrow(subset(cor.dat,class=='amphibian'))
nrow(subset(cor.dat,class=='reptile'))
nrow(subset(cor.dat,class=='bird'))
nrow(subset(cor.dat,class=='mammal'))


# Other data for paper ---------------------------------------------------------
# Check latitudinal range of classes
dat <- readRDS('Data/BergmannsRule_data_forCorrelations_20211102.rds')
am <- subset(dat, class== 'amphibian')
re <- subset(dat, class== 'reptile')
bi <- subset(dat, class== 'bird')
ma <- subset(dat, class== 'mammal')

summary(am)



range(am$y) # ~48
range(re$y) # ~56
range(bi$y) # ~125
range(ma$y) # ~133

# The latitudinal range of amphibians and reptiles is narrower than that of
# mammals and birds. Mammals and birds cover ~2x the latitudinal range. This
# may explain why environmental factors did not explain much of the body size 
# variation for ectotherms. It is likely that body size variation in ectotherms
# is instead governed by factors at a smaller scale.

boxplot(am$y, re$y, bi$y, ma$y)

# Check sample size
results <- readRDS('Results/BergmannsRule_results_correlations_20211102.rds')
nrow(subset(results, env.var == 'tavg')) # total species
nrow(subset(results, env.var == 'tavg' & class == 'amphibian'))
nrow(subset(results, env.var == 'tavg' & class == 'reptile'))
nrow(subset(results, env.var == 'tavg' & class == 'bird'))
nrow(subset(results, env.var == 'tavg' & class == 'mammal'))

sum(subset(results, env.var == 'tavg')$freq) # total cells
sum(subset(results, env.var == 'tavg' & class == 'amphibian')$freq)
sum(subset(results, env.var == 'tavg' & class == 'reptile')$freq)
sum(subset(results, env.var == 'tavg' & class == 'bird')$freq)
sum(subset(results, env.var == 'tavg' & class == 'mammal')$freq)

# Full dataset only including species in analysis
dat <- readRDS('Data/BergmannsRule_data_final_20211031.rds')
dat <- subset(dat, speciesname %in% results$speciesname) # include species in results

nrow(dat)
nrow(subset(dat, class == 'amphibian'))
nrow(subset(dat, class == 'reptile'))
nrow(subset(dat, class == 'bird'))
nrow(subset(dat, class == 'mammal'))

