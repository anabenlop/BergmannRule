##############################################################
# Authors: 
# Erin Henry, Ana Benitez-Lopez (@anabenlop)
# Email: erinhenry55@gmail.com, abenitez81@gmail.com


# Script to create figure 2 in the Bergmann's rule paper
# Made 18 October 2020
# Adapted on 16 December 2021

# clean environment
rm(list = ls())

# Packages and working directory -----------------------------------------------
library(metafor)
library(ggplot2)
library(ggpubr)
library(rphylopic)
library(raster)  # intersect()
# library(rworldmap)
# library(rworldxtra)

# Figure 2: Meta-analysis -----------------------------------------------------

### birds ####

# get data
bird.mods <- readRDS("Results/BergmannsRule_results_MA_birds_phylo_nonphylo.rds")
names(bird.mods)

# create table of values for plot
model <- data.frame(beta = unlist(lapply(bird.mods,"[","beta")),
                    env.var = c("MT","MaxT","NPP","NPPsd"),
                    ci.lb = unlist(lapply(bird.mods,"[","ci.lb")),
                    ci.ub = unlist(lapply(bird.mods,"[","ci.ub")),
                    n = unlist(sapply(bird.mods,"[","k")))
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)
model$hyp <- c("Heat conservation","Heat dissipation","Resource availability",
               "Starvation resistance")

# remove redundant variables
# model <- subset(model, env.var!="MP" & env.var!="MinT" & env.var!="PET")

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


### mammals ####

# get data
mam.mods <- readRDS("Results/BergmannsRule_results_MA_mammals_phylo_nonphylo.rds")
names(mam.mods)

# create table of values for plot
model <- data.frame(beta = unlist(lapply(mam.mods,"[","beta")),
                    env.var = c("MT","MaxT","NPP","NPPsd"),
                    ci.lb = unlist(lapply(mam.mods,"[","ci.lb")),
                    ci.ub = unlist(lapply(mam.mods,"[","ci.ub")),
                    n = unlist(sapply(mam.mods,"[","k")))
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)
model$hyp <- c("Heat conservation","Heat dissipation","Resource availability",
               "Starvation resistance")

# remove redundant variables
# model <- subset(model,env.var!="MP" & env.var!="MinT" & env.var!="PET")

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
  theme_classic() +
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
  xlim(-.2,.2) + 
  #scale_x_continuous(breaks=c(-.1,0, .1)) +
  add_phylopic(img,1,.17,4,ysize=1) 
pm


### reptiles ###

# get data
rep.mods <- readRDS("Results/BergmannsRule_results_MA_reptiles_phylo_nonphylo.rds")
names(rep.mods)

# create table of values for plot
model <- data.frame(beta = unlist(lapply(rep.mods,"[","beta")),
                    env.var = c("NPP","NPPsd"),
                    ci.lb = unlist(lapply(rep.mods,"[","ci.lb")),
                    ci.ub = unlist(lapply(rep.mods,"[","ci.ub")),
                    n = unlist(sapply(rep.mods,"[","k")))
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)
model$hyp <- c("Resource availability","Seasonality")

# keep variables for resource availability and seasonality
# model <- subset(model,env.var=="NPP"|env.var=="NPPsd")


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
amph.mods <-  readRDS("Results/BergmannsRule_results_MA_amphibians_phylo_nonphylo.rds")
names(amph.mods)

# create table of values for plot
model <- data.frame(beta = unlist(lapply(amph.mods,"[","beta")),
                    env.var = c("MP","NPP","NPPsd"),
                    ci.lb = unlist(lapply(amph.mods,"[","ci.lb")),
                    ci.ub = unlist(lapply(amph.mods,"[","ci.ub")),
                    n = unlist(sapply(amph.mods,"[","k")))
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)
model$hyp <- c("Water conservation","Resource availability","Seasonality")

# keep variables for resource availability, seasonality, and water conservation hypotheses
# model <- subset(model,env.var!="MT"&env.var!="MaxT"&env.var!="MinT"&env.var!="PET")

# reorder env.vars
model$env.var <- factor(model$env.var, levels = c("NPPsd","NPP","MeanP"))  

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
  xlim(-.4,.4) +
  add_phylopic(img,1,.35,3.2,ysize=.6)
pa


### put all together ###
ggarrange(pa,pr,pb,pm,ncol=2,nrow=2,
          labels="auto",hjust=-2,vjust=1.53,
          common.legend = F) # save as 1100 x 700

# save
ggsave(filename='Figures/Figure_2.png', 
       width=180, height=120, units = 'mm', dpi=600)

# End of script ####