##############################################################
# Authors: 
# Erin Henry, Ana Benitez-Lopez (@anabenlop)
# Email: erinhenry55@gmail.com, ana.benitez@ugr.es

# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# www.anabenitezlopez.com
# Department of Environmental Science, Radboud University (the Netherlands)
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, Spain)
# Department of Zoology, University of Granada (UGR, Spain)

# First created in 17 May 2022

##############################################################
# Description of script and instructions                  ####
##############################################################

# Script to create figure env variation for the paper

# Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Uncovering the environmental drivers 
# of intraspecific body size variation in terrestrial vertebrates. 

##############################################################

# load libraries
library(metafor)
library(ggplot2)
library(ggpubr)

# clean environment
rm(list= ls())

# Figure 4: Meta-regressions testing the effect of environmental variation ----------------------

### Migration meta-regression panel birds ---------
# get migration models
env.mods <- readRDS("Results/BergmannsRule_results_MR_bird_env.rds")

# Get sample sizes for figure legend
birds <- read.csv("Data/birddata_ph.csv", stringsAsFactors = F)

# dataframe of models
predictors <- data.frame(sd.tavg = seq(from = min(birds$sd.tavg), to = max(birds$sd.tavg), length.out = 1000),
                         sd.tmax = seq(from = min(birds$sd.tmax), to = max(birds$sd.tmax), length.out = 1000),
                         sd.npp = log10(seq(from = min(birds$sd.npp), to = max(birds$sd.npp), length.out = 1000)),
                         sd.npp.sd = log10(seq(from = min(birds$sd.npp.sd), to = max(birds$sd.npp.sd), length.out = 1000)))

#predict for each variable depicting env variation.
names(env.mods)

i <- 2

colpalette <-c("#4477AA","#EE6677","#228833","#AA3377")
tag <-c("a","b","c","d")
df_b <- list() # to store the dataframes

for (i in 1:length(predictors)) {
  
  df_b[[i]] <- predict(env.mods[[i]], newmods = cbind(predictors[,names(env.mods)[i]]), addx = T)
  df_b[[i]] <- data.frame(df_b[[i]])
  df_b[[i]]$pred <- transf.ztor(df_b[[i]]$pred)
  df_b[[i]]$ci.lb <- transf.ztor(df_b[[i]]$ci.lb)
  df_b[[i]]$ci.ub <- transf.ztor(df_b[[i]]$ci.ub)
  df_b[[i]]$var <- rep(names(env.mods)[i], rep = 1000)
  colnames(df_b[[i]])[8] <- "value"
  
}

df_b2 <- rbind(df_b[[1]],df_b[[2]],df_b[[3]],df_b[[4]])
df_b2$var <- factor(df_b2$var,levels = names(env.mods))

B <- ggplot()+ 
            geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+
            theme_bw(base_size=18) +
  #   # annotation_custom(rasterGrob(sil_M,
  #   #                              x = unit(0.14, "npc"),
  #   #                              y = unit(0.15, "npc"),
  #   #                              width = unit(0.22,"npc"),
  #   #                              height = unit(0.27,"npc")),
  #   #                             -Inf, Inf, -Inf, Inf) +
  #   
  #   # geom_point(aes(names(env.mods)[i],z.cor.yi), colour= "#0072B2",size = size[1:1545],shape=20, alpha=I(.3)) +
  #   # scale_shape_identity()+
  geom_line(data=df_b2, aes(x=value, y= pred, color = var), size = 1.2)+
  geom_ribbon(data=df_b2, aes(x=value, ymin=ci.lb,ymax=ci.ub, fill = var),alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black")) +
  xlab("environmental variability (sd)")+
  ylab("Spearman's r")+
  scale_color_manual(values=colpalette) +
  scale_fill_manual(values=colpalette) +
  facet_wrap(~var,nrow=2,ncol=2, scales = "free") + 
  theme(strip.text.x = element_blank(), # remove facet labels
        strip.background = element_blank()) 
B

### Save figure
ggsave(filename='Figures/Figure_4_birds.png', 
       width=260, height=160, units = 'mm', dpi=600)

### Migration meta-regression panel mammals ---------
# get migration models
env.mods <- readRDS("Results/BergmannsRule_results_MR_mam_env.rds")

# Get sample sizes for figure legend
mammals <- read.csv("Data/mamdata_ph.csv", stringsAsFactors = F)

# dataframe of models
predictors <- data.frame(sd.tavg = seq(from = min(mammals$sd.tavg), to = max(mammals$sd.tavg), length.out = 1000),
                         sd.tmax = seq(from = min(mammals$sd.tmax), to = max(mammals$sd.tmax), length.out = 1000),
                         sd.npp = log10(seq(from = min(mammals$sd.npp), to = max(mammals$sd.npp), length.out = 1000)),
                         sd.npp.sd = log10(seq(from = min(mammals$sd.npp.sd), to = max(mammals$sd.npp.sd), length.out = 1000)))

#predict for each variable depicting env variation.
names(env.mods)

colpalette <-c("#4477AA","#EE6677","#228833","#AA3377")
tag <-c("a","b","c","d")
df_b <- list() # to store the dataframes

for (i in 1:length(predictors)) {
  
  df_b[[i]] <- predict(env.mods[[i]], newmods = cbind(predictors[,names(env.mods)[i]]), addx = T)
  df_b[[i]] <- data.frame(df_b[[i]])
  df_b[[i]]$pred <- transf.ztor(df_b[[i]]$pred)
  df_b[[i]]$ci.lb <- transf.ztor(df_b[[i]]$ci.lb)
  df_b[[i]]$ci.ub <- transf.ztor(df_b[[i]]$ci.ub)
  df_b[[i]]$var <- rep(names(env.mods)[i], rep = 1000)
  colnames(df_b[[i]])[8] <- "value"
  
}

df_b2 <- rbind(df_b[[1]],df_b[[2]],df_b[[3]],df_b[[4]])
df_b2$var <- factor(df_b2$var,levels = names(env.mods))

B <- ggplot()+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+
  theme_bw(base_size=18) +
  #   # annotation_custom(rasterGrob(sil_M,
  #   #                              x = unit(0.14, "npc"),
  #   #                              y = unit(0.15, "npc"),
  #   #                              width = unit(0.22,"npc"),
  #   #                              height = unit(0.27,"npc")),
  #   #                             -Inf, Inf, -Inf, Inf) +
  #   
  #   # geom_point(aes(names(env.mods)[i],z.cor.yi), colour= "#0072B2",size = size[1:1545],shape=20, alpha=I(.3)) +
  #   # scale_shape_identity()+
  geom_line(data=df_b2, aes(x=value, y= pred, color = var), size = 1.2)+
  geom_ribbon(data=df_b2, aes(x=value, ymin=ci.lb,ymax=ci.ub, fill = var),alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black")) +
  xlab("environmental variability (sd)")+
  ylab("Spearman's r")+
  scale_color_manual(values=colpalette) +
  scale_fill_manual(values=colpalette) +
  facet_wrap(~var,nrow=2,ncol=2, scales = "free") + 
  theme(strip.text.x = element_blank(), # remove facet labels
        strip.background = element_blank())

B

### Save figure
ggsave(filename='Figures/Figure_4_mammals.png', 
       width=260, height=160, units = 'mm', dpi=600)

### Migration meta-regression panel amphibians ---------
# get migration models
env.mods <- readRDS("Results/BergmannsRule_results_MR_amph_env.rds")

# load data
amphibians <- read.csv("Data/amphdata_ph.csv", stringsAsFactors = F)

#predict for each variable depicting env variation.
names(env.mods)

# dataframe of models
predictors <- data.frame(
                        # sd.tavg = seq(from = min(amphibians$sd.tavg), to = max(amphibians$sd.tavg), length.out = 1000),
                        sd.prec = seq(from = min(amphibians$sd.prec), to = max(amphibians$sd.prec), length.out = 1000),
                         sd.npp = log10(seq(from = min(amphibians$sd.npp), to = max(amphibians$sd.npp), length.out = 1000)),
                         sd.npp.sd = log10(seq(from = min(amphibians$sd.npp.sd), to = max(amphibians$sd.npp.sd), length.out = 1000)))


colpalette <-c("#228833","#AA3377", "#66CCEE")
# tag <-c("a","b","c","d")
df_b <- list() # to store the dataframes

for (i in 1:length(predictors)) {
  
  df_b[[i]] <- predict(env.mods[[i]], newmods = cbind(predictors[,names(env.mods)[i]]), addx = T)
  df_b[[i]] <- data.frame(df_b[[i]])
  df_b[[i]]$pred <- transf.ztor(df_b[[i]]$pred)
  df_b[[i]]$ci.lb <- transf.ztor(df_b[[i]]$ci.lb)
  df_b[[i]]$ci.ub <- transf.ztor(df_b[[i]]$ci.ub)
  df_b[[i]]$var <- rep(names(env.mods)[i], rep = 1000)
  colnames(df_b[[i]])[8] <- "value"
  
}

df_b2 <- rbind(df_b[[1]],df_b[[2]], df_b[[3]])
df_b2$var <- factor(df_b2$var,levels = names(env.mods))

B <- ggplot()+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+
  theme_bw(base_size=18) +
  #   # annotation_custom(rasterGrob(sil_M,
  #   #                              x = unit(0.14, "npc"),
  #   #                              y = unit(0.15, "npc"),
  #   #                              width = unit(0.22,"npc"),
  #   #                              height = unit(0.27,"npc")),
  #   #                             -Inf, Inf, -Inf, Inf) +
  #   
  #   # geom_point(aes(names(env.mods)[i],z.cor.yi), colour= "#0072B2",size = size[1:1545],shape=20, alpha=I(.3)) +
  #   # scale_shape_identity()+
  geom_line(data=df_b2, aes(x=value, y= pred, color = var), size = 1.2)+
  geom_ribbon(data=df_b2, aes(x=value, ymin=ci.lb,ymax=ci.ub, fill = var),alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black")) +
  xlab("environmental variability (sd)")+
  ylab("Spearman's r")+
  scale_color_manual(values=colpalette) +
  scale_fill_manual(values=colpalette) +
  facet_wrap(~var,nrow=2,ncol=2, scales = "free") + 
  theme(strip.text.x = element_blank(), # remove facet labels
        strip.background = element_blank())

B

### Save figure
ggsave(filename='Figures/Figure_4_amphibians.png', 
       width=260, height=160, units = 'mm', dpi=600)

### Migration meta-regression panel reptiles ---------
# get migration models
env.mods <- readRDS("Results/BergmannsRule_results_MR_rept_env.rds")

# load data
reptiles <- read.csv("Data/reptdata_ph.csv", stringsAsFactors = F)

#predict for each variable depicting env variation.
names(env.mods)

# dataframe of models
predictors <- data.frame(
  # sd.tavg = seq(from = min(reptiles$sd.tavg), to = max(reptiles$sd.tavg), length.out = 1000),
  #  sd.tmax = seq(from = min(reptiles$sd.tmax), to = max(reptiles$sd.tmax), length.out = 1000),
  sd.npp = log10(seq(from = min(reptiles$sd.npp), to = max(reptiles$sd.npp), length.out = 1000)),
  sd.npp.sd = log10(seq(from = min(reptiles$sd.npp.sd), to = max(reptiles$sd.npp.sd), length.out = 1000)))


colpalette <-c("#228833","#AA3377")
# tag <-c("a","b","c","d")
df_b <- list() # to store the dataframes

for (i in 1:length(predictors)) {
  
  df_b[[i]] <- predict(env.mods[[i]], newmods = cbind(predictors[,names(env.mods)[i]]), addx = T)
  df_b[[i]] <- data.frame(df_b[[i]])
  df_b[[i]]$pred <- transf.ztor(df_b[[i]]$pred)
  df_b[[i]]$ci.lb <- transf.ztor(df_b[[i]]$ci.lb)
  df_b[[i]]$ci.ub <- transf.ztor(df_b[[i]]$ci.ub)
  df_b[[i]]$var <- rep(names(env.mods)[i], rep = 1000)
  colnames(df_b[[i]])[8] <- "value"
  
}

df_b2 <- rbind(df_b[[1]],df_b[[2]])
df_b2$var <- factor(df_b2$var,levels = names(env.mods))

B <- ggplot()+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+
  theme_bw(base_size=18) +
  #   # annotation_custom(rasterGrob(sil_M,
  #   #                              x = unit(0.14, "npc"),
  #   #                              y = unit(0.15, "npc"),
  #   #                              width = unit(0.22,"npc"),
  #   #                              height = unit(0.27,"npc")),
  #   #                             -Inf, Inf, -Inf, Inf) +
  #   
  #   # geom_point(aes(names(env.mods)[i],z.cor.yi), colour= "#0072B2",size = size[1:1545],shape=20, alpha=I(.3)) +
  #   # scale_shape_identity()+
  geom_line(data=df_b2, aes(x=value, y= pred, color = var), size = 1.2)+
  geom_ribbon(data=df_b2, aes(x=value, ymin=ci.lb,ymax=ci.ub, fill = var),alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black")) +
  xlab("environmental variability (sd)")+
  ylab("Spearman's r")+
  scale_color_manual(values=colpalette) +
  scale_fill_manual(values=colpalette) +
  facet_wrap(~var,nrow=2,ncol=2, scales = "free") + 
  theme(strip.text.x = element_blank(), # remove facet labels
        strip.background = element_blank())

B

### Save figure
ggsave(filename='Figures/Figure_4_reptiles.png', 
       width=260, height=80, units = 'mm', dpi=600)

