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
B <- list()

for (i in 1:length(names(env.mods))) {
  
  df_b <- predict(env.mods[[i]], newmods = cbind(predictors[,i]), addx = T)
  df_b <- data.frame(df_b)
  colnames(df_b)[10] <- names(env.mods)[i]
  df_b$pred <- transf.ztor(df_b$pred)
  df_b$ci.lb <- transf.ztor(df_b$ci.lb)
  df_b$ci.ub <- transf.ztor(df_b$ci.ub)

  #calculate size of points based on sampling variance only
# wi    <- 1/sqrt(birds$z.cor.vi)
# size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

# import silhouette
# raster format
# sil_M <- readPNG("Silhouettes/PhyloPic.72f2f997.Steven-Traver.Cervus-elaphus.png")

B[[i]] <- ggplot()+ 
        geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+
                  theme_bw(base_size=18) +
  # annotation_custom(rasterGrob(sil_M,
  #                              x = unit(0.14, "npc"),
  #                              y = unit(0.15, "npc"),
  #                              width = unit(0.22,"npc"),
  #                              height = unit(0.27,"npc")),
  #                             -Inf, Inf, -Inf, Inf) +
  
  # geom_point(aes(names(env.mods)[i],z.cor.yi), colour= "#0072B2",size = size[1:1545],shape=20, alpha=I(.3)) +
  # scale_shape_identity()+
 
  
  geom_line(data=df_b,aes(df_b[,names(env.mods[i])],pred), colour = colpalette[i], size = 1.2)+
  geom_ribbon(data=df_b, aes(x=df_b[,names(env.mods[i])], ymin=ci.lb,ymax=ci.ub),fill = colpalette[i], alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black")) + 
  xlab(names(env.mods[i]))+
  ylab("Spearman's r")+
  # scale_x_continuous(breaks=seq(round(min(df_b[,names(env.mods[i])]), digits = 1), max(df_b[,names(env.mods[i])]), 3)) +
  # scale_y_continuous(breaks=seq(round(min(df_b[,"pred"]), digits = 2), round(max(df_b[,"pred"]), digits = 2), 0.05)) +
  labs(tag = tag[i])
}

B[[1]]
B[[2]]
B[[3]]
B[[4]]
