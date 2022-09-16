################################################################################

# Author(s): Erin Henry (erin.henry@wur.nl)

# Created on 15 September 2022

# Description:

# Creates supplementary figures: S4.4 - S4.7 for the paper:

# Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Unveiling the 
# environmental drivers of intraspecific body size variation in terrestrial 
# vertebrates. 

# Creates histograms of correlation coefficients, including one figure per 
# class, each with one histogram per environmental variable tested. 
# The weighting of the correlation coefficients is represented by the inverse of 
# the sampling variance (1/vi). These values were split into quartiles:
# "0-25%", "25-50%", "50%-75%", and "75-100%", and illustrated by different
# shades.


# Packages and working directory -----------------------------------------------
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(metafor)
library(RColorBrewer)
library(rphylopic)
library(dplyr)


# Function to split sampling variance (vi) into quantiles ----------------------

q <- function(cor.dat) {
  
  #find quantile of weights (inverse of sampling variance, vi)
  cor.dat$inv.vi <- cor.dat$z.cor.vi ^ -1 # calculate inverse of sampling variance
  quant.vi <- quantile(cor.dat$inv.vi)

  # Code sampling variance into levels
  cor.dat <- cor.dat %>%
    mutate(
      w.quant = case_when(
        inv.vi <= quant.vi[[2]] ~ "0-25%",  
        inv.vi > quant.vi[[2]] & inv.vi <= quant.vi[[3]] ~ "25-50%",
        inv.vi > quant.vi[[3]] & inv.vi <= quant.vi[[4]] ~ "50-75%",
        inv.vi > quant.vi[[4]] & inv.vi <= quant.vi[[5]] ~ "75-100%",
      )
  )
  cor.dat$w.quant <- factor(cor.dat$w.quant, levels = c("0-25%","25-50%","50-75%","75-100%"))
  return(cor.dat)
}


# Define theme and other arguments to be used in each figure -------------------

# standard theme
my.theme <- theme(axis.title=element_text(size=14),
                  plot.title = element_text(size=18),
                  legend.title = element_text(size=12),
                  legend.text=element_text(size=10),
                  axis.text=element_text(size=12))

# misc arguments
my.args <- list(geom_histogram(breaks=seq(-1, 1, by=.2),colour = 'white',alpha=.9),
                xlim(-1,1),
                theme_classic(),
                labs(x="Spearman's r",y="Frequency"),
                geom_vline(xintercept=0,color="gray40",size = .5),
                scale_y_continuous(expand = expansion(mult = c(0, 0.05))))

# Color sets -------------------------------------------------------------------
# Environmental variables are represented by the same colors in each figure.
# Here, differences in weighting are represented by shades of each color.

# colors tavg
c.tavg <- scale_fill_manual(name="Quartile 1/vi",
                            values = c('#91b2d3', '#6794c2','#4477aa','#335a80'))

# colors tmax
c.tmax <- scale_fill_manual(name="Quartile 1/vi",
                            values = c('#fad0d5', '#f49ba6','#EE6677','#e83148'))

# colors npp
c.npp <- scale_fill_manual(name="Quartile 1/vi",
                            values = c('#4dd363', '#2eb745','#228833','#165921'))

# colors npp.sd
c.npp.sd <- scale_fill_manual(name="Quartile 1/vi",
                            values = c('#d77bb0','#ca4e95','#AA3377','#7d2557'))

# colors prec
c.prec <- scale_fill_manual(name="Quartile 1/vi",
                            values = c('#d0effa', '#9bdef4', '#66CCEE', '#31bae8'))



# Fig. S4.4: Amphibians --------------------------------------------------------
# get data
am <- read.csv("amphdata_ph.csv")

# split into dataframes based on env.var
am.list <- split(am, f = am$env.var)
am.list <- lapply(am.list, FUN=q)

# number of species for figure legend
nrow(am.list$npp)

# get silhouette
a.img <- image_data("4679516b-405b-444f-974d-9775876716e2", size = "512")[[1]]

# plots
aprec <- ggplot(data = am.list$prec, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.prec + theme(legend.position="right")

anpp <- ggplot(data = am.list$npp, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.npp + theme(legend.position="right") +
  add_phylopic(a.img,1,x=.85,y=6.5,ysize=1.5) # include phylopic

anpp.sd <- ggplot(data = am.list$npp.sd, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.npp.sd + theme(legend.position="right")

# get legend from original histogram
ap <- ggplot(data = am.list$npp, aes(x = corr.coeff,fill=env.var)) + 
  my.theme + my.args + theme(legend.position="bottom") +
  scale_fill_manual(name=element_blank(),
                    breaks=c('prec', 'npp','npp.sd'),
                    values=c('prec'='#66CCEE','npp'='#228833', 'npp.sd' = '#AA3377'),
                    labels=c('MeanP', 'NPP', 'NPPsd'))
a.leg <- get_legend(ap)
rm(ap)

# create one plot with 3 panels and common legend
aplot <- grid.arrange(arrangeGrob(aprec,anpp,anpp.sd,nrow=2), 
             a.leg, 
             nrow=2,heights=c(10, 1))
# save
#ggsave(plot=aplot,filename='BergmannsRule_Figure_S4.4.png', 
#       width=180, height=120, units = 'mm', dpi=600)

# remove objects
rm(a.img,am,am.list,aprec,anpp,anpp.sd)

###


# Fig. S4.5: Reptiles ----------------------------------------------------------
# get data
re <- read.csv("reptdata_ph.csv")

# split into dataframes based on env.var
re.list <- split(re, f = re$env.var)
re.list <- lapply(re.list, FUN=q)

# number of species for figure legend
nrow(re.list$npp)

# get silhouette
r.img <- image_data("0e2a08ed-13a1-4b9e-a047-ef4045e7d88f", size = "512")[[1]]

# plots
rnpp <- ggplot(data = re.list$npp, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.npp

rnpp.sd <- ggplot(data = re.list$npp.sd, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.npp.sd +
  add_phylopic(r.img,1,x=.85,y=15,ysize=3) # include phylopic

# get legend from original histogram
rp <- ggplot(data = re.list$npp, aes(x = corr.coeff,fill=env.var)) + 
  my.theme + my.args + theme(legend.position="bottom") +
  scale_fill_manual(name=element_blank(),
                    breaks=c('npp','npp.sd'),
                    values=c('npp'='#228833', 'npp.sd' = '#AA3377'),
                    labels=c('NPP', 'NPPsd'))
r.leg <- get_legend(rp)
rm(rp)

# create one plot with 2 panels and common legend
rplot <- grid.arrange(arrangeGrob(rnpp,rnpp.sd,ncol=2), 
                      r.leg, 
                      nrow=2,heights=c(10, 1))
# save
#ggsave(plot=rplot,filename='BergmannsRule_Figure_S4.5.png', 
#       width=180, height=80, units = 'mm', dpi=600)

# remove objects
rm(r.img,re,re.list,rnpp,rnpp.sd)

###


# Fig. S4.6: Birds -------------------------------------------------------------
# get data
bi <- read.csv("birddata_ph.csv")

# split into dataframes based on env.var
bi.list <- split(bi, f = bi$env.var)
bi.list <- lapply(bi.list, FUN=q)

# number of species for figure legend
nrow(bi.list$npp)

# get silhouette
b.img <- image_data("19f3f55c-f942-464b-a61f-27794b9000f7", size = "512")[[1]]

# plots
btavg <- ggplot(data = bi.list$tavg, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.tavg

btmax <- ggplot(data = bi.list$tmax, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.tmax +
  add_phylopic(b.img,1,x=.8,y=300,ysize=60) # orig size = 1

bnpp <- ggplot(data = bi.list$npp, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.npp 

bnpp.sd <- ggplot(data = bi.list$npp.sd, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.npp.sd

# get legend from original histogram
bp <- ggplot(data = bi.list$npp, aes(x = corr.coeff,fill=env.var)) + 
  my.theme + my.args + theme(legend.position="bottom") +
  scale_fill_manual(name=element_blank(),
                    breaks=c('tavg', 'tmax', 'npp','npp.sd'),
                    values=c('tavg'='#4477AA', 'tmax'='#EE6677', 
                             'npp'='#228833', 'npp.sd' = '#AA3377'),
                    labels=c('MeanT','MaxT', 'NPP', 'NPPsd'))
b.leg <- get_legend(bp)
rm(bp)

# create one plot with 4 panels and common legend
bplot <- grid.arrange(arrangeGrob(btavg,btmax,bnpp,bnpp.sd,ncol=2), 
                      b.leg, 
                      nrow=2,heights=c(10, 1))
# save
#ggsave(plot=bplot,filename='BergmannsRule_Figure_S4.6.png', 
#       width=180, height=120, units = 'mm', dpi=600)

# remove objects
rm(b.img,bi,bi.list,btavg,btmax,bnpp,bnpp.sd)

###


# Fig. S4.7: Mammals -----------------------------------------------------------
# get data
ma <- read.csv("mamdata_ph.csv")

# split into dataframes based on env.var
ma.list <- split(ma, f = ma$env.var)
ma.list <- lapply(ma.list, FUN=q)

# number of species for figure legend
nrow(ma.list$npp)

# get silhouette
m.img <- image_data("bb553480-e37f-4236-8c69-ce9fa8116b39", size = "512")[[1]]

# plots
mtavg <- ggplot(data = ma.list$tavg, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.tavg

mtmax <- ggplot(data = ma.list$tmax, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.tmax +
  add_phylopic(m.img,1,x=.8,y=120,ysize=30) # include phylopic

mnpp <- ggplot(data = ma.list$npp, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.npp 

mnpp.sd <- ggplot(data = ma.list$npp.sd, aes(x = corr.coeff,fill=w.quant)) + 
  my.theme + my.args + c.npp.sd

# get legend from original histogram
mp <- ggplot(data = ma.list$npp, aes(x = corr.coeff,fill=env.var)) + 
  my.theme + my.args + theme(legend.position="bottom") +
  scale_fill_manual(name=element_blank(),
                    breaks=c('tavg', 'tmax', 'npp','npp.sd'),
                    values=c('tavg'='#4477AA', 'tmax'='#EE6677', 
                             'npp'='#228833', 'npp.sd' = '#AA3377'),
                    labels=c('MeanT','MaxT', 'NPP', 'NPPsd'))
m.leg <- get_legend(mp)
rm(mp)

# create one plot with 4 panels and common legend
mplot <- grid.arrange(arrangeGrob(mtavg,mtmax,mnpp,mnpp.sd,ncol=2), 
                      m.leg, 
                      nrow=2,heights=c(10, 1))
# save
#ggsave(plot=mplot,filename='BergmannsRule_Figure_S4.7.png', 
#       width=180, height=120, units = 'mm', dpi=600)

# remove objects
rm(m.img,ma,ma.list,mtavg,mtmax,mnpp,mnpp.sd)

###
