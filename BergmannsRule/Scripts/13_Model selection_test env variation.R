##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 12th of December 2021

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script fits the phylogenetic meta-analysis and metaregressions with environmental variability 
# and performs model selection to test whetherspecies exposed to a wider ranfe of variation in environmental
# conditions are more likely to adhere to the tested hypothesis for mammals for the paper: 


# Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Uncovering the environmental drivers 
# of intraspecific body size variation in terrestrial vertebrates. 


##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
# library(ggplot2)
# library(ggpubr)
library(tictoc)
# library(png)
# library(ggimage)
# library(grid)
# library(rsvg)
# library(grImport2)

#clean memory
rm(list=ls())

##############################################################
# Importing datasets                                      ####
##############################################################
# 1. Amphibians ---------------------------------------------------------------------

#Load data
amphibians_ph <- read.csv("Data/amphdata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/Phylogeny/amph_phylo_cor.Rdata") #amph_phylo_cor

# vector of environmental variables
env.vars <- c('prec','npp','npp.sd')

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = amph_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# for loop running a meta-analysis for each environmental variable
# models are run with "ML" instead of "REML" so that model selection based on BIC, LRT or similar is meaningful
tic("Run phylo meta-analysis in a loop")
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('amph.',env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(amphibians_ph, env.var == env.vars[i]),
                random = RE, R = phylocor, method = "ML"))
}

toc()

# save results
amph.ma <- list(amph.prec, amph.npp, amph.npp.sd)
names(amph.ma) <- c('prec','npp','npp.sd')

saveRDS(amph.ma,
        'Results/BergmannsRule_results_MA_amphibians_phylo_nonphylo_ML.rds')

# env variation models
# npp model with env variation
amph.npp.env <- rma.mv(yi = z.cor.yi,
                       V = z.cor.vi,
                       data = amphibians_ph,
                       subset = env.var=="npp",
                       mods = ~ log10(sd.npp),
                       random = RE, R = phylocor, method = "ML")
summary(amph.npp.env) # negative effects, non significant

# save results
# saveRDS(amph.npp.env,"Results/BergmannsRule_results_MR_amph_npp_env.rds")
# rm(amph.npp.env)

# npp.sd model with env variation
amph.npp.sd.env <- rma.mv(yi = z.cor.yi,
                          V = z.cor.vi,
                          data = amphibians_ph,
                          subset = env.var=="npp.sd",
                          mods = ~ log10(sd.npp.sd),
                          random = RE, R = phylocor,method = "ML")
summary(amph.npp.sd.env) # no clear support for the seasonality hypothesis

# save results
# saveRDS(amph.npp.sd.env,"Results/BergmannsRule_results_MR_amph_nppsd_env.rds")
# rm(amph.npp.sd.env)

# prec model with env variation
amph.prec.env <- rma.mv(yi = z.cor.yi,
                        V = z.cor.vi,
                        data = amphibians_ph,
                        subset = env.var=="prec",
                        mods = ~ sd.prec,
                        random = RE, R = phylocor, method = "ML")
summary(amph.prec.env) # tendency towards smaller size as species are exposed to more variation in precipitation

# save results
# saveRDS(amph.prec.env,"Results/BergmannsRule_results_MR_amph_prec_env.rds")
# rm(amph.prec.env)

# # load fitted models
# amph.npp.env <- readRDS('Results/BergmannsRule_results_MR_amph_npp_env.rds')
# amph.npp.sd.env <- readRDS('Results/BergmannsRule_results_MR_amph_nppsd_env.rds')
# amph.prec.env <- readRDS('Results/BergmannsRule_results_MR_amph_prec_env.rds')

# save results in list
am.mr.env <- list(amph.npp.env,amph.npp.sd.env,amph.prec.env) 
names(am.mr.env) <- c("sd.npp","sd.npp.sd","sd.prec")

saveRDS(am.mr.env,'Results/BergmannsRule_results_MR_amph_env_ML.rds')

anova(amph.ma$prec, am.mr.env$sd.prec) # model with env var slightly better
anova(amph.ma$npp, am.mr.env$sd.npp) # null model better
anova(amph.ma$npp.sd, am.mr.env$sd.npp.sd) # null model better

# 2. Reptiles ---------------------------------------------------------------------

#Load data
reptiles_ph <- read.csv("Data/reptdata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/rept_phylo_cor.Rdata") #rept_phylo_cor

# vector of environmental variables
env.vars <- c('npp','npp.sd')

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = rept_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# for loop running a meta-analysis for each environmental variable
tic("Run phylo meta-analysis in a loop")
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('rept.',env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(reptiles_ph, env.var == env.vars[i]),
                random = RE, R = phylocor, method ="ML"))
}

toc()

# save results
rept.ma <- list(rept.npp, rept.npp.sd)
names(rept.ma) <- c('npp','npp.sd')

saveRDS(rept.ma,
        'Results/BergmannsRule_results_MA_reptiles_phylo_nonphylo_ML.rds')

# env variation models
# npp model with env variation
rept.npp.env <- rma.mv(yi = z.cor.yi,
                       V = z.cor.vi,
                       data = reptdata,
                       subset = env.var=="npp",
                       mods = ~ log10(sd.npp),
                       random = RE, R = phylocor)
summary(rept.npp.env) # tendency to positive effect (large CI) with variation in npp across the gradient of body size records

# save results
saveRDS(rept.npp.env,"Results/BergmannsRule_results_MR_rept_npp_env.rds")
rm(rept.npp.env)

# npp.sd model with env variation
rept.npp.sd.env <- rma.mv(yi = z.cor.yi,
                          V = z.cor.vi,
                          data = reptdata,
                          subset = env.var=="npp.sd",
                          mods = ~ log10(sd.npp.sd),
                          random = RE, R = phylocor)
summary(rept.npp.sd.env) # no clear support for the seasonality hypothesis

# save results
saveRDS(rept.npp.sd.env,"Results/BergmannsRule_results_MR_rept_nppsd_env.rds")
rm(rept.npp.sd.env)

# load fitted models
rept.npp.env <- readRDS('Results/BergmannsRule_results_MR_rept_npp_env.rds')
rept.npp.sd.env <- readRDS('Results/BergmannsRule_results_MR_rept_nppsd_env.rds')

# save results in list
re.mr.env <- list(rept.npp.env,rept.npp.sd.env) 
names(re.mr.env) <- c("sd.npp","sd.npp.sd")

saveRDS(re.mr.env,'Results/BergmannsRule_results_MR_rept_env.rds')



