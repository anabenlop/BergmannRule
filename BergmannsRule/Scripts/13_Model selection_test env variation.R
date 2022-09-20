##############################################################
# Authors: 
# Erin Henry, Ana Benitez-Lopez (@anabenlop)
# Email: erinhenry55@gmail.com, abenitez81@gmail.com, ana.benitez@ugr.es
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# https://www.anabenitezlopez.com/

##############################################################
# Description of script and instructions
##############################################################

# This script fits the phylogenetic meta-analysis and meta-regressions with environmental variability 
# and performs model selection to test whether species exposed to a wider range of variation in environmental
# conditions are more likely to adhere to the tested hypothesis for mammals 

##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
library(tictoc)


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

# saveRDS(amph.ma,
#         'Results/BergmannsRule_results_MA_amphibians_phylo_nonphylo_ML.rds')

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

# read results
amph.ma <- readRDS('Results/BergmannsRule_results_MA_amphibians_phylo_nonphylo_ML.rds')
am.mr.env <- readRDS('Results/BergmannsRule_results_MR_amph_env_ML.rds')

anova(amph.ma$prec, am.mr.env$sd.prec) # model with env var slightly better
anova(amph.ma$npp, am.mr.env$sd.npp) # null model better
anova(amph.ma$npp.sd, am.mr.env$sd.npp.sd) # null model better

# 2. Reptiles ---------------------------------------------------------------------

#Load data
reptiles_ph <- read.csv("Data/reptdata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/Phylogeny/rept_phylo_cor.Rdata") #rept_phylo_cor

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
                       data = reptiles_ph,
                       subset = env.var=="npp",
                       mods = ~ log10(sd.npp),
                       random = RE, R = phylocor, 
                       method = "ML")
summary(rept.npp.env) # tendency to positive effect (large CI) with variation in npp across the gradient of body size records

# save results
saveRDS(rept.npp.env,"Results/BergmannsRule_results_MR_rept_npp_env.rds")
rm(rept.npp.env)

# npp.sd model with env variation
rept.npp.sd.env <- rma.mv(yi = z.cor.yi,
                          V = z.cor.vi,
                          data = reptiles_ph,
                          subset = env.var=="npp.sd",
                          mods = ~ log10(sd.npp.sd),
                          random = RE, R = phylocor,
                          method = "ML")
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

saveRDS(re.mr.env,'Results/BergmannsRule_results_MR_rept_env_ML.rds')

# read results
rept.ma <- readRDS('Results/BergmannsRule_results_MA_reptiles_phylo_nonphylo_ML.rds')
re.mr.env <- readRDS('Results/BergmannsRule_results_MR_rept_env_ML.rds')

anova(rept.ma$npp, re.mr.env$sd.npp) # null model better, model with env variation tendency to positive effects of NPP
anova(rept.ma$npp.sd, re.mr.env$sd.npp.sd) # null model better

# 3. Birds ---------------------------------------------------------------------

#Load data
birds_ph <- read.csv("Data/birddata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/Phylogeny/bird_phylo_cor.Rdata") #bird_phylo_cor

# vector of environmental variables
env.vars <- c('tavg','tmax','npp','npp.sd')

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = bird_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# for loop running a meta-analysis for each environmental variable
tic("Run phylo meta-analysis in a loop")
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('bird.',env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(birds_ph, env.var == env.vars[i]),
                random = RE, R = phylocor, method ="ML"))
}

toc()

# save results
bird.ma <- list(bird.tavg, bird.tmax, bird.npp, bird.npp.sd)
names(bird.ma) <- env.vars 

saveRDS(bird.ma,
        'Results/BergmannsRule_results_MA_bird_phylo_nonphylo_ML.rds')

# env variation models
# tavg model with env variation
bird.tavg.env <- rma.mv(yi = z.cor.yi,
                       V = z.cor.vi,
                       data = birds_ph,
                       subset = env.var=="tavg",
                       mods = ~ log10(sd.tavg),
                       random = RE, R = phylocor, 
                       method = "ML")
summary(bird.tavg.env) # 

# save results
saveRDS(bird.tavg.env,"Results/BergmannsRule_results_MR_bird_tavg_env.rds")
rm(bird.tavg.env)

# tmax model with env variation
bird.tmax.env <- rma.mv(yi = z.cor.yi,
                        V = z.cor.vi,
                        data = birds_ph,
                        subset = env.var=="tmax",
                        mods = ~ log10(sd.tmax),
                        random = RE, R = phylocor, 
                        method = "ML")
summary(bird.tmax.env) # 

# save results
saveRDS(bird.tmax.env,"Results/BergmannsRule_results_MR_bird_tmax_env.rds")
rm(bird.tmax.env)

# npp model with env variation
bird.npp.env <- rma.mv(yi = z.cor.yi,
                       V = z.cor.vi,
                       data = birds_ph,
                       subset = env.var=="npp",
                       mods = ~ log10(sd.npp),
                       random = RE, R = phylocor, 
                       method = "ML")
summary(bird.npp.env) # tendency to positive effect (large CI) with variation in npp across the gradient of body size records

# save results
saveRDS(bird.npp.env,"Results/BergmannsRule_results_MR_bird_npp_env.rds")
rm(bird.npp.env)

# npp.sd model with env variation
bird.npp.sd.env <- rma.mv(yi = z.cor.yi,
                          V = z.cor.vi,
                          data = birds_ph,
                          subset = env.var=="npp.sd",
                          mods = ~ log10(sd.npp.sd),
                          random = RE, R = phylocor,
                          method = "ML")
summary(bird.npp.sd.env) # no clear support for the seasonality hypothesis

# save results
saveRDS(bird.npp.sd.env,"Results/BergmannsRule_results_MR_bird_nppsd_env.rds")
rm(bird.npp.sd.env)

# load fitted models
bird.tavg.env <- readRDS('Results/BergmannsRule_results_MR_bird_tavg_env.rds')
bird.tmax.env <- readRDS('Results/BergmannsRule_results_MR_bird_tmax_env.rds')
bird.npp.env <- readRDS('Results/BergmannsRule_results_MR_bird_npp_env.rds')
bird.npp.sd.env <- readRDS('Results/BergmannsRule_results_MR_bird_nppsd_env.rds')

# save results in list
bi.mr.env <- list(bird.tavg.env, bird.tmax.env, bird.npp.env,bird.npp.sd.env) 
names(bi.mr.env) <- c("sd.tavg", "sd.tmax","sd.npp","sd.npp.sd")

saveRDS(bi.mr.env,'Results/BergmannsRule_results_MR_bird_env_ML.rds')

# read results
bird.ma <- readRDS('Results/BergmannsRule_results_MA_bird_phylo_nonphylo_ML.rds')
bi.mr.env <- readRDS('Results/BergmannsRule_results_MR_bird_env_ML.rds')

anova(bird.ma$tavg, bi.mr.env$sd.tavg) 
anova(bird.ma$tmax, bi.mr.env$sd.tmax) 
anova(bird.ma$npp, bi.mr.env$sd.npp) 
anova(bird.ma$npp.sd, bi.mr.env$sd.npp.sd) 

# 4. Mammals -----------------------------------------------------------------------
#Load data
mammals_ph <- read.csv("Data/mamdata_ph.csv", stringsAsFactors = F)

# loading phylogenetic matrixes 
load("Data/Phylogeny/mam_phylo_cor.Rdata") #mam_phylo_cor

# vector of environmental variables
env.vars <- c('tavg','tmax','npp','npp.sd')

# define phylo vcov matrix and random effects
phylocor<-list(speciesname  = mam_phylo_cor)
RE = list( ~1|speciesname, ~1|SPID)

# for loop running a meta-analysis for each environmental variable
tic("Run phylo meta-analysis in a loop")
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('mam.',env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(mammals_ph, env.var == env.vars[i]),
                random = RE, R = phylocor, method ="ML"))
}

toc()

# save results
mam.ma <- list(mam.tavg, mam.tmax, mam.npp, mam.npp.sd)
names(mam.ma) <- env.vars 

saveRDS(mam.ma,
        'Results/BergmannsRule_results_MA_mam_phylo_nonphylo_ML.rds')

# env variation models
# tavg model with env variation
mam.tavg.env <- rma.mv(yi = z.cor.yi,
                        V = z.cor.vi,
                        data = mammals_ph,
                        subset = env.var=="tavg",
                        mods = ~ log10(sd.tavg),
                        random = RE, R = phylocor, 
                        method = "ML")
summary(mam.tavg.env) # 

# save results
saveRDS(mam.tavg.env,"Results/BergmannsRule_results_MR_mam_tavg_env.rds")
rm(mam.tavg.env)

# tmax model with env variation
mam.tmax.env <- rma.mv(yi = z.cor.yi,
                        V = z.cor.vi,
                        data = mammals_ph,
                        subset = env.var=="tmax",
                        mods = ~ log10(sd.tmax),
                        random = RE, R = phylocor, 
                        method = "ML")
summary(mam.tmax.env) # 

# save results
saveRDS(mam.tmax.env,"Results/BergmannsRule_results_MR_mam_tmax_env.rds")
rm(mam.tmax.env)

# npp model with env variation
mam.npp.env <- rma.mv(yi = z.cor.yi,
                       V = z.cor.vi,
                       data = mammals_ph,
                       subset = env.var=="npp",
                       mods = ~ log10(sd.npp),
                       random = RE, R = phylocor, 
                       method = "ML")
summary(mam.npp.env) # tendency to positive effect (large CI) with variation in npp across the gradient of body size records

# save results
saveRDS(mam.npp.env,"Results/BergmannsRule_results_MR_mam_npp_env.rds")
rm(mam.npp.env)

# npp.sd model with env variation
mam.npp.sd.env <- rma.mv(yi = z.cor.yi,
                          V = z.cor.vi,
                          data = mammals_ph,
                          subset = env.var=="npp.sd",
                          mods = ~ log10(sd.npp.sd),
                          random = RE, R = phylocor,
                          method = "ML")
summary(mam.npp.sd.env) # no clear support for the seasonality hypothesis

# save results
saveRDS(mam.npp.sd.env,"Results/BergmannsRule_results_MR_mam_nppsd_env.rds")
rm(mam.npp.sd.env)

# load fitted models
mam.tavg.env <- readRDS('Results/BergmannsRule_results_MR_mam_tavg_env.rds')
mam.tmax.env <- readRDS('Results/BergmannsRule_results_MR_mam_tmax_env.rds')
mam.npp.env <- readRDS('Results/BergmannsRule_results_MR_mam_npp_env.rds')
mam.npp.sd.env <- readRDS('Results/BergmannsRule_results_MR_mam_nppsd_env.rds')

# save results in list
ma.mr.env <- list(mam.tavg.env, mam.tmax.env, mam.npp.env,mam.npp.sd.env) 
names(ma.mr.env) <- c("sd.tavg", "sd.tmax","sd.npp","sd.npp.sd")

saveRDS(ma.mr.env,'Results/BergmannsRule_results_MR_mam_env_ML.rds')

# read results
mam.ma <- readRDS('Results/BergmannsRule_results_MA_mam_phylo_nonphylo_ML.rds')
ma.mr.env <- readRDS('Results/BergmannsRule_results_MR_mam_env_ML.rds')

anova(mam.ma$tavg, ma.mr.env$sd.tavg) 
anova(mam.ma$tmax, ma.mr.env$sd.tmax) 
anova(mam.ma$npp, ma.mr.env$sd.npp) 
anova(mam.ma$npp.sd, ma.mr.env$sd.npp.sd) 

# End of script ----

