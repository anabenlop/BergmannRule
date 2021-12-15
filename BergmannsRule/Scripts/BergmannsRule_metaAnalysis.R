## Meta-analysis for Bergmann's Rule Project##
## Created on 15 September 2020 ##
## Modified on 02 November 2021 ##

# This script takes the correlations calculated in the script
# 02_BergmannsRule_correlationAnalysis.R and performs a mixed-effects
# meta-analysis with order and family as random effects.
# The model estimates a grand mean correlation coefficient for each tested
# hypothesis (see Methods).

# Packages and working directory -----------------------------------------------

library(metafor)

# setwd('D:/BergmannsRule_upload/BergmannsRule_toSend')

# clean environment
rm(list=ls())

# 1. Amphibians ----------------------------------------------------------------

# read in correlation results
results <- readRDS('Results/BergmannsRule_results_correlations_20211114.rds')

# subset results for amphibians
amphibians <- subset(results, class == 'amphibian')

# vector of environmental variables
env.vars <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

# for loop running a meta-analysis for each environmental variable
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('am.',env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(amphibians, env.var == env.vars[i]),
                random = list(~1|family/speciesname)))
}

# save results for amphibians
am.ma <- list(am.tavg, am.tmin, am.tmax, am.prec, 
              am.pet, am.npp, am.npp.sd)
names(am.ma) <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

saveRDS(am.ma,
        'Results/BergmannsRule_results_MA_amphibians_20211215.rds')

# remove objects
#rm(list=ls(pattern="am"))
#rm(env.vars,i)


# 2. Reptiles ------------------------------------------------------------------

# read in correlation results
results <- readRDS('Results/BergmannsRule_results_correlations_20211114.rds')

# subset results for reptiles
reptiles <- subset(results, class == 'reptile')

# vector of environmental variables
env.vars <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

# for loop running a meta-analysis for each environmental variable
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('rep.',env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(reptiles, env.var == env.vars[i]),
                random = list(~1|family/speciesname)))
}

# save results for reptiles
re.ma <- list(rep.tavg, rep.tmin, rep.tmax, rep.prec, 
              rep.pet, rep.npp, rep.npp.sd)
names(re.ma) <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

#saveRDS(re.ma,
#        'Results/BergmannsRule_results_MA_reptiles_20211115.rds')

# remove objects
#rm(list=ls(pattern="rep."))
#rm(env.vars,i)


# 3. Mammals -------------------------------------------------------------------

# read in correlation results
results <- readRDS('Results/BergmannsRule_results_correlations_20211114.rds')

# subset results for mammals
mammals <- subset(results, class == 'mammal')

# read elton traits dataset and remove marine mammmals
elton_mam <- read.csv("Data/EltonTraits_Mammals_taxid.csv", header = T, stringsAsFactors = F)

mammals <- left_join(mammals, elton_mam[,c("Scientific", "ForStrat.Value")], by = c("speciesname" = "Scientific"))
mammals <-mammals[mammals$ForStrat.Value != "M",]

# vector of environmental variables
env.vars <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

# for loop running a meta-analysis for each environmental variable
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('mam.',env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(mammals, env.var == env.vars[i]),
                random = list(~1|order/family/speciesname)))
}

# save results
ma.ma <- list(mam.tavg, mam.tmin, mam.tmax, mam.prec, 
              mam.pet, mam.npp, mam.npp.sd)
names(ma.ma) <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

#saveRDS(ma.ma,
#        'Results/BergmannsRule_results_MA_mammals_2021214.rds')

# remove objects
#rm(list=ls(pattern="mam."))
#rm(env.vars,i)


# 4. Birds ---------------------------------------------------------------------

# read in correlation results
results <- readRDS('Results/BergmannsRule_results_correlations_20211114.rds')

# subset results for birds
birds <- subset(results, class == 'bird')

# read elton traits dataset and remove marine mammmals
elton_bird <- read.csv("Data/BirdFuncDat.csv", header = T, stringsAsFactors = F)

birds <- left_join(birds, elton_bird[,c("Scientific", "PelagicSpecialist")], by = c("speciesname" = "Scientific"))
birds <-birds[c(birds$PelagicSpecialist == 0 | is.na(birds$PelagicSpecialist)),]
birds <- birds[birds$family != "Pelecanidae", ]

# vector of environmental variables
env.vars <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

# for loop running a meta-analysis for each environmental variable
for(i in 1:length(env.vars)){
  print(i)
  assign(paste0('bir.', env.vars[i]),
         rma.mv(yi = z.cor.yi,
                V = z.cor.vi,
                data = subset(birds, env.var == env.vars[i]),
                random = list(~1|order/family/speciesname)))
}

# save results
bi.ma <- list(bir.tavg, bir.tmin, bir.tmax, bir.prec, 
              bir.pet, bir.npp, bir.npp.sd)
names(bi.ma) <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

saveRDS(bi.ma,
       'Results/BergmannsRule_results_MA_birds_20211214.rds')

# remove objects
#rm(list=ls(pattern="bird."))
#rm(env.vars,i)
