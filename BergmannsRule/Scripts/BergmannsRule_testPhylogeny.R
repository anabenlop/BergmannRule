## Script to test for phylogenetic signal in meta-analyses ##
## Calculates Pagel's lambda for each class: Amphibians, reptiles, mammals, birds
## Modified on 14 November 2021 ##

# Packages and working directory -----------------------------------------------
library(phytools)
library(plyr)
library(metafor)

# setwd('D:/BergmannsRule_upload/BergmannsRule_toSend')

# 1. Amphibians ----------------------------------------------------------------

# data
amph.tree <- read.tree("Data/amph_shl_dates.tre")
amphibians <- readRDS("Results/BergmannsRule_results_correlations_20211114.rds")
amphibians <- subset(amphibians,class=='amphibian')
amphibians$Species_ph <- gsub(" ", "_", trimws(amphibians$speciesname))

# models
amph.mods <- readRDS("Results/BergmannsRule_results_MA_amphibians_20211115.rds")
env.vars <- unique(amphibians$env.var)

for(i in 1:length(env.vars)){
  print(i)
  
  # subset correlation results
  am <- subset(amphibians,env.var==env.vars[i])
  
  #exclude species in the tree that are not in your dataset
  drops<-amph.tree$tip.label[!amph.tree$tip.label %in% am$Species_ph]
  amph.tree<-drop.tip(amph.tree, drops)
  
  if(i==1){length(am[!am$Species_ph %in% amph.tree$tip.label,"Species_ph"])}
  
  #Fit model MOD
  res <- resid(amph.mods[[i]])
  names(res) <- am$Species_ph #name the residuals with the species names (with the underscore)
  assign(paste0("lam.",env.vars[i]),
         phylosig(amph.tree, res, method="lambda", test=TRUE)) #if significant, means there is phylogenetic autocorrelation
}

amph.lambda <- list(lam.tavg,lam.tmin,lam.tmax,lam.prec,lam.pet,lam.npp,lam.npp.sd)
names(amph.lambda) = c("tavg","tmin","tmax","prec","pet","npp","npp.sd")

# save
#saveRDS(amph.lambda,"Results/BergmannsRule_results_phylLambda_amphibians.rds")
rm(amph.tree,amphibians,amph.mods,drops,res,lam.tavg,lam.tmin,lam.tmax,
   lam.prec,lam.pet,lam.npp,lam.npp.sd,am)



# 2. Reptiles ------------------------------------------------------------------

# data
rep.tree <- read.tree("Data/rep_phyl_Pyron.txt") #Pyron & Burbrink 2014
reptiles <- readRDS("Results/BergmannsRule_results_correlations_20211114.rds")
reptiles <- subset(reptiles, class=="reptile")
reptiles$Species_ph <- gsub(" ", "_", trimws(reptiles$speciesname))

# models
rep.mods <- readRDS("Results/BergmannsRule_results_MA_reptiles_20211115.rds")
env.vars <- unique(reptiles$env.var)

# calculate lambda
for(i in 1:length(env.vars)){
  
  # subset correlation results
  re <- subset(reptiles, env.var==env.vars[i])
  
  #exclude species in the tree that are not in your dataset
  drops <- rep.tree$tip.label[!rep.tree$tip.label %in% re$Species_ph]
  rep.tree<-drop.tip(rep.tree, drops)
  
  if(i==1){length(re[!re$Species_ph %in% rep.tree$tip.label,"Species_ph"])}
  
  #Fit model MOD
  res<-resid(rep.mods[[i]])
  names(res)<-re$Species_ph #name the residuals with the species names (with the underscore)
  assign(paste0("lam.",env.vars[i]),
         phylosig(rep.tree, res, method="lambda", test=TRUE)) #if significant, means there is phylogenetic autocorrelation
}

rep.lambda <- list(lam.tavg,lam.tmin,lam.tmax,lam.prec,lam.pet,lam.npp,lam.npp.sd)
names(rep.lambda) <- c("tavg","tmin","tmax","prec","pet","npp","npp.sd")

# save
#saveRDS(rep.lambda,"Results/BergmannsRule_results_phylLambda_reptiles.rds")
rm(rep.tree,reptiles,rep.mods,re,drops,res,lam.tavg,lam.tmin,
   lam.tmax,lam.prec,lam.pet,lam.npp,lam.npp.sd)



# 3. Mammals -------------------------------------------------------------------

# data
mam.tree <- read.nexus("Data/Mammalian philogeny_Fritz.txt")
mam.tree <- mam.tree$mammalST_MSW05_bestDates
mammals <- readRDS("Results/BergmannsRule_results_correlations_20211114.rds")
mammals <- subset(mammals, class=="mammal")
mammals$Species_ph <- gsub(" ", "_", trimws(mammals$speciesname))

# models
mam.mods <- readRDS("Results/BergmannsRule_results_MA_mammals_20211115.rds")
env.vars <- unique(mammals$env.var)

for(i in 1:length(env.vars)){
  
  # subset data by env.var
  ma <- subset(mammals, env.var==env.vars[i])
  
  #exclude species in the tree that are not in your dataset
  drops <- mam.tree$tip.label[!mam.tree$tip.label %in% ma$Species_ph]
  mam.tree <- drop.tip(mam.tree,drops)
  
  if(i==1){length(ma[!ma$Species_ph %in% mam.tree$tip.label,"Species_ph"])}
  
  #Fit model MOD
  res<-resid(mam.mods[[i]])
  names(res)<-ma$Species_ph #name the residuals with the species names (with the underscore)
  assign(paste0("lam.",env.vars[i]),
         phylosig(mam.tree, res, method="lambda", test=TRUE)) #if significant, means there is phylogenetic autocorrelation
}

mam.lambda <- list(lam.tavg,lam.tmin,lam.tmax,lam.prec,lam.pet,lam.npp,lam.npp.sd)
names(mam.lambda) <- c("tavg","tmin","tmax","prec","pet","npp","npp.sd")

# save
#saveRDS(mam.lambda,"Results/BergmannsRule_results_phylLambda_mammals.rds")

# save mam.tree to run phylogenetic meta-analysis
write.tree(mam.tree, file = "Data/mam.tree.tre")

# save mammal data
write.csv(mammals, file = "Data/mammals.csv", row.names = F)

rm(mam.tree,mammals,mam.mods,drops,res,lam.tavg,lam.tmin,lam.tmax,lam.prec,
   lam.pet,lam.npp,lam.npp.sd,ma)



# 4. Birds ---------------------------------------------------------------------

# get models and correlation results (note: correlation results come from 
# a different file that includes synonyms in the species_ph column)
bird.mods <- readRDS("Results/BergmannsRule_results_MA_birds_20211214.rds")
birds <- readRDS("Results/BergmannsRule_results_correlationsBirds_20211115.rds")
birds <- subset(birds,class=='bird')

# read elton traits dataset and remove marine mammmals
elton_bird <- read.csv("Data/BirdFuncDat.csv", header = T, stringsAsFactors = F)

birds <- left_join(birds, elton_bird[,c("Scientific", "PelagicSpecialist")], by = c("speciesname" = "Scientific"))
# birds$PelagicSpecialist <- ifelse(is.na(birds$PelagicSpecialist),0,birds$PelagicSpecialist)

birds <-birds[c(birds$PelagicSpecialist == 0 | is.na(birds$PelagicSpecialist)),]
birds <- birds[birds$family != "Pelecanidae", ]

# get phylogenetic tree
bird.tree.orig <- read.tree("Data/AllBirdsHackett1.tre")
bird.tree <- as.phylo(bird.tree.orig[[1]]) # use first tree

#exclude species in the tree that are not in your dataset
drops<-bird.tree$tip.label[!bird.tree$tip.label %in% birds$Species_ph]
bird.tree<-drop.tip(bird.tree, drops)

# env variables
env.vars <- unique(birds$env.var)

# i <- 1
for(i in 1:length(env.vars)){
  print(i)
  
  # subset results
  bi <- subset(birds, env.var == env.vars[i])
  
  #exclude species in the tree that are not in your dataset
  drops<-bird.tree$tip.label[!bird.tree$tip.label %in% bi$Species_ph]
  bird.tree<-drop.tip(bird.tree, drops)
  
  if(i==1){length(bi[!bi$Species_ph %in% bird.tree$tip.label,"Species_ph"])}
  
  #Fit model MOD
  res<-resid(bird.mods[[i]])
  names(res)<-bi$Species_ph #name the residuals with the species names (with the underscore)
  assign(paste0("lam.",env.vars[i]),
         phylosig(bird.tree, res, method="lambda", test=TRUE)) #if significant, means there is phylogenetic autocorrelation
}

bird.lambda <- list(lam.tavg,lam.tmin,lam.tmax,lam.prec,lam.pet,lam.npp,lam.npp.sd)
names(bird.lambda) = c("tavg","tmin","tmax","prec","pet","npp","npp.sd")

rm(lam.tavg,lam.tmin,lam.tmax,lam.prec,lam.pet,lam.npp,lam.npp.sd,bird.lambda,drops,
   i,k,res,bird.tree)

# save
saveRDS(bird.lambda,"Results/BergmannsRule_results_phylLambda_birds.rds")  

# save mam.tree to run phylogenetic meta-analysis
write.tree(bird.tree, file = "Data/bird.tree.tre")

# save mammal data
write.csv(birds, file = "Data/birds.csv", row.names = F)


# *Test first 10 trees for birds -----------------------------------------------

# get models and correlation results (note: correlation results come from 
# a different file that includes synonyms in the species_ph column)
bird.mods <- readRDS("Results/BergmannsRule_results_MA_birds_20211214.rds")
birds <- readRDS("Results/BergmannsRule_results_correlationsBirds_20211115.rds")
birds <- subset(birds,class=='bird')

# read elton traits dataset and remove marine mammmals
elton_bird <- read.csv("Data/BirdFuncDat.csv", header = T, stringsAsFactors = F)

birds <- left_join(birds, elton_bird[,c("Scientific", "PelagicSpecialist")], by = c("speciesname" = "Scientific"))
# birds$PelagicSpecialist <- ifelse(is.na(birds$PelagicSpecialist),0,birds$PelagicSpecialist)

birds <-birds[c(birds$PelagicSpecialist == 0 | is.na(birds$PelagicSpecialist)),]
birds <- birds[birds$family != "Pelecanidae", ]

# get phylogenetic tree
bird.tree.orig <- read.tree("Data/AllBirdsHackett1.tre")
# bird.tree <- bird.tree.orig[[1]] # use first tree

# Code tests lambda calculated using different trees
for(k in 1:10){
  print(paste0("k = ",k))

# get tree
bird.tree <- bird.tree.orig[[k]]

# env variables
env.vars <- unique(birds$env.var)

for(i in 1:length(env.vars)){
  print(i)
  
  # subset results
  bi <- subset(birds, env.var == env.vars[i])
  
  #exclude species in the tree that are not in your dataset
  drops<-bird.tree$tip.label[!bird.tree$tip.label %in% bi$Species_ph]
  bird.tree<-drop.tip(bird.tree, drops)
  
  if(i==1){length(bi[!bi$Species_ph %in% bird.tree$tip.label,"Species_ph"])}
  
  #Fit model MOD
  res<-resid(bird.mods[[i]])
  names(res)<-bi$Species_ph #name the residuals with the species names (with the underscore)
  assign(paste0("lam.",env.vars[i]),
         phylosig(bird.tree, res, method="lambda", test=TRUE)) #if significant, means there is phylogenetic autocorrelation
}

bird.lambda <- list(lam.tavg,lam.tmin,lam.tmax,lam.prec,lam.pet,lam.npp,lam.npp.sd)
names(bird.lambda) = c("tavg","tmin","tmax","prec","pet","npp","npp.sd")

assign(paste0("bird.lambda",k),
       bird.lambda)

rm(lam.tavg,lam.tmin,lam.tmax,lam.prec,lam.pet,lam.npp,lam.npp.sd,drops,
   i,k,res,bird.tree)
}

