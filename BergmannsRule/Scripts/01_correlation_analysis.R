##############################################################
# Authors: 
# Erin Henry, Ana Benitez-Lopez (@anabenlop)
# Email: erinhenry55@gmail.com, abenitez81@gmail.com

## Correlation Analysis for Bergmann's Rule project ##
## Created in October 2020 ##
# Adapted on 24 December 2021


# This script calculates correlation coefficients for each species.
# For each species, we tested the correlation of 
# mean body size (per 5min grid cell) with each environmental variable. 
# The correlation coefficients are then transformed
# into Fisher's z-scores for the meta-analysis (see Methods).

# Output is a data.frame containing each species, taxonomic information, 
# correlation coefficients between body size and an environmental variable, 
# the environmental variable tested, Fishers z-scores and associated confidence
# intervals and sampling variance.

# Packages and working directory -----------------------------------------------
library(plyr)
library(metafor)
library(wCorr)

# setwd("D:/BergmannsRule_upload")

# Clean environment
rm(list = ls())

# 1. Read in prepared data -----------------------------------------------------

# Dataset with mean body size per 5 minute grid cell for each species
dat <- readRDS('Data/BergmannsRule_data_forCorrelations_20211114.rds')
head(dat)
str(dat)

# Fix some species names that are synonyms
dat[dat$speciesname == "Speotyto cunicularia", "speciesname"] <- "Athene cunicularia"
dat[dat$speciesname == "Parus bicolor", "speciesname"] <- "Baeolophus bicolor"
dat[dat$speciesname == "Myiothlypis rivularis", "speciesname"] <- "Basileuterus rivularis"
dat[dat$speciesname == "Nyctea scandiaca", "speciesname"] <- "Bubo scandiacus"
dat[dat$speciesname == "Pycnonotus virens", "speciesname"] <- "Eurillas virens"
dat[dat$speciesname == "Dendragapus canadensis", "speciesname"] <- "Falcipennis canadensis"
dat[dat$speciesname == "Pipra coronata", "speciesname"] <- "Lepidothrix coronata"
dat[dat$speciesname == "Ceryle alcyon", "speciesname"] <- "Megaceryle alcyon"
dat[dat$speciesname == "Otus asio", "speciesname"] <- "Megascops asio"
dat[dat$speciesname == "Parkesia motacilla", "speciesname"] <- "Seiurus motacilla"
dat[dat$speciesname == "Ammodramus sandwichensis", "speciesname"] <- "Passerculus sandwichensis"
dat[dat$speciesname == "Columba fasciata", "speciesname"] <- "Patagioenas fasciata"
dat[dat$speciesname == "Hirundo pyrrhonota", "speciesname"] <- "Petrochelidon pyrrhonota"
dat[dat$speciesname == "Ceyx picta", "speciesname"] <- "Ispidina picta"
dat[dat$speciesname == "Halcyon sanctus", "speciesname"] <- "Todiramphus sanctus"
dat[dat$speciesname == "Microtus mogollonensis", "speciesname"] <- "Microtus mexicanus"
dat[dat$speciesname == "Dicrostonyx kilangmiutak", "speciesname"] <- "Dicrostonyx groenlandicus"

# save data with fixed species names
write.csv(dat, "Data/Bergmanns_bodysize.csv", row.names = F)

# 2. Perform correlations for all species --------------------------------------

# vector of environmental variables
env.vars <- c('tavg','tmin','tmax','prec','pet','npp','npp.sd')

# correlations for all species (between env.var and body mass/SVL)
for (i in 1:length(env.vars)) {
  
  print(i)
  
  # make summary table of all species in data
  results <- plyr::count(dat, vars = c('speciesname','class','order','family'))
  
  # empty column for environmental variable and correlation coefficient
  results$env.var <- env.vars[i]
  results$corr.coeff <- NA
  
  # correlations for each species - weightedCorr()
  for(j in 1:nrow(results)) {
    
    if (j%%100 == 0) {print(j)} # benchmark
    
    # species subset
    sub <- subset(dat, speciesname == results[j,'speciesname'])
    
    # ensure sample size is numeric
    sub$Sample_size <- as.numeric(sub$Sample_size)
    
    # NA Sample size -> 1
    sub[which(is.na(sub$Sample_size)),'Sample_size'] <- 1 
    
    # define x for correlations (current env. predictor)
    x <- env.vars[i]
    
    # define y: size variable (ectotherms -> length; endotherms -> mass)
    if(results[j,'class'] == 'amphibian' | results[j,'class'] == 'reptile') {
      
      y <- 'avg.length'
      
    } else {
      
      y <- 'avg.mass'
      
    }
    
    # weighted correlation (weight = log10(Sample_size + 1))
    results[j,'corr.coeff'] <- wCorr::weightedCorr(sub[,x],
                                              sub[,y],
                                              method='spearman',
                                              weights=log10(sub$Sample_size+1))
  }
  
  # class as factor
  results$class <- factor(results$class)
  
  # add to output table
  if(i==1){
    corr.results <- results
    
  } else {
    corr.results <- rbind(corr.results, results)
  }
}


# 3. Calculate Fisher's z-scores -----------------------------------------------

# If correlations equal 1 or -1, convert to .999 or -.999. Otherwise escalc 
# is not able to calculate the z-scores.

# correlation = 1 -> 0.999999
corr.results[which(corr.results$corr.coeff == 1), 
                    'corr.coeff'] <- 0.9999999

# correlation = -1 -> -0.999999
corr.results[which(corr.results$corr.coeff == -1), 
                    'corr.coeff'] <- -0.9999999

# Some rows with correlation coefficients equal to 1 or -1 are not registering.
# Change these manually:

corr.results[which(corr.results$speciesname == 'Alophoixus ochraceus' & 
                   corr.results$env.var == 'prec'), 
             'corr.coeff'] <-  0.9999999

corr.results[which(corr.results$speciesname == 'Muscicapa muttui' & 
                     corr.results$env.var == 'tavg'), 
             'corr.coeff'] <-  0.9999999

corr.results[which(corr.results$speciesname == 'Phoenicurus ochruros' & 
                     corr.results$env.var == 'tavg'), 
             'corr.coeff'] <-  0.9999999

corr.results[which(corr.results$speciesname == 'Phoenicurus ochruros' & 
                     corr.results$env.var == 'tmax'), 
             'corr.coeff'] <-  0.9999999

corr.results[which(corr.results$speciesname == 'Picus canus' & 
                     corr.results$env.var == 'prec'), 
             'corr.coeff'] <-  0.9999999

corr.results[which(corr.results$speciesname == 'Inezia caudata' & 
                     corr.results$env.var == 'npp'), 
             'corr.coeff'] <-  0.9999999

corr.results[which(corr.results$speciesname == 'Phoenicurus ochruros' & 
                     corr.results$env.var == 'npp'), 
             'corr.coeff'] <-  0.9999999

corr.results[which(corr.results$speciesname == 'Chlidonias niger' & 
                     corr.results$env.var == 'pet'), 
             'corr.coeff'] <-  0.9999999

corr.results[which(corr.results$speciesname == 'Limosa haemastica' & 
                     corr.results$env.var == 'npp.sd'), 
             'corr.coeff'] <-  0.9999999

corr.results[which(corr.results$speciesname == 'Uraeginthus bengalus' & 
                     corr.results$env.var == 'prec'), 
             'corr.coeff'] <-  0.9999999

corr.results[which(corr.results$speciesname == 'Uraeginthus bengalus' & 
                     corr.results$env.var == 'npp'), 
             'corr.coeff'] <-  0.9999999

corr.results[which(corr.results$speciesname == 'Tyrannus couchii' & 
                     corr.results$env.var == 'prec'), 
             'corr.coeff'] <-  -0.9999999


# calculate z-scores
z.cor <- metafor::escalc("ZCOR",
                         ri = corr.results$corr.coeff,
                         ni = corr.results$freq) # weight by number of cells
corr.results$z.cor.yi = z.cor$yi # Fisher's z-score
corr.results$z.cor.vi = z.cor$vi # sampling variance

# 4. Save results --------------------------------------------------------------

saveRDS(corr.results,
        'Results/BergmannsRule_results_correlations_20211224.rds')

write.csv(corr.results,
          'Results/BergmannsRule_results_correlations_202111224.csv')


# End of script -----