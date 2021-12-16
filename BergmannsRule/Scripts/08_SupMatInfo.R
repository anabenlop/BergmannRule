##############################################################
# Authors: 
# Erin Henry
# Email: erinhenry55@gmail.com

# Script to create Tables in Supp. Mat, Appendix, and explore data for  the Bergmann's rule paper
# Made 18 October 2020
# Adapted on 16 December 2021


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