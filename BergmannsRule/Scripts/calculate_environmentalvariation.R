##############################################################
# Authors: 
# Erin Henry
# Email: erinhenry55@gmail.com


## Description ##

# This script calculates environmental variation for each species. It returns
# a table containing each species in the data set, the standard deviation,
# and the range of each environmental variable.

# Script created on 5 May 2022

############################################################## 

# Working directory ------------------------------------------------------------

# working directory
# setwd("M:/Bergmann's Rule revisions/environmental variation")


# Read in data -----------------------------------------------------------------
dat <- read.csv("Bergmanns_bodysize.csv")


# Prep output ------------------------------------------------------------------

# create output table
envvar <- data.frame(species = unique(dat$speciesname))

# vector of environmental variable names
env <- c("tavg", "tmin", "tmax", "prec","pet","npp","npp.sd")


# Calculate environmental variation --------------------------------------------

for(i in 1:length(envvar$species)) { # for each species:
  
  # subset of species i
  spsub <- subset(dat,speciesname == envvar$species[i])

  # calculate sd for each env. variable
  for(j in 1:length(env)){
    print(j)
    sdev <- sd(spsub[,env[j]])
    envvar[i,paste0("sd.",env[j])] <- sdev
  }
  
  # calculate range of each env. variable
  for(j in 1:length(env)){
    rng <- max(spsub[,env[j]]) - min(spsub[,env[j]])
    envvar[i,paste0("rng.",env[j])] <- rng
  }
}

# save -------------------------------------------------------------------------
#write.csv(envvar,"Bergmann_envvariation.csv")

