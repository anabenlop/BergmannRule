# BergmannRule

This repository contains all the data and scripts to reproduce the results of the manuscript: 

Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Unveiling the environmental drivers of intraspecific body size variation in terrestrial vertebrates. Global Ecology and Biogeography

The data used for this study are available at https://dx.doi.org/10.6084/m9.figshare.17042993 and in this GitHub repository. 

The main dataset is stored in the file: BergmannsRule_data_forCorrelations_20211114.rds. This dataset contains the raw data downloaded from VertNet and the data pulled from the literature. 

Following the sequence of scripts in the repository from 01_correlation_analysis.R to 14_FigS2.R, these data are cleaned, correlation coefficients and phylogenetic correlation matrixes are calculated, the data for birds and mammals are merged with migratory info, and all the phylogenetic meta-analyses are run. There are some scripts to produce the figures of the paper, or to calculate some summary stats. All scripts are commented.

The final datasets used in the analyses are stored in the folder Final data in the github repository. These files are: mamdata_ph.csv, birddata_ph.csv, reptdata_ph.csv, amphdata_ph.csv. birddata_ph_tarsus.csv corresponds to the dataset using tarsus length as an alternative metric for double checking the consistency in our results.

Some intermediate files are used in the workflow:

Migratory behaviour
SpeciesList3_1_migbehav_v2_0.csv corresponds to the dataset on bird migratrory behaviour by Eyres et al. (2017), used to classify bird species as migratory or not. Nomadic species species were removed from the analysis using the script 05_join_bird_migrationinfo.R

Diet
B_traits_guild.csv and M_traits_guild.csv are based on the EltonTraits database (Wilman, H. et al. EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals: Ecological Archives E095-178. Ecology 95, 2027-2027 (2014)) and used to classify species as carnivores or not. For reptiles we used diet information from Scharf, I. et al. Late bloomers and baby boomers: ecological drivers of longevity in squamates and the tuatara. Global Ecol. Biogeogr. 24, 396-405 (2015). and Meiri, S. Traits of lizards of the world: Variation around a successful evolutionary design. Global Ecol. Biogeogr. 27, 1168-1172 (2018).
