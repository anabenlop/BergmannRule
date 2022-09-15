# BergmannRule

This repository contains all the data and scripts to reproduce the results of the manuscript: 

Henry, E., Santini, L., Huijbregts, M. A. J., Benítez-López, A. Unveiling the environmental drivers of intraspecific body size variation in terrestrial vertebrates. Global Ecology and Biogeography.

The data used for this study are available at https://dx.doi.org/10.6084/m9.figshare.17042993 and in this GitHub repository (https://github.com/anabenlop/BergmannRule/). 

The main dataset is stored in the file: BergmannsRule_data_forCorrelations_20211114.rds. This dataset contains the raw data downloaded from VertNet and the data pulled from the literature. 

Following the sequence of scripts in the repository from 01_correlation_analysis.R to 14_FigS2.R, these data are cleaned, correlation coefficients and phylogenetic correlation matrixes are calculated, the data for birds and mammals are merged with migratory info, and all the phylogenetic meta-analyses are run. There are some scripts to produce the figures of the paper, or to calculate some summary stats. All scripts are commented.

The final datasets to run the phylogenetic meta-analyses are stored in the Data folder:
- Amphibians: amphdata_ph.csv,  
- Reptiles: reptdata_ph.csv, 
- Herps: herpdata_ph.csv
- Birds: birddata_ph.csv,  
- Mammals: mamdata_ph.csv,  

The final datasets to run the phylogenetic meta-regressions to test adherence to the hypotheses of sedentary and migratory species stored in the Data folder:
- Birds: birds_ph_mig.csv
- Mammals: mammals_ph_mig.csv

Some intermediate files are used in the workflow:

**Migratory behaviour**

SpeciesList3_1_migbehav_v2_0.csv corresponds to the dataset on bird migratrory behaviour by Eyres et al. (2017), used to classify bird species as migratory or not. Species synonyms were matched with the package taxize. Nomadic species species were removed from the analysis. These analyses were run using the script 05_join_bird_migrationinfo.R

The datasets: Soriano-Redondo_mig_behaviour.csv .... are used to classify mammal species as migratory or not. Species synonyms were matched with the package taxize. These analyses were run using the script 05_join_mammal_migrationinfo.R

**Diet**

BirdFuncDat.csv and EltonTraits_Mammals_taxid.csv are based on the EltonTraits database (Wilman, H. et al. EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals: Ecological Archives E095-178. Ecology 95, 2027-2027 (2014)) and used to classify species as carnivores or not. 

**Family-level phylogenies**
Family-level phylogenies were built to assess the taxonomic coverage of our datasets (see Fig. S4.1 in Appendix S4). Species-level trees were pruned and collapsed to family level using the scripts 11_Obtain_family-level_tree.R and 12_Taxonomic_coverage.R. Species-level families were obtained from: 

- Amphibians: Jetz and Pyron (2018).
- Reptiles: Tonini et al. (2016).
- Birds: AVONET, Tobias et al. (2022), based on the global bird phylogeny by Jetz et al. (2012).
- Mammals: PHYLACINE, Faurby et al. (2018).




