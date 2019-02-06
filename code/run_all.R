#Run all scripts

keep <- c("parameters", "sheardata", "wsr_results", "size_class", "concentrationdata", "estimated_parameters")

setwd("C:\\Users\\Bearkey\\Documents\\honors_thesis\\code")
source("01_parametric_particle_sizes.R")
source("02_shear_velocity.R")
source("03_location_comparison.R")
source("04_settling.R")
source("07_treatment_comparisons.R")
rm(list=ls()[!(ls() %in% keep)])
source("06_shields.R")
source("08_mass_balance.R")
source("09_ks.R")