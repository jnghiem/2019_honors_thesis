#Run all scripts

keep <- c("parameters", "sheardata", "size_class", "concentrationdata")

setwd("C:\\Users\\Bearkey\\Documents\\honors_thesis\\code")
source("01_parametric_particle_sizes.R")
source("02_shear_velocity.R")
source("03_concentration_ts.R")
source("04_settling.R")
source("05_rouse.R")
rm(list=ls()[!(ls() %in% keep)])
