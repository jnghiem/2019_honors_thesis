#02 Shear Velocity Calculation
#This script calculates the shear velocity at various sampled points in the flume using the law of the wall.
#Conditions are 30hZ flow rate

library(R.matlab)
library(data.table)
library(dplyr)
library(stringr)
library(magrittr)

setwd("C:\\Users\\Bearkey\\Documents\\honors_thesis\\code") #setting the working directory
source("01_parametric_particle_sizes.R") #calling the script to estimate the particle size distribution
files <- list.files("C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\vectrino profile experiment 181205", pattern="\\.mat$", full.names=TRUE) #all data files
z0 <- (out[1]+3*out[2])*10^(-6) #rough estimate of the characteristic length of bed roughness, taken to be the particle diameter three standard deviations above the mean from the estimated distribution

metadata <- fread("C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\vectrino profile experiment 181205\\positions.csv", data.table=FALSE) %>%
  select(1:4) #reading in metadata
names(metadata) <- c("dowels", "relative_position", "height", "number") #renaming to better names
metadata <- mutate(metadata, height=height/100) #converting cm to m

sheardata <- data.frame() #initializing a data frame to store shear velocity calculations
for (i in 1:length(files)) { #for loop to iterate over each data file
  file <- files[i]
  fileno <- str_match(file, pattern=" ([[:digit:]]+)\\.mat$")[,2] %>% as.numeric()
  row <- filter(metadata, number==fileno)
  data <- readMat(file)
  velocities <- unlist(data$Data[[3]]) %>% unname()
  u <- mean(velocities)
  add_df <- data.frame(row, shearv=unname(0.4/(u*log(row[,"height"]/z0)))) #calculation according to the law of the wall (units are m/s)
  sheardata <- rbind(sheardata, add_df)
}
sheardata <- arrange(sheardata, number) #rearranging the rows for better readibility
