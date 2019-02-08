#02 Shear Velocity Calculation
#This script calculates the shear velocity at various sampled points in the flume using the law of the wall.
#Conditions are 30 Hz flow rate, corresponding to a discharge of 0.013564 m3/s.

library(R.matlab)
library(data.table)
library(dplyr)
library(stringr)
library(magrittr)

files <- list.files("C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\vectrino profile experiment 181205", pattern="\\.mat$", full.names=TRUE) #all data files
z0 <- (parameters[1]+2.5*parameters[2])*10^(-6) #rough estimate of the characteristic length of bed roughness, taken to be the particle diameter two standard deviations above the mean from the estimated distribution

metadata <- fread("C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\vectrino profile experiment 181205\\positions.csv", data.table=FALSE) %>%
  select(1:4) #reading in metadata
names(metadata) <- c("dowels", "relative_position", "height", "number") #renaming to better names
metadata <- mutate(metadata, height=height/100) #converting cm to m

sheardata <- data.frame() #initializing a data frame to store shear velocity calculations
temperature <- c() #the temperature
for (i in 1:length(files)) { #for loop to iterate over each data file
  file <- files[i]
  fileno <- str_match(file, pattern=" ([[:digit:]]+)\\.mat$")[,2] %>% as.numeric()
  row <- filter(metadata, number==fileno)
  data <- readMat(file)
  velocities <- unlist(data$Data[[3]]) %>% unname()
  u <- mean(velocities)
  temperature <- c(temperature, data$Data[[26]] %>% unlist() %>% unname())
  add_df <- data.frame(row, shearv=unname((0.4*u)/log(row[,"height"]/z0))) #calculation according to the law of the wall (units are m/s)
  sheardata <- rbind(sheardata, add_df)
}
sheardata <- arrange(sheardata, number) #rearranging the rows for better readibility
mean_temp <- mean(temperature) #mean of temperatures