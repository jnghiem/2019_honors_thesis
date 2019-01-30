#03 Peristaltic Pump Sediment Concentration Time Series
#This script processes the data from peristaltic pump samples for use in the settling velocity computation.

library(data.table)
library(dplyr)
library(magrittr)

#Processing the data
data <- fread("C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\1019pumpdata.csv", data.table=FALSE) %>% #reading in the data
  select(-2, -6, -7, -9, -10) #selecting the important fields
names(data) <- c("time", "location", "height", "vol", "mass") #giving more succinct field names
data <- data %>%
  filter(!is.na(time)) %>% #removing the blanks
  mutate(height=height/100, vol=vol/1000, mass=mass) %>% #height in m, volume in L, mass in g
  mutate(mvc=mass/vol) #mass concentration in g/L, same as kg/m3

concentrationdata <- data %>%
  select(height, location, time, mvc) %>% #selecting the important fields
  arrange(height, location, time) #rearranging the rows for readibility