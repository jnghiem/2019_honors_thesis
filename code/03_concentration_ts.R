#03 Peristaltic Pump Sediment Concentration Time Series
#This script processes the data from peristaltic pump samples for use in the settling velocity computation.

library(data.table)
library(dplyr)
library(ggplot2)
library(magrittr)

#Processing the data
data <- fread("C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\1115pumpdata.csv", data.table=FALSE) %>% #reading in the data
  select(-2, -6, -7, -9, -10) #selecting the important fields
names(data) <- c("time", "location", "height", "vol", "mass") #giving more succinct field names
times <- rep(seq(from=-600, by=300, length.out=nrow(data)/6), times=6)
data <- data %>%
  #filter(!is.na(time)) %>% #removing the blanks
  arrange(location, height) %>%
  mutate(time=times, height=height/100, vol=vol/1000, mass=mass) %>% #height in m, volume in L, mass in g
  mutate(mvc=mass/vol) #mass concentration in g/L, same as kg/m3

concentrationdata <- data %>%
  select(height, location, time, mvc) %>% #selecting the important fields
  filter(time>=0) %>%
  arrange(height, location, time) #rearranging the rows for readibility

ggplot(concentrationdata, aes(time, mvc))+geom_line(aes(col=location))+facet_grid(height~.)

                                                    