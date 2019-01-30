library(data.table)
library(dplyr)
library(magrittr)

data <- fread("C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\1019pumpdata.csv", data.table=FALSE) %>%
  select(-2, -6, -7, -9, -10)
names(data) <- c("time", "location", "height", "vol", "mass")
density <- 1300 #density for 40/100 and 60/200 walnut shell in kg/m3
data <- data %>%
  filter(!is.na(time)) %>%
  mutate(height=height/100, vol=vol/1000000, mass=mass/1000) %>% #height in m, volume in m3, mass in kg
  mutate(mvc=mass/vol) #mass concentration in kg/m3

data <- filter(data, location=="U", height==0.05) %>%
  select(time, mvc) %>%
  arrange(time)
refh <- 0.05

