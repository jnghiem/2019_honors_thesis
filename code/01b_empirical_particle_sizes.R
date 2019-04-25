#01b Particle Size Distribution from LISST Portable
#This script plots data of walnut shell size distribution obtained from LISST Portable.

#Packages
library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)

#Reading in data
data <- fread("C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\LISST portable size distn experiment 190422\\WSD101.ASC", data.table=FALSE)
sizes <- !(names(data) %>% as.numeric() %>% is.na()) #find the columns with distribution data
distn <- data[sizes] #extract the data
total_conc <- data[,"Total Volume Concentration"][2] %>% as.numeric() #total volume concentration, to normalize data points
d84 <- data[,"D84"][2] %>% as.numeric() #D84, for use with the shear velocity estimation
d50 <- data[,"D50"][2] %>% as.numeric()
sizedistn <- data.frame(median_size=names(distn) %>% as.numeric(), #data frame with size distribution data
                        vol_concentration=t(distn[2,])[,1] %>% as.numeric()) %>%
  mutate(density=vol_concentration/total_conc)
cumsizedistn <- mutate(sizedistn, cumulative=cumsum(density)) #data frame with cumulative size distribution
#A logarithmic plot of particle size distribution
#png("C:\\Users\\Bearkey\\Documents\\honors_thesis\\images\\pdf_empirical.png", width=1500, height=1000, res=300)
ggplot(sizedistn, aes(median_size, density))+
  geom_point(size=0.8)+
  geom_line()+
  scale_x_log10()+
  annotation_logticks()+
  theme_bw()+
  theme(axis.title=element_blank(), panel.grid.minor=element_blank())
#dev.off()
#A plot of cumulative size distribution
#png("C:\\Users\\Bearkey\\Documents\\honors_thesis\\images\\cdf_empirical.png", width=1500, height=1000, res=300)
ggplot(cumsizedistn, aes(median_size, cumulative))+
  geom_point(size=0.8)+
  geom_line()+
  geom_hline(yintercept=c(0.5, 0.84), lty=2)+
  scale_x_log10()+
  annotation_logticks(sides="b")+
  theme_bw()+
  theme(axis.title=element_blank(), panel.grid.minor.x=element_blank())
#dev.off()

#Stokes law to estimate settling velocity of particle at the critical threshold for initiation of motion
stokes <- function(d) {
  drho <- 300 #kg/m^3, density contrast
  dyn <- 1.002*10^(-3)
  g <- 9.8 #m/s^2, gravitational acceleration
  r2 <- ((d*10^(-6))/2)^2
  return((2/9)*(drho/dyn)*g*r2)
}
curve(stokes(x)*6000, from=0, 16.18)
