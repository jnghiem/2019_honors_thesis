#04 Settling Velocity Computation

library(magrittr)
library(dplyr)

setwd("C:\\Users\\Bearkey\\Documents\\honors_thesis\\code") #setting the working directory
source("01_parametric_particle_sizes.R") #calling the script to estimate the particle size distribution

#Defining constants
P <- 2.5 #Powers roundness index, 2.5 is slightly angular/subangular
csf <- 0.5 #Corey shape factor, 0.5 is intermediary between flatness and roundness

#Create size classes and find their means
clnum <- 32 #set the number of classes
bound <- round(out[1]*2)/clnum #computes the size of each class
int_func <- function(x) x*dnorm(x, mean=out[1], sd=out[2]) #a function to integrate to find the expectation
class_mean <- function(up) { #a function to integrate and then normalize to find the mean size within certain bounds
  num <- integrate(int_func, lower=bound*(up-1), upper=bound*up)$value
  den <- pnorm(bound*up, mean=out[1], sd=out[2])-pnorm(bound*(up-1), mean=out[1], sd=out[2])
  return(num/den)
}

dstar_func <- function(d) ((rho_s-rho_w)*9.81*d^3)/(rho_w*kv^2) #a function to compute D* (from Dietrich 1982)

size_class <- data.frame(class=1:clnum, upper=rep(bound, times=clnum) %>% cumsum(), mean=sapply(1:clnum, class_mean)) %>% #data by class number and mean particle diameter
  mutate(dstar=dstar_func(mean*10^(-6)), lower=upper-bound) %>%
  select(class, lower, upper, mean, dstar)

