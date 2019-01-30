#04 Settling Velocity Computation

library(magrittr)
library(dplyr)

#Defining constants
rho_w <- 998.2071 #approximate water density at 20 degrees centigrade in kg/m3
rho_s <- 1300 #approximate walnut shell density in kg/m3
kv <- 1.0023*10^(-6) #approximate kinematic viscosity of water at 20 degrees centigrade in m2/s
P <- 2 #Powers roundness index, 2 is slightly angular/subangular
csf <- 0.5 #Corey shape factor, 0.5 is intermediary between flatness and roundness

#Create size classes and find their means
clnum <- 32 #set the number of classes
bound <- round(parameters[1]*2)/clnum #computes the size of each class
int_func <- function(x) x*dnorm(x, mean=parameters[1], sd=parameters[2]) #a function to integrate to find the expectation
class_mean <- function(up) { #a function to integrate and then normalize to find the mean size within certain bounds
  num <- integrate(int_func, lower=bound*(up-1), upper=bound*up)$value
  den <- pnorm(bound*up, mean=parameters[1], sd=parameters[2])-pnorm(bound*(up-1), mean=parameters[1], sd=parameters[2])
  return(num/den)
}

dstar_func <- function(d) ((rho_s-rho_w)*9.81*d^3)/(rho_w*kv^2) #a function to compute D* (from Dietrich 1982)
wstar_func <- function(dstar) { #a function to compute W* (from Dietrich 1982)
  if (dstar<0.05) {
    return((dstar^2)/5832)
  } else if (dstar>5*10^9) {
    return(NA)
  } else {
    r1 <- -3.76715+1.92944*log10(dstar)-(0.09815*(log10(dstar))^2)-(0.00575*(log10(dstar))^3)+(0.00056*(log10(dstar))^4)
  }
  r2 <- (log10(1-((1-csf)/0.85)))-(((1-csf)^2.3)*tanh(log10(dstar)-4.6))+(0.3*(0.5-csf)*(log10(dstar)-4.6)*(1-csf)^2)
  r3exp <- 1+((3.5-P)/2.5)
  r3 <- (0.65-((csf/2.83)*tanh(log10(dstar)-4.6)))^r3exp
  return(r3*10^(r1+r2))
}
ws_func <- function(wstar) { #a function to compute Ws (from Dietrich 1982)
  ftr <- ((rho_s-rho_w)*9.81*kv)/rho_w
  return((ftr*wstar)^(1/3))
}

size_class <- data.frame(class=1:clnum, upper=rep(bound, times=clnum) %>% cumsum(), mean=sapply(1:clnum, class_mean)) %>% #data by class number and mean particle diameter
  mutate(dstar=dstar_func(mean*10^(-6)), lower=upper-bound) %>% #applying a function to compute D*
  select(class, lower, upper, mean, dstar)

size_class <- data.frame(size_class, wstar=sapply(size_class$dstar, wstar_func)) %>% #applying a function to compute W*
  mutate(ws=ws_func(wstar)) #applying a function to compute settling velocity ws
plot(ws~class, data=size_class) #plotting the relationship between size class and settling velocity