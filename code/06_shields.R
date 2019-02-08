#06 Shields Criterion Analysis

library(ggplot2)
library(magrittr)
library(dplyr)

#Defining constants
sv <- 0.003742604 #shear velocity in m/s using the ADV experiments from 12/5/2018 in the upstream reach
kv <- 1.0023*10^(-6) #approximate kinematic viscosity of water at 20 degrees centigrade in m2/s
rho_w <- 998.2071 #approximate water density at 20 degrees centigrade in kg/m3
rho_s <- 1300 #approximate walnut shell density in kg/m3
bs <- rho_w*sv^2 #computed bed shear stress for flume at 30 hZ (upstream section)

re <- function(d) (sv*d)/kv #a function to compute the particle Reynolds number
Rep <- function(re) re*sqrt((rho_s-rho_w)/rho_w) #a function to compute the transformed particle Reynolds number
taustar <- function(Rep) { #a function to compute the critical Shields stress
  if (Rep<3.16) {
    return(0.135*Rep^(-0.261))
  } else {
    y <- Rep^(-0.6)
    return((0.22*y)+(0.06*exp(-17.77*y)))
  }
}
ts_plot <- function(Rep) sapply(Rep, taustar) #a function to vectorize the function taustar

#A plot of transformed particle Reynolds number vs time
ggplot(data.frame(x=seq(10, 10^4, length.out=1000)))+
  stat_function(fun=ts_plot, aes(x=x))+
  coord_trans(x="log10")+
  ylim(0.025, 0.075)+
  labs(x="transformed particle Reynolds number", y="critical Shields stress")+
  scale_x_continuous(breaks=10^c(1:4), labels=10^c(1:4))+
  theme_bw()
taub <- function(d, taustar) (rho_s-rho_w)*9.81*d*taustar 

#A plot of bed shear stress as a function of particle diameter
shields <- function(d) { #a function to compute bed shear stress for initial motion in N/m2
  ts <- d %>%
    re() %>%
    Rep() %>%
    taustar()
  return((rho_s-rho_w)*9.81*d*ts)
}
shields_vectorized <- function(d) sapply(d, shields)

ggplot(data.frame(x=seq(0, 330*10^(-6), length.out=300)))+
  stat_function(fun=shields_vectorized, aes(x=x))+
  labs(x=expression(Particle~diameter~(mu~m)), y=expression(Critical~bed~shear~stress~(N/m^2)))+
  scale_x_continuous(labels=function(x) x/10^(-6))+
  geom_hline(yintercept=bs, lty=2, col="red")+ #assuming that this upstream bed shear stress is representative
  theme_bw()

#Determine the intersection of the diameter vs shear stress curve and the observed bed shear stress
##Only applicable if the particle diameter is less than 1490 microns:
shields_inv <- function(tau) {
  first <- (sqrt((rho_s-rho_w)/rho_w)*(sv/kv))^0.261
  second <- tau/(0.135*(rho_s-rho_w)*9.81)
  return((first*second)^(1/0.739))
}
shields_inv(bs)*10^6 #the particle diameter in microns that corresponds to the threshold of incipient motion, about 15 microns

png("C:\\Users\\Bearkey\\Documents\\honors_thesis\\images\\shields.png", width=1500, height=1000, res=300)
ggplot(data.frame(x=seq(0, 330*10^(-6), length.out=300)))+
  stat_function(fun=shields_vectorized, aes(x=x))+
  labs(x=expression(Particle~diameter~(mu~m)), y=expression(Critical~bed~shear~stress~(N/m^2)))+
  scale_x_continuous(labels=function(x) x/10^(-6))+
  geom_vline(xintercept=shields_inv(bs), col="red", lty=2)+
  geom_hline(yintercept=bs, lty=2, col="red")+ #assuming that this upstream bed shear stress is representative
  theme_bw()+
  theme(axis.title=element_blank())
dev.off()

#Determine the portion of sediment whose critical shear stress lies below the observed shear stress
pnorm(shields_inv(bs)*10^6, mean=parameters[1], sd=parameters[2])*200 #a very small number (in g)
