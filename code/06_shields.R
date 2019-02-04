#06 Shields Criterion Analysis

#Defining constants
sv <- 1.4284171 #shear velocity im m/s using the ADV experiments from 12/5/2018 in the upstream reach
kv <- 1.0023*10^(-6) #approximate kinematic viscosity of water at 20 degrees centigrade in m2/s
rho_w <- 998.2071 #approximate water density at 20 degrees centigrade in kg/m3
rho_s <- 1300 #approximate walnut shell density in kg/m3

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
  theme(panel.background=element_rect(fill="white"), panel.grid.major=element_line(color="gray"), panel.border=element_rect("black", fill=NA))
taub <- function(d, taustar) (rho_s-rho_w)*9.81*d*taustar #bed shear stress for initial motion in N/m2
