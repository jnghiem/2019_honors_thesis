#Need to run the mass balance script (need the variable total_settled)
source("C:\\Users\\Bearkey\\Documents\\honors_thesis\\code\\08_mass_balance.R")
rm(list=ls()[ls()!="total_settled"])
#Need to run the treatment comparison script (need the run data)
source("C:\\Users\\Bearkey\\Documents\\honors_thesis\\code\\07_treatment_comparison.R")
rm(list=ls()[!(ls() %in% c("data", "total_settled", "estimated_parameters"))])

m0 <- 200 #initial mass added to flume (g)
u <- 0.057 #average flow velocity (m/s)
dc <- 0.003175 #collector diameter (m)
lc <- 1450 #approximate stem density (m^-2)
kv <- 1.0023*10^(-6) #approximate kinematic viscosity of water at 20 degrees centigrade in m2/s

ks_comp <- function(t, m0, ms, k) { #based on Jordan's result
  numerator <- ms*k
  denom <- m0*(1-exp(-k*t))
  return(numerator/denom)
}
ks <- ks_comp(max(data$time), m0, settled_test+settled_outside, estimated_parameters[2,"k"])
kc <- estimated_parameters[2,"k"]-ks #subtract to obtain kc
ece_comp <- function(kc, u, dc, lc) kc/(u*dc*lc) #a function to calculate effective capture efficiency (mks units)
ece <- ece_comp(kc, u, dc, lc)*100 #effective capture efficiency (as a percentage)

#Compare to Fauria et al. [2015] model for 7209 m^-2 stem density and biofilm
ece_empir <- function(Rec, R) 2.06*Rec^(-1.14)*R^0.65 #Rec is collector Reynolds number and R is ratio of particle diameter to collector diameter

ece_empir((u*dc)/kv, (165.08*10^(-6))/dc)
