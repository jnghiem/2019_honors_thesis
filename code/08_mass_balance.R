#08 Mass Balance

library(fields)
library(magrittr)
library(dplyr)
library(data.table)
library(raster)

#Upstream region
rho_w <- 998.2071 #approximate water density at 20 degrees centigrade in kg/m3
rho_s <- 1300 #approximate walnut shell density in kg/m3
kv <- 1.0023*10^(-6) #approximate kinematic viscosity of water at 20 degrees centigrade in m2/s
P <- 2 #Powers roundness index, 2 is slightly angular/subangular
csf <- 0.5 #Corey shape factor, 0.5 is intermediary between flatness and roundness

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

ws_d <- function(d) {
  d %>%
    dstar_func() %>%
    wstar_func() %>%
    ws_func() %>%
    return()
}
ws_vectorized <- function(d) sapply(d, ws_d)

curve(ws_vectorized, from=0, to=330*10^(-6))

#Test section
st <- fread("C:\\Users\\Bearkey\\Documents\\honors_thesis\\data\\sediment_trap\\11152019_sediment_trap.csv", data.table=FALSE) %>%
  filter(mass>0)
##Thin-plate spline interpolation
x_res <- 3
y_res <- 3
sp <- Tps(x=st[,c("x", "y")], Y=st$mass/(pi*0.65^2))
eval <- expand.grid(x=seq(x_res/2, by=x_res, length.out=195/x_res), 
                   y=seq(y_res/2, by=y_res, length.out=60/y_res)) %>%
  as.matrix()
predicted <- data.frame(eval, pred=predict.Krig(sp, x=eval)) %>%
  arrange(desc(y)) %>%
  mutate(pred=pred/1000)
ras <- raster(nrows=60/y_res, ncols=195/x_res, vals=predicted$pred*x_res*y_res) #interpolated raster of sediment trap data (each pixel is 3 cm x 3 cm)
plot(ras)

total_settled <- sum(predicted$pred*x_res*y_res) #total settled mass in the test section in g

##Simple averaging
total_settled <- mean(st$mass) %>%
  divide_by(1000*pi*0.65^2) %>%
  multiply_by(60*195)