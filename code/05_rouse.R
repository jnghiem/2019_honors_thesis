#05 Rouse Equation Analysis
#This script is based on a Rouse sediment profile analysis, drawing on Lupker et al. 2011.

library(magrittr)
library(dplyr)

loc <- "U" #location, either upstream (U) or downstream (D)
H <- 0.4
refh <- 0.05

#Rouse number computation
sv <- filter(sheardata, dowels=="F", height==0.07) #the shear velocity in m/s
if (loc=="U") {
  sv <- filter(sv, relative_position<0)[1,"shearv"]
} else if (loc=="D") {
  sv <- filter(sv, relative_position>0)[1,"shearv"]
}
  
rn_func <- function(ws) ws/(0.4*sv)

size_class <- size_class %>%
  mutate(rn=rn_func(ws))

spec_conc <- filter(concentrationdata, location==loc, height==refh)

rouse_eq <- function(z, refcon, rn, H=0.4, refh=0.05) {
  inner_num <- (H-z)/z
  inner_den <- (H-refh)/refh
  inner <- (inner_num/inner_den)^rn
  return(refcon*(inner^rn))
}

rp <- vector("list", length=nrow(size_class))
for (i in 1:nrow(size_class)) {
  rp_t <- vector("list", length=nrow(spec_conc))
  rn <- force(filter(size_class, class==i)[1,"rn"])
  for (j in 1:nrow(spec_conc)) {
    rp_t[[j]] <- function(z) {
      inner_num <- (H-z)/z
      inner_den <- (H-refh)/refh
      inner <- (inner_num/inner_den)^rn
      refcon <- force(spec_conc[j,"mvc"])
      return(refcon*(inner^rn))
    }
  }
  rp[[i]] <- rp_t
}
plot(x=sapply(seq(from=0, to=0.4, by=0.001), rp[[1]][[15]]), y=seq(from=0, to=0.4, by=0.001), type="l")
