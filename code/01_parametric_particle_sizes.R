#01 Parametric (Normal Distribution) Particle Size Estimation from Sieve Analysis
#This script fits a normal distribution CDF to sieve analysis data. In effect, it is a parametric estimate of particle size distribution.
library(magrittr)
library(dplyr)
library(ggplot2)

#Input (information from Composition Materials Sieve Analysis)
size <- data.frame(microns=c(420, 297, 250, 74), passing=c(1, 0.985, 0.95, 0.04)) #for 60/200 walnut shell
#size <- data.frame(microns=c(595, 420, 149, 74), passing=c(1, 0.975, 0.06, 0.0075)) #for 40/100 walnut shell

#Setup
size <- rbind(size, c(0, 0))
x <- size$microns #assigning the columns of size to individual vectors
y <- size$passing
loss <- function(par, x) sum((pnorm(x, par[1], par[2])-y)^2) #defining a quadratic loss function

#Estimate some guesses
param_est <- function(pair_df) { #a function to estimate mean and standard deviation of normal distn given two points on the CDF
  pair <- pair_df$microns
  prob <- pair_df$passing
  sigma <- diff(pair)/diff(qnorm(prob))
  mu <- pair[1]-qnorm(prob[1])*sigma
  return(c(mu, sigma))
}
size_trunc <- filter(size, !(passing %in% c(0, 1))) #removing the 0 and 1 probability points because they are not well-defined
uniq <- combn(nrow(size_trunc), 2, simplify=FALSE) #all unique combinations of rows from edited data frame

estimates <- data.frame() #initializing a data frame to store parameter estimates for each pair
for (i in 1:length(uniq)) { #for loop to compute parameter estimates
  estimates <- size_trunc[uniq[[i]],] %>%
    param_est() %>%
    rbind(estimates)
}
names(estimates) <- c("mu", "sigma") #renaming columns with better names

parameters <- optim(mapply(median, estimates), loss, lower=floor(mapply(min, estimates)), upper=ceiling(mapply(max, estimates)), method="L-BFGS-B", x=x)$par
#optimize by minimizing the RSS using the initial values and bounds from the rough estimates

#Plot the CDF
plot_func <- function(x) pnorm(x, mean=parameters[1], sd=parameters[2]) #a function to plot the CDF
curve(plot_func, from=0, to=max(x), xlab=expression(Particle~diameter~(mu~m)), ylab="Probability") #plotting the analytical CDF
abline(h=c(0, 0.5, 0.84, 1), lty=2) #plotting horizontal lines marking 0, D50 (median), D84, and 1
points(size, col="red", cex=2) #plotting the sieve analysis data

#Plot the PDF
pdf_func <- function(x) dnorm(x, mean=parameters[1], sd=parameters[2]) #a function to plot the PDF
curve(pdf_func, from=0, to=max(x), xlab=expression(Particle~diameter~(mu~m)), ylab="Density") #plotting the analytical PDF
abline(v=parameters[1], col="red", lty=2) #plotting a vertical line for the mean
abline(h=0, v=0, col="gray") #plotting gray lines for the axes
abline(v=c(abs(diff(parameters)), sum(parameters)), col="blue", lty=2) #plotting the intervals within one standard deviation
abline(v=c(parameters[1]-2*parameters[2], parameters[1]+2*parameters[2]), col="purple", lty=2) #plotting the intervals within two standard deviations
text(x=parameters[1], y=pdf_func(parameters[1]), labels=paste0("Mean = ", round(parameters[1], 4)), pos=4) #plotting label for the mean

rss <- sum((y-plot_func(x))^2) #computing the residual sum of squares
tss <- sum((y-mean(y))^2) #computing the total sum of squares
rsq <- 1-(rss/tss) #computing the coeffcient of determination to assess model fit

#Making a better plot
##Plot of the CDF
#png("C:\\Users\\Bearkey\\Documents\\honors_thesis\\images\\cdf.png", width=1500, height=1000, res=300)
ggplot(data.frame(x=seq(0, 420, length.out=300)))+
  stat_function(fun=plot_func, aes(x=x))+
  geom_point(data=filter(size, passing!=0), aes(x=microns, y=passing), col="red", size=2, inherit.aes=FALSE)+
  geom_hline(yintercept=c(0.5, 0.84), col="blue", lty=2)+
  theme_bw()+
  theme(axis.title=element_blank())
#dev.off()
##Plot of the PDF
ggplot(data.frame(x=seq(0, 330, length.out=300)))+
  stat_function(fun=pdf_func, aes(x=x))+
  geom_vline(xintercept=c(parameters[1]-parameters[2], sum(parameters)), col="blue", lty=2)+
  geom_vline(xintercept=parameters[1], lty=2)+
  theme_bw()+
  theme(axis.title=element_blank())