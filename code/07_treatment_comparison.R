#07 Treatment Comparison
library(magrittr)
library(dplyr)
library(ggplot2)

#Read in the pump data
files <- list.files("C:\\Users\\Bearkey\\Documents\\honors_thesis\\data\\pump", full.names=TRUE)
data <- lapply(files, readRDS)
treatments <- rep(seq(from=1, to=length(files)), times=sapply(data, nrow))
data <- plyr::rbind.fill(data) %>%
  mutate(treatment=treatments)

time_limit <- data %>%
  group_by(treatment) %>%
  summarise(max_time=max(time)) %>%
  magrittr::extract(,2) %>%
  min()

data <- filter(data, time<=time_limit)

#Plot for a quick visualization
ggplot(data, aes(time, mc, col=as.character(treatment)))+
  geom_point()+
  geom_vline(xintercept=c(5400, 6000), col="orange")

estimated_parameters <- data.frame()
#Estimate concentration parameters
conc_est <- function(tm) { #a function to return a function with the estimated parameters
  data_fit <- filter(data, treatment==tm)
  coefs <- lm(log_mc~time, data=data_fit) %>% coef() %>% unname()
  estimated_parameters <<- rbind(estimated_parameters, data.frame(initial=exp(coefs[1]), k=-coefs[2], treatment=tm))
  return(function(x) exp(coefs[1])*exp(coefs[2]*x))
}
conc_fn <- lapply(1:length(files), conc_est) #a list of functions to model concentration, one for each experiment

#Compute the coefficient of determination
rsq_fn <- function(tm) {
  data_fit <- filter(data, treatment==tm) %>%
    mutate(fitted_mc=conc_fn[[tm]](time))
  rss <- sum((data_fit$mc-data_fit$fitted_mc)^2)
  tss <- sum((data_fit$mc-mean(data_fit$mc))^2)
  return(1-(rss/tss))
}
rsq <- sapply(1:length(files), rsq_fn)

#Integrate to find the mean difference between curves from 5400 s to 6000 s (final 10 min of experiment)
intgr_fn <- function(initial, k, t0=5400, tf=time_limit) {
  coef <- initial/k
  first <- exp(-k*t0)
  second <- exp(-k*tf)
  return(coef*(first-second))
}

tm_intgr <- c()
for (i in 1:length(files)) {
  params <- filter(estimated_parameters, treatment==i)
  tm_intgr <- c(tm_intgr, intgr_fn(params[,"initial"], params[,"k"]))
}

tm_mean_diff <- -diff(tm_intgr)/(time_limit-5400) #the time-averaged difference in concentration between the experiments for the final 10 min

diff_inf <- estimated_parameters %>%
  mutate(inf_int=initial/k) %>%
  magrittr::extract(,"inf_int") %>%
  diff() %>%
  abs() #the difference of the concentrations integrated from t=0 to infinity, not really interpretable but maybe a useful metric?

#Final visualization
#png("C:\\Users\\Bearkey\\Documents\\honors_thesis\\images\\timeseries.png", width=1500, height=1000, res=300)
ggplot(data, aes(time, mc, col=as.character(treatment)))+
  geom_point(show.legend=FALSE, alpha=0.2)+
  #geom_vline(xintercept=c(5400, 6000), col="purple", lty=2)+
  stat_function(data=data.frame(x=seq(0, 6000, length.out=300)), fun=conc_fn[[1]], col="red", aes(x=x), inherit.aes=FALSE)+
  stat_function(data=data.frame(x=seq(0, 6000, length.out=300)), fun=conc_fn[[2]], col="blue", aes(x=x), inherit.aes=FALSE)+
  labs(x="Time (s)", y="Mass concentration (g/L)")+
  scale_color_manual(values=c("red", "blue"))+
  theme_bw()+
  theme(axis.title=element_blank())
#dev.off()

#Find the mean of paired differences of upstream and downstream concentrations
data %>%
  filter(treatment==2) %>%
  select(-log_mc) %>%
  tidyr::spread(key=location, value=mc) %>%
  mutate(diff=U-D) %>%
  extract(,"diff") %>%
  mean()
