library(dplyr)
library(magrittr)
library(ggplot2)
library(data.table)

settings <- data.frame(
  bins=c(1:14, 15:24),
  st=c(4700, 5900),
  ed=c(Inf, 7100)
)

file <- "C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\Flume -- LISST -- Sept 2018-\\20181115\\11_15_particle_size_dist.asc"
correction <- TRUE #set this to TRUE if you want to remove all data points that have a transmission outside of the optimal range of 0.1 to 0.995
time_start <- 0 #set start time for the data analysis to exclude initialization
time_end <- Inf #set end time for data analysis

###Data Processing Section###
bins <- 1:32
raw <- read.table(file, header=FALSE) #%>% mutate(V41=V41/100) #reading the .asc file as a table

laser_reference <- any(raw[,36]==0) #laser reference values should all be much greater than 0

if (laser_reference==TRUE) {
  warning("Laser reference has a 0 value")
}

transmission <- any(c(raw[,41]<=0, raw[,41]>=1)) #it is physically impossible to have a transmission not in the range 0 to 1

if (transmission==TRUE) {
  warning("Transmission has value(s) outside of the range 0 to 1")
}

if (correction==TRUE) {
  raw <- raw %>%
    filter(V36>0, V41>=0.1, V41<=0.995)
}

time_converted_sec <- raw %>%
  mutate(min=as.numeric(ifelse(V40>=100, gsub("[[:digit:]]{2}$", "", V40), 0))) %>%
  mutate(sec=as.numeric(ifelse(V40>=1000, gsub("^[[:digit:]]{2}", "", V40), ifelse(V40<1000 & V40>=100, gsub("^[[:digit:]]", "", V40), V40)))) %>%
  mutate(sec=sec+(min*60)) %>%
  select(-min) #conversion of time units to seconds

diff <- c(0, diff(time_converted_sec[,"sec"]))
diff[diff<0] <- diff[diff<0]+3600
seconds <- cumsum(diff)+1

time_converted_final <- time_converted_sec %>%
  mutate(sec=seconds) #setting seconds to start at 0 and count increasing (as opposed to cyclically)

gathered_measurements <- time_converted_final %>%
  tidyr::gather(key=bin, value=measurement, 1:32) %>%
  mutate(bin=as.numeric(gsub("^V", "", bin))) %>%
  select(sec, bin, measurement) %>% #combining particle concentration data into a single field
  filter(sec>=time_start, sec<=time_end) %>%
  mutate(logc=log(measurement))

grp1 <- gathered_measurements %>%
  filter(bin<=14, sec>=4700) %>%
  mutate(sec=sec-min(sec))

grp2 <- gathered_measurements %>%
  filter(bin>=15, bin<=24, sec>=5900, sec<=7100) %>%
  mutate(sec=sec-min(sec))

data_ed <- rbind(grp1, grp2)

k <- c()
for (i in 1:length(unique(data_ed$bin))) {
  mod <- lm(logc~sec, data=filter(data_ed, bin==i))
  k <- c(k, -coef(mod)[2])
}
k <- unname(k)

bin_info <- fread("C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\lisst bins.tsv", data.table=FALSE)[,c(1, 8:10)]
names(bin_info) <- c("bin", "lower", "upper", "mid")

maxc <- filter(gathered_measurements, bin<=24, sec>=3300, sec<=3420) %>%
  group_by(bin) %>%
  summarise(maxc=mean(measurement)) %>%
  left_join(bin_info, by="bin") %>%
  as.data.frame()

ggplot(maxc, aes(mid, maxc))+geom_point()+coord_trans(x="log10")
ggplot(maxc[-c(1:3, 7, 10),], aes(mid, maxc))+geom_point()+coord_trans(x="log10")

#Try a smoothing spline
maxc_log <- maxc[-c(1:3, 7, 10),] %>%
  mutate(mid_log=log10(mid))

ss <- smooth.spline(x=as.matrix(select(maxc_log, mid_log, maxc)), lambda=0.0001)
plot_ss <- function(x) predict(ss, x)$y
plot(maxc~mid_log, data=maxc_log, col="red")
curve(plot_ss, from=min(log10(maxc$lower)), to=max(log10(maxc$upper)), add=TRUE)
abline(v=log10(165))

maxc <- maxc %>%
  mutate(predc=plot_ss(log10(mid)), k=k) %>%
  mutate(predc=ifelse(predc<0, 0, predc), k=ifelse(k<0, 0, k))

