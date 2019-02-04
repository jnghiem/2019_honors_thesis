#08 Location Comparison

#Wilcoxon signed rank test with just the dowel and grease treatment
wsr <- function(file) {
  data <- readRDS(file) %>%
    select(-log_mc) %>%
    tidyr::spread(key=location, value=mc) %>%
    mutate(diff=U-D) %>%
    filter(!is.na(diff))
  return(wilcox.test(x=data$U, y=data$D, alternative="greater", paired=TRUE))
}
wsr("C:\\Users\\Bearkey\\Documents\\honors_thesis\\data\\pump/11152018pumpdata.rds") #differences are insignificant

#Wilcoxon signed rank test with just the no dowel treatment
wsr("C:\\Users\\Bearkey\\Documents\\honors_thesis\\data\\pump/10192018pumpdata.rds") #differences are insignificant

#Wilcoxon signed rank test with both experiments
files <- list.files("C:\\Users\\Bearkey\\Documents\\honors_thesis\\data\\pump", full.names=TRUE)
data <- lapply(files, readRDS)
treatments <- rep(seq(from=1, to=length(files)), times=sapply(data, nrow))
data <- plyr::rbind.fill(data) %>%
  mutate(treatment=treatments) %>%
  select(-log_mc) %>%
  tidyr::spread(key=location, value=mc) %>%
  mutate(diff=U-D) %>%
  filter(!is.na(diff))
wilcox.test(x=data$U, y=data$D, alternative="greater", paired=TRUE) #differences are insignificant