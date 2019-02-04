library(dplyr)
library(magrittr)
library(data.table)

file <- c("C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\1115pumpdata.csv",
          "C:\\Users\\Bearkey\\Documents\\Ecogeomorphic_Flume\\esdlflume\\data\\raw\\1205pumpdata.csv")
data <- data.frame()
for (i in 1:length(file)) {
  add_df <- fread(file[i], data.table=FALSE) %>%
    select(2:4, 9)
  names(add_df) <- c("time", "location", "height", "mc")
  add_df <- add_df %>%
    filter(time>=3, is.finite(mc), !is.nan(mc)) %>%
    tidyr::spread(key=location, value=mc) %>%
    select(U, D) %>%
    filter(!is.na(U), !is.na(D))
  data <- rbind(data, add_df)
}

wilcox.test(x=data$U, y=data$D, alternative="greater", paired=TRUE)
#difference between upstream and downstream seem insignificant