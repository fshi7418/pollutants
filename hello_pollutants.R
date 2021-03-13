data_dir <-  r"(C:\Users\Frank Shi\Documents\FrankS\Waterloo\pollutants)"
setwd(data_dir)
pollutants <- read.csv(file='pollutants.csv')

summary(pollutants)
