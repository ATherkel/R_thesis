# install.packages("TSA")
# install.packages("xts")


## ---- Time_series_analysis ----

setwd("V:/commercial/dk/BI/Afdeling/Anders/Speciale")
source("Functions_Codan.R")
source_work("Mean_excess.R")

library("xts")
library("TSA")


barplot(udg,names.arg = data_pos$SKADATO)

par(mfrow = c(1,2))
acf(udg,main = "ACF of claim sizes")
pacf(udg,main = "PACF of claim sizes")
par(mfrow = c(1,1))

udg_subset <- armasubsets(udg,3,3,"Claim")
plot(udg_subset)



udg_armafit <- arma(udg,c(2,1))
udg_armafit

udg_sim <- arima.sim(udg_armafit,length(udg))
plot(udg_sim)

udg_diff <- diff(udg)

ts.plot(udg_diff)

par(mfrow = c(1,2))
acf(udg_diff,main = "ACF of claim sizes")
pacf(udg_diff,main = "PACF of claim sizes")
par(mfrow = c(1,1))

plot(armasubsets(udg_diff,5,5,"Claim"))
udg_diff_armafit <- arma(udg,c(3,3))

udg_diff_sim <- arima.sim(udg_diff_armafit,length(udg_diff))
ts.plot(udg_diff_sim)


