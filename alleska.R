library("actuar")
library("fitdistrplus")
#library("VGAM")
library("evir")


setwd("V:/commercial/dk/BI/Afdeling/Anders/Speciale")
data <- read.csv("brandproperty.csv")


## Remove claims that are either NA or zero
data_pos <- data[data$skadeudgift>0 & !is.na(data$skadeudgift),]

## Simplify reference
udg <- data_pos$skadeudgift

## Source function calls
source("Functions_Codan.R")

## Histograms
hist(udg,500,xlab = "Claim size")
hist(udg[udg > 1e6],50,main = "Histogram of udg > 1e6",xlab = "Claim size")
hist(udg[udg > quantile(udg,probs = 0.95)],50,
     main = "Histogram of udg > 95% quantile",
     xlab = "Claim size")



## Data looks log-gamma!
hist(log(udg),20, xlab = "Log claim size")



## ---- fit ----

if("VGAM" %in% search()) detach("package:VGAM",unload = TRUE)

## Lognormal fit
fitln <- fitdist(udg,"lnorm")

loglike_pareto <- function(par,data){
    n <- length(data)
    a <- par[1]
    b <- par[2]
    
    loglike <- n*log(a)+a*n*log(b)-sum((a+1)*log(data+b))
    -loglike
}

## Pareto fit
fitpar <- optim(c(1,1),loglike_pareto,data = udg,lower = c(1e-6,1e-6),method = "L-BFGS-B")




fitpar2 <- c(min(udg),length(udg)/(sum(log(udg))-log(min(udg))))

fitgam <- MASS::fitdistr(log(udg),"gamma")

fitlgam <- fitdist(udg,"lgamma",start = list(shapelog=1,ratelog = 1))

fitexp <- MASS::fitdistr(log(udg/min(udg))[log(udg/min(udg))>0],"exponential")



## ---- qqplots ----



par(mfrow = c(1,2))
myqqplot(actuar::qpareto,fitpar$par,udg)
myqqplot(actuar::qpareto,c(0.95,1),udg,TRUE)

myqqplot(VGAM::qpareto,c(398,1),udg)

myqqplot(qlpareto,fitpar$par,log(udg))



## Data is not Pareto (log should then be exponential)
myqqplot(qexp,fitexp$estimate,log(udg/min(udg)))

## Best fit: Is data loggamma distributed?
## Det er det smukkeste fit jeg nogensinde har set
par(mfrow = c(1,2))
myqqplot(qlgamma,fitlgam$estimate,udg)
myqqplot(qgamma,fitgam$estimate,log(udg))


## Log normal
myqqplot(qlnorm,fitln$estimate,udg)




par(mfrow = c(1,1))
evir::meplot(udg)
evir::meplot(log(udg))

set.seed(9748)
par(mfrow = c(2,2))
sapply(rep(1e5,4), function(n){
    evir::meplot(rlgamma(n,fitgam$estimate[1],fitgam$estimate[2]),
                 main = "Mean excess Log gamma sim")
    abline(a = 0,b = fitgam$estimate[2])
    invisible(NULL)
})





fitlgam$estimate


hill(data_pos$skadeudgift,end = 100)
abline(h = 1.3)
abline(h = 2)


qqplot(actuar::qpareto(ppoints(500),2,1), data_pos$skadeudgift)
qqline(data_pos$skadeudgift, distribution = function(p) actuar::qpareto(p,2,1))


meplot2()








