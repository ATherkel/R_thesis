## Mean excess plotting and other extreme value analytic stuff
## such as PBdH plot and PoT method

# setwd("~/1_Speciale/R_thesis/")




## ---- data_load ----

bachfire <- readRDS("~/Year4/blok4/BachSkade/brand_cleaned.rds")
udg <- with(bachfire,udg[udg>0])
#udg[udg == min(udg)] <- 1

sudg <- sort(udg,decreasing = TRUE)



## ---- functions ----
## Functions that are supposed to be generic are written in the
## functions.R file. Contains e.g. own modified versions of plot
## amongst others.

source("Functions_Codan.R")




## ---- Pre_plot ----
## What does the data look like, no assumptions made whatsoever?

par(mfrow = c(1,2))
barplot(udg,main = "Claim sizes")
hist(log(udg),100,main = "Histogram of log claim sizes",
     xlab = "log claim",freq = FALSE)
par(mfrow = c(1,1))



## ---- me_plot ----
par(mfrow = c(1,2))
meplot3(udg);abline(v = v <- 5.5e5)
meplot3(udg,from = v)
par(mfrow = c(1,1))



## ---- fit ----
if("package:VGAM" %in% search()) detach("package:VGAM",unload = TRUE)
library("actuar")
library("fitdistrplus")
fitpar <- fitdist(udg,distr = "pareto")
fitlg <- fitdist(udg,distr = "lgamma",start = c(1,1))
fitln <- fitdist(udg,distr = "lnorm")
fitcauchy <- fitdist(udg,distr = "cauchy")

## Log expenses
fitgam <- fitdist(log(udg), distr = "gamma")
fitnorm <- fitdist(log(udg), distr = "norm")


qqplot(qexp(ppoints(500)),a <- rpareto(1e3,2,1))
b <- fitdist(a,"exp")

plot(b)
qqplot2(qexp,b$estimate,a)

## ---- fit_qqplots ----

# par(mfrow = c(2,2))
qqcomp(fitpar,fitcol = "black", 
       main = "Q-Q plot Pareto",addlegend = FALSE)
qqcomp(fitlg, fitcol = "black", 
       main = "Q-Q plot Log gamma",addlegend = FALSE)
qqcomp(fitgam,fitcol = "black", 
       main = "Q-Q plot Gamma (log claim)",addlegend = FALSE)
qqcomp(fitln, fitcol = "black", 
       main = "Q-Q plot Log normal",addlegend = FALSE)
# par(mfrow = c(1,1))

par(mfrow = c(2,2))
qqplot2(qpareto,fitpar$estimate,udg)
qqplot2(qlgamma,fitlg$estimate,udg)
qqplot2(qlnorm,fitln$estimate,udg)
qqplot2(qcauchy,fitcauchy$estimate,udg)

par(mfrow = c(1,2))
qqplot2(qgamma,fitgam$estimate,log(udg))
qqplot2(qnorm,fitnorm$estimate,log(udg))
par(mfrow = c(1,1))

## ---- test ----

a <- fitgam$estimate[1]
b <- fitgam$estimate[2]

par(mfrow = c(1,1))
hist(log(udg),100,freq = FALSE)
curve((dgamma(x,fitgam$estimate[1],fitgam$estimate[2])),add = TRUE,col = 4)
curve((dnorm(x,fitnorm$estimate[1],fitnorm$estimate[2])),add = TRUE,col  =2)
legend("topright",legend = c("Gamma fit","Gaussian fit"),col = c(3,2),lty = 1)
par(mfrow = c(1,1))

c(samplemean = mean(log(udg)), fitmean = unname(a/b))
c(samplevar = var(log(udg)), fitvar = unname(a/b^2))

## simulation
sim <- lapply(rep(length(udg),3), rgamma, a,b)

par(mfrow = c(2,2))
hist(log(udg),100,main = "Histogram log data")
sapply(sim, hist, breaks = 100, main = "Histogram gamma fit simulation")
par(mfrow = c(1,1))


## ---- pot ----
## peaks over threshold (ch. 6.5.1 EKM / page 352)

mecdf <- ecdf((udg))
plot(mecdf)
excesscdf <- function(u,data = udg){
    sum(data>u)/length(data)
}
excesscdf(1e5)

n <- 0

n <- n+1
sim <- rnorm(1e7)
simval <- abs(sim[abs(sim)>3])

par(mfrow = c(2,1))
evir::meplot(simval,main = 1,omit = 0)
asd <- par()
VGAM::meplot(simval,ylim = asd$usr[3:4])




        

## ---- lgamma_me_sim ----

par <- c(shapelog = 34.89786, ratelog = 3.11133)

n <- 1e5 ## Simulation length
m <- 8   ## Number of simulation paths



lgamsim <- lapply(rep(n,m),actuar::rlgamma,
                  shapelog = par[1], ratelog = par[2])


par(mfrow = c(2,2))
sapply(seq_along(lgamsim)[1:4], function(x) {
    meplot2(lgamsim[[x]])
    curve(x/(par[1]-1),add = TRUE)
})
par(mfrow = c(1,1))




par(mfrow = c(1,2))
qqplot2(qlgamma,par,udg,namespace = "actuar")
qqplot2(qgamma,par,log(udg))
par(mfrow = c(1,1))






