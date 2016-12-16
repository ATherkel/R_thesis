## Mean excess plotting and other extreme value analytic stuff
## such as PBdH plot and PoT method

# setwd("~/1_Speciale/R_thesis/")


library("sfsmisc")

## ---- data_load ----

# bachfire <- readRDS("~/Year4/blok4/BachSkade/brand_cleaned.rds")
# udg <- with(bachfire,udg[udg>0])
#udg[udg == min(udg)] <- 1



# setwd("V:/commercial/dk/BI/Afdeling/Anders/Speciale")
setwd("c:/Users/xatk/Desktop/Speciale/R")
data <- read.csv("../brandproperty.csv")

## Sort data on claim date to catch potential dependencies
data$SKADATO <- as.Date(data$SKADATO,"%d/%m/%Y")
data <- data[order(data$SKADATO),]

## Remove claims that are either NA or zero
data_pos <- data[data$skadeudgift>0 & !is.na(data$skadeudgift),]

## Simplify reference
udg <- data_pos$skadeudgift


sudg <- sort(udg,decreasing = TRUE)


## ENDREAD

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
meplot3(udg,axes = 0,ylab = "")
meplot3(udg,panel.first = {abline(v = v <- 5e6, lty = 2,col = 2)},
        add = TRUE,axes = 0,ylab = "")
eaxis(1)
eaxis(2)
box()
meplot3(udg,from = v,main = bquote(bold("Mean Excess above" ~ .(pretty10(v)))),axes = 0,
        ylab = "")
eaxis(1)
eaxis(2)
box()
par(mfrow = c(1,1))





## ---- fit ----
if("package:VGAM" %in% search()) detach("package:VGAM",unload = TRUE)
library("actuar")
library("fitdistrplus")

pos <- rep(1e-2,2)

## fitpar is naive! It uses ALL data, when in reality, we only
## suspect the *tail* of the distribution to be pareto-like
fitpar <- fitdist(udg,distr = "pareto",
                  start = c(shape = 1, scale = 1),lower = pos)
fitlg <- fitdist(udg,distr = "lgamma",
                 start = c(shapelog = 1,ratelog = 1),lower = pos)
fitln <- fitdist(udg,distr = "lnorm")
fitwe <- fitdist(udg,distr = "weibull",lower = pos)



## Log expenses
fitgam <- fitdist(log(udg), distr = "gamma",lower = pos)
fitnorm <- fitdist(log(udg), distr = "norm")




## Test cases:
## Pareto
optim(c(2,1),parloglike,data = udg, method = "L-BFGS-B",
      lower = pos)
## log gamma
optim(c(2,1),lgloglike,data = udg, method = "L-BFGS-B",
      lower = pos)

fitpar2 <- fitdist(udg[udg>1e6],"pareto",
                   start = c(shape = 1, scale = 1),lower = pos)
fitpar2




qqplot3(qpareto,fitpar2$estimate,udg[udg>5e6],omit = 2)



## ---- fit_qqplots ----

# # par(mfrow = c(2,2))
# qqcomp(fitpar,fitcol = "black", 
#        main = "Q-Q plot Pareto",addlegend = FALSE)
# qqcomp(fitlg, fitcol = "black", 
#        main = "Q-Q plot Log gamma",addlegend = FALSE)
# qqcomp(fitgam,fitcol = "black", 
#        main = "Q-Q plot Gamma (log claim)",addlegend = FALSE)
# qqcomp(fitln, fitcol = "black", 
#        main = "Q-Q plot Log normal",addlegend = FALSE)
# par(mfrow = c(1,1))

par(mfrow = c(2,2))
qqplot2(qpareto,fitpar$estimate,udg)
qqplot2(qlgamma,fitlg$estimate,udg)
qqplot2(qlnorm,fitln$estimate,udg)
qqplot2(qweibull,fitwe$estimate,udg)

par(mfrow = c(1,2))
qqplot2(qgamma,fitgam$estimate,log(udg))
qqplot2(qnorm,fitnorm$estimate,log(udg))
par(mfrow = c(1,1))

## ---- Histlog_gaussgamma ----

a <- fitgam$estimate[1]
b <- fitgam$estimate[2]

par(mfrow = c(1,1))
hist(log(udg),100,freq = FALSE)
curve((dgamma(x,fitgam$estimate[1],fitgam$estimate[2])),add = TRUE,col = 4)
curve((dnorm(x,fitnorm$estimate[1],fitnorm$estimate[2])),add = TRUE,col  =2)
par(xpd=TRUE)
legend("topright",legend = c("Gamma","Gaussian"),col = c(4,2),lty = 1,
       bty = "n",adj = c(0, 0.6), cex=0.75)
par(xpd=FALSE)
par(mfrow = c(1,1))

c(samplemean = mean(log(udg)), fitmean = unname(a/b))
c(samplevar = var(log(udg)), fitvar = unname(a/b^2))



## ---- Histogram_gammafit_logdata ----


## simulation
sim <- lapply(rep(length(udg),3), rgamma, a,b)

par(mfrow = c(2,2))
hist(log(udg),100,main = "Histogram log data")
sapply(sim, hist, breaks = 100, main = "Histogram gamma fit simulation")
par(mfrow = c(1,1))


## ---- pot_test ----
## peaks over threshold (ch. 6.5.1 EKM / page 352)

par(mfrow = c(1,2))
barplot(udg)
abline(h = bigval <- 2e6,lty = 2)

udgpot <- ifelse(udg>bigval,udg-bigval,0)
barplot(udgpot)

sum(udgpot>bigval)
par(mfrow = c(1,1))

hist(udgpot[udgpot>0],50,freq = FALSE)

fitpar3 <- fitdist(udgpot[udgpot>0],"pareto",start = c(shape = 1,scale = 1),
        lower = pos)
curve(dpareto(x,fitpar3$estimate[1],fitpar3$estimate[2]),add = TRUE)


udgpot2 <- sample(udgpot,size = 1e6,replace = TRUE)

fitpar4 <- fitdist(udgpot2[udgpot2>0],"pareto",start = c(shape = 1,scale = 1),
                   lower = pos)





n <- 0

n <- n+1
sim <- rnorm(1e7)
simval <- abs(sim[abs(sim)>3])

par(mfrow = c(2,1))
evir::meplot(simval,main = 1,omit = 0)
asd <- par()
VGAM::meplot(simval,ylim = asd$usr[3:4])




        

## ---- lgamma_me_sim ----

set.seed(507) ## my room number

par <- c(shapelog = 34.89786, ratelog = 3.11133)

n <- 1e4 ## Simulation length
m <- 8   ## Number of simulation paths



lgamsim <- lapply(rep(n,m),actuar::rlgamma,
                  shapelog = par[1], ratelog = par[2])


if(!exists("opar")) opar <- par()
par(mfrow=c(2,2),oma=c(0,0,2,0))
par(mar = c(3,4,2,2)+0.1)
sapply(lgamsim[1:4],meplot3)
title(main = "Log gamma simulations",outer = TRUE)
par(opar)

## ---- evt_plots ----
# plot(hillest(udg,length(udg)*0.1),type = "l",ylim = c(0,4),yaxt = "n")
# abline(h = xi <- 1.3,lty = 2,col = 2)
# 
# legend(legend = bquote(xi==.(xi)),"topright",bty = "n",col = 2,lty = 2,xpd = TRUE,inset = c(0.05,0))

par(mfrow = c(2,1))
hill2(udg,ylim = c(0.5,2),xlim = c(15,length(udg)*0.2),auto.scale = FALSE)
abline(h = xi <- 1.1,lty = 2,col = 2)
legend(legend = bquote(xi==.(xi)),"topright",bty = "n",col = 2,lty = 2,xpd = TRUE,inset = c(0.05,0))

hill2(udg,option = "quantile",p = 0.992,end = length(udg)*0.2)

par(opar)
par(mfrow = c(2,1))
plot(pickand(udg,10:430),type = "l",main = "Pickands plot")
abline(h = xi,lty = 2,col = 2)
plot(dedh(udg,10:430),type = "l",ylim = c(-2,1.3),main = "DEdH plot")
abline(h = c(xi,0),lty = 2,col = c(2,4))




## ---- qqplot_bestfit ----

par(mfrow = c(1,2))
qqplot2(qlgamma,par,udg,namespace = "actuar",omit = 3)
qqplot2(qgamma,par,log(udg))
par(mfrow = c(1,1))






