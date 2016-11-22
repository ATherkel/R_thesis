
## ---- hill_output ----

library("actuar")

set.seed(123)
dat <- lapply(1:100,function(x) rpareto(1e4,2,1))
# dat <- lapply(1:100,function(x) runif(1e3)^(-1/2))

dat_hill <- lapply(1:length(dat), function(i) {
    hillest(dat[[i]],1000,10)
})

dat_min <- min(sapply(dat_hill,min))
dat_max <- max(sapply(dat_hill,max))
dat_min <- 1.4
dat_max <- 3

plot2(dat_hill[[1]],type = "l", ylim = c(dat_min,dat_max),
      panel.first = abline(h = 2,col = "red", lty = 2))

cols <- heat.colors(min(length(dat_hill),10)-1)
sapply(2:min(length(dat_hill),10), function(i){
    par(new = TRUE)
    plot(dat_hill[[i]],ylim = c(dat_min,dat_max),
         axes = FALSE,xlab = "",ylab = "",type = "l",
         col = cols[i])
})

vals <- sapply(1:length(dat_hill), function(i){
    points <- dat_hill[[i]][c(200,400,600)]
    names(points) <- c(200,400,600)
    points
})



boxplot(t(vals), main = "Boxplots of Hill value at different levels from the tail")




## ---- mean_excess_output ----

set.seed(1234)
nsim <- 1e5
sim_func <- function(distr = c("exp","lognorm","pareto","gamma"),
                     simlen = 1e5,numsim = 1e2){
    rexp2 <- function(n) rexp(n)
    rlnorm2 <- function(n) rlnorm(n)
    rpareto2 <- function(n) rpareto(n,3,4/5)
    rgamma2 <- function(n) rgamma(7,1)
    
    lapply(setNames(distr,distr), function(d){
        Vectorize(get(paste0("r",distr,"2")))(rep(simlen,numsim))
    })
}
sim <- sim_func()



excessout <- lapply(1:ncol(sim$exp),function(simnum,...){
    meanexcess(data = sim$exp[,simnum])
},plotout = FALSE)

# excessout$pareto <- meanexcess(sim$pareto)
#meanexcess(sim2 <- rexp(1e6),to = max(sim2))


# excess_min <- min(sapply(excessout,min,na.rm = TRUE))
# excess_max <- max(sapply(excessout,max,na.rm = TRUE))
excess_min <- 0
excess_max <- 20

par(mfrow = c(2,2))
sapply(names(excessout), function(name){
    plot2(excessout[[name]], main = paste("Mean excess of", name))
})
par(mfrow = c(1,1))


plot2(excessout[[1]],ylim = c(excess_min,excess_max))

sapply(2:length(excessout), function(i){
    par(new = TRUE)
    plot(excessout[[i]],ylim = c(excess_min,excess_max),
         axes = FALSE, xlab = "", ylab = "", col = i)
})
legend("topright",legend = names(excessout),col = 1:length(excessout),
       pch = 1)


asd <- rlnorm(1e5)
par(mfrow = c(1,2))
a <- meanexcess(asd,from = 0, to = par()$usr[2],length.out = 5e2,xaxs = "i",
                xlim = par()$usr[1:2],ylim = par()$usr[3:4],type = "l",lty = 1)


par(mfrow = c(1,1))
meplot(asd,main = "meplot",xlim = c(0,70),ylim = c(0,20),xaxs = "i",yaxs = "i")
usr <- par()$usr
par(new = TRUE)
plot(a,xlim = usr[1:2],ylim = usr[3:4],col = 2,
     xaxs = "i",yaxs = "i",type = "l")
abline(v = 43,lty = 2)
abline(h = mean(ewq)-43,lty = 2,col = "grey")
qwe <- sort(asd)
qwe[qwe>43] ->ewq
sum(ewq)/43



## ---- eventually_pareto ----

tails <- sapply(c("runif","rnorm","rexp"),function(rfunc){
    dat <- get(rfunc)(1e7)
    dat[which(dat>=quantile(dat,0.99999))]
})
qqplot(qpareto(ppoints(length(tails[,1])),1,10),tails[,1])
