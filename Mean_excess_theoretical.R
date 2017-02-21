
# setwd("~/1_Speciale/R_thesis/")
source("Functions_Codan.R")
library(actuar)
library(dplyr)


set.seed(125)
sim <- rpareto(1e4, p <- 2,p[2] <- 1e5)


meplot3(sim)
meplot3(qpareto(ppoints(10000),p[1],p[2]),type = "l")

meplot3(qlnorm(ppoints(1e5)),type = "l")
meplot3(qexp(ppoints(1e5)),type = "l")
meplot3(qgamma(ppoints(1e5),3,1),type = "l")
meplot3
evir::meplot(qexp(ppoints(1e3)))


hist(qexp(ppoints(1e3)),100)


library(fExtremes)
fExtremes::mePlot(qexp(ppoints(1e3)))
mePlot(rbeta(1e3,2,2))


me <- function(u,dist,par,upper = Inf,plot = TRUE,...){
    pdist_temp <- get(paste0("p",dist))
    ddist_temp <- get(paste0("d",dist))
    
    
    
    if(!missing(par)){
        if(dist == "pareto" & par[1] <=1) 
            stop(paste("mean excess not defined for chosen parameters"))
        if(length(par) == 1){
            pdist <- function(q,...) pdist_temp(q,par,...)
            ddist <- function(x,...) ddist_temp(x,par,...)
        } else {
            pdist <- function(q,...) pdist_temp(q,par[1],par[2],...)
            ddist <- function(x,...) ddist_temp(x,par[1],par[2],...)
        }
    } else {
        pdist <- pdist_temp
        ddist <- ddist_temp
    }
    pdist <- Vectorize(pdist)
    ddist <- Vectorize(ddist)
    

    integrand <- function(x,u) (x-u)*ddist(x)
    
    outf <- function(u) round(1/(pdist(u,lower.tail = FALSE)) * 
                                  integrate(integrand,u,upper,u)$value,3)
    
    out <- cbind(u,sapply(u,outf))
    
    if(plot == TRUE){
        plot(out,...)
        invisible(out)
    } else if(plot == "lines"){
        lines(out,...)
    } else {
        out
    }
}

{
    reg <- seq(1,50,0.25)
    parm <- list(par = c(25,1),
                 exp = 1,
                 gamma = c(3,1),
                 lgamma = c(35,35),
                 lnorm = c(-500,20))
    
    # layout(t(matrix(c(1,1,1,2),4,2)))
    par(mar = c(5,4,4,13)+.1)
    me(reg,"pareto",parm[["par"]],type = "l")
    me(reg,"exp",parm[["exp"]],col = 2,plot = "lines")
    me(reg,"gamma",parm[["gamma"]],col = 3,plot = "lines")
    me(reg,"lgamma",parm[["lgamma"]],col = 4,plot = "lines")
    # par(new = TRUE)
    me(reg,"lnorm",parm[["lnorm"]],col = 5,plot = FALSE) -> asd
    
    ## Legend
    legn <- c("Pareto","Exponential","Gamma","Log-gamma","Log-normal")
    leg <- paste0(legn,"(",sapply(parm,paste,collapse = ","),")")
    legend("topright",legend = c("Dist+parameters",leg),
           inset = c(-0.35,0),xpd=NA,cex = 1,col = 0:5, lty = 1,
           text.font = c(2,rep(1,5)))
    # plot.new()
    title(main = "Mean excess function")
}



## ---- ggplot2_version ----

# library(ggplot2)
# df <- data.frame(u = 1:50, par = me(1:50,"pareto",c(25,1),plot = FALSE)[,2])
# df$exp <- me(1:50,"exp",plot = FALSE)[,2]
# df$gamma <- me(1:50,"gamma",c(3,1),plot = FALSE)[,2]
# df$lgamma <- me(1:50,"lgamma",c(35,35),plot = FALSE)[,2]
# df$lnorm <- me(1:50,"lnorm",c(-500,20),plot = FALSE)[,2]
# 
# mytheme <- theme_classic(base_size = 12) +
#     theme(panel.border = element_rect(fill = NA,size = 0.5),
#           plot.title = element_text(hjust = 0.5))
#     
# 
# ggplot(df,aes(u)) +
#     mytheme + 
#     ggtitle("asd") + 
#     # theme(plot.margin = unit(c(5,4,4,2)+0.1))+
#     geom_line(aes(y = par),col = 1) +
#     geom_line(aes(y = exp),col = 2) + 
#     geom_line(aes(y = gamma), col = 3) +
#     geom_line(aes(y = lgamma),col = 4) 
#     #geom_line(aes(y = lnorm),col = 5)








meplot4(qpareto(ppoints(500),10,1e5))
abline(a = 11000,b = 0.11)
set.seed(1)
meplot4(sim <- rpareto(1e4,10,10),type = "l")
me(0:10,"pareto",c(10,10),plot = "lines")




plot(meplot4(qexp(ppoints(500))),type = "l")
plot(meplot4(qpareto(ppoints(5000),4,2),omit = 10),type = "l")

