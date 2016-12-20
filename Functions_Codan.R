
## ---- functions ----

plot2 <- function(...){
    plot(...,axes = FALSE,mgp = c(2.5,1,0))
    axis(1, tck = -.015, labels = NA)
    axis(2, tck = -.015, labels = NA)
    axis(1, lwd = 0, line = -.4)
    axis(2, lwd = 0, line = -.4, las = 2)
    box()
}


## Hill estimator
hillest <- function(data,length,first = 1){
    ## Equation (6.25) EKM
    ## data: input dataset
    ## length: number of observations to include
    
    sdata <- sort(data,decreasing = TRUE)
    
    sapply(seq(from = first, to = length-first+1),function(len){
        mean(log(sdata[1:len])) - log(sdata[len])
    })^(-1)
}


which2 <- function(u,data){
    lapply(u, function(x) which(data>x))
}

## Empirical mean excess
meanexcess <- function(data, from = 0, to = max(data)*0.4,
                       length.out = 1e2,plotout = TRUE,...){
    ## Equation (6.6) EKM
    ## data: input dataset
    ## from, to and length.out: sequence parameters for 
    ## the excess levels
    ## plotout: if TRUE, outputs plot. Otherwise outputs values
    ## ...: Additional plot parameters
    
    u_vals <- seq(from,to,length.out = length.out)
    
    delta_n <- which2(u_vals,data)
    
    # duplicates <- sapply(2:length(delta_n), function(i){
    #     if(length(delta_n[[i]]) == length(delta_n[[i-1]])){
    #         i
    #     }
    # })
    
    # delta_n[unlist(duplicates)] <- NA
    
    
    out <- sapply(seq(delta_n), function(i){
        sum(data[delta_n[[i]]] - u_vals[i]) / length(delta_n[[i]])
    })
    
    if(plotout){
        dots <- list(...)
        if(!("ylab" %in% names(dots))){
            dots$ylab <- "e(u)"
        }
        if(!("xlab" %in% names(dots))){
            dots$xlab <- "u"
        }
        do.call(plot2,c(list(x = u_vals,y = out),dots))
        # plot2(u_vals, out,pch = 20,...=unlist(dots))
    }
    invisible(cbind(u = u_vals, "e(u)" = out))
}


## ---- me_plot_def ----

meplot2 <- function (data, from = 0,omit = 3, labels = TRUE, ...) 
{
    data <- as.numeric(data)
    myrank <- function(x, na.last = TRUE) {
        ranks <- sort.list(sort.list(x, na.last = na.last))
        if (is.na(na.last)) 
            x <- x[!is.na(x)]
        for (i in unique(x[duplicated(x)])) {
            which <- x == i & !is.na(x)
            ranks[which] <- max(ranks[which])
        }
        ranks
    }
    data <- sort(data)
    n.excess <- unique(floor(length(data) - myrank(data)))
    points <- unique(data)
    nl <- length(points)
    n.excess <- n.excess[-nl]
    points <- points[-nl]
    
    lower <- which(points > from)[1]
    
    
    excess <- cumsum(rev(data))[n.excess] - n.excess * points
    y <- excess/n.excess
    xx <- points[lower:(nl - omit)]
    yy <- y[lower:(nl - omit)]
    plot(xx, yy, xlab = "", ylab = "", ...)
    if (labels) 
        title(xlab = "Threshold", ylab = "Mean Excess")
    invisible(list(x = xx, y = yy))
}



## Configured VGAM::meplot to allow for omition of the last points
## and anything below a lower bound. 
meplot3 <- function(y, from = 0,omit = 3, labels = TRUE,
                    main = "Mean Excess Plot", 
                    xlab = "Threshold", ylab = "Mean Excess", lty = c(2,1,2),
                    conf = 0.95, col = c("blue", "black", "blue"), 
                    type = c("l","p","l"),pch = 1,...){
    if (!VGAM::is.Numeric(y)) 
        stop("bad input for argument 'y'")
    n <- length(y)
    sy <- sort(y)
    dsy <- rev(sy)
    me <- rev(cumsum(dsy))/(n:1) - sy
    me2 <- rev(cumsum(dsy^2))
    var <- (me2 - (n:1) * (me + sy)^2)/(n:1)
    ci <- qnorm((1 + conf)/2) * sqrt(abs(var))/sqrt(n:1)
    ci[length(ci)] <- NA
    
    lower <- which(sy > from)[1]
    upper <- n - omit
    bound <- lower + seq_len(max(upper-lower+1,0)) - 1
    
    mymat <- cbind(me - ci, me, me + ci)[bound,]
    sy <- (sy - sqrt(.Machine$double.eps))[bound]
    matplot(sy, mymat, main = main, xlab = xlab, ylab = ylab, 
            lty = lty, col = col, type = type, pch = pch,...)
    invisible(list(threshold = sy, meanExcess = me, plusminus = ci))
}

## ---- mainqqplot ----

qqplot2 <- function(distr, data,omit = 3, conf = 0.95,
                    namespace = NULL,log = FALSE, ...){
    ## distr the theoretical distribution quantile function
    ## fit parameters fitted (vector)
    ## data the empirical distribution
    ## ... parameter values for distribution
    
    
    qdistr <- get(paste0("q",distr))
    ddistr <- get(paste0("d",distr))
    
    thepoints <- ppoints(500)
    n <- length(thepoints)
        
    main <- paste("QQ-plot of data vs", distr)
    ylab <- deparse(substitute(data))
    xlab <- paste("Theoretical quantiles of",distr)
    
    qdistrpar <- function(p) qdistr(p,...)
        
    out <- qqplot4(qdistr(thepoints,...),data,
                   main = main,xlab = xlab,ylab = ylab,omit = omit)
    
    line <- qqline2(data,distribution = qdistrpar)
    ## Standard error of order statistics (asymptotic)
    fit.line <- line$int + line$slope * qdistr(thepoints, ...)
    
    se <- (line$slope/ddistr(qdistrpar(thepoints), ...))*
        sqrt(thepoints * (1 - thepoints)/n)
    
    err <- se * qnorm(1-(1-conf)/2) 
    upper <- fit.line + err
    lower <- fit.line - err
    
    lines(qdistrpar(thepoints), upper, lty = 2, col = 2)
    lines(qdistrpar(thepoints), lower, lty = 2, col = 2)

    invisible(list(out,lower = lower, upper = upper))
}

qqplot2("pareto",rpareto(1e4,3,5),shape = 3,scale = 5,omit = 1)

## ---- additional ----

qqplot4 <- function (x, y, plot.it = TRUE, xlab = deparse(substitute(x)), 
                     ylab = deparse(substitute(y)), omit = 3,...) {
    sx <- sort(x)
    sy <- sort(y)
    lenx <- length(sx)
    leny <- length(sy)
    if (leny < lenx) 
        sx <- approx(1L:lenx, sx, n = leny)$y
    if (leny > lenx) 
        sy <- approx(1L:leny, sy, n = lenx)$y
    sx <- sx[1:(length(sx)-omit)]
    sy <- sy[1:(length(sy)-omit)]
    
    if (plot.it) 
        plot(sx,sy,xlab = xlab, ylab = ylab, ...)
    invisible(list(x = sx, y = sy))
}



qqline2 <- function (y, datax = FALSE, distribution = qnorm, omit = 3,
                     probs = c(0.25, 0.75), qtype = 7, ...) {
    stopifnot(length(probs) == 2, is.function(distribution))
    y <- y[-tail(y,3)]
    y <- quantile(y, probs, names = FALSE, type = qtype, na.rm = TRUE)
    x <- distribution(probs)
    if (datax) {
        slope <- diff(x)/diff(y)
        int <- x[1] - slope * y[1]
    }
    else {
        slope <- diff(y)/diff(x)
        int <- y[1] - slope * x[1]
    }
    abline(int, slope, ...)
    invisible(list(int = int, slope = slope))
}



## ---- distributions ----
## Log pareto quantile
qlpareto <- function(p,shape, scale){
    sapply(p,function(x){
        log(qpareto(x,shape,scale))
    })
}

parloglike <- function(par,data){
    alpha <- par[1]
    beta <- par[2]
    n <- length(data)
    
    out <- n*log(alpha)+alpha*n*log(beta) - (alpha+1)*sum(log(beta+data))
    #     out <- sum(log(dpareto(udg,alpha,beta)))
    -out
}

lgloglike <- function(par,data){
    alpha <- par[1]
    beta <- par[2]
    n <- length(data)
    
    out <- n*alpha*log(beta) - n * lgamma(alpha) + 
        sum((alpha-1)*log(log(data))) - sum((beta+1)*log(data))
    -out
}



## Allow for 2*10^5 notation in plots
pretty10 <- function(x){
    digits <- 7
    eT <- floor(log10(abs(x)) + 10^-digits)
    mT <- signif(x/10^eT, digits)
    substitute(mT %*% 10^eT)
}




## ---- evir_myversion ----
hill2 <- function (data, option = c("alpha", "xi", "quantile"), start = 15, 
                   end = NA, reverse = FALSE, p = NA, ci = 0.95, auto.scale = TRUE, 
                   labels = TRUE, ...) 
{
    if(!"package:evir" %in% search()) library(evir)
    dots <- list(...)
    data <- as.numeric(data)
    ordered <- rev(sort(data))
    ordered <- ordered[ordered > 0]
    n <- length(ordered)
    option <- match.arg(option)
    if ((option == "quantile") && (is.na(p))) 
        stop("Input a value for the probability p")
    if ((option == "quantile") && (p < 1 - start/n)) {
        cat("Graph may look strange !! \n\n")
        cat(paste("Suggestion 1: Increase `p' above", 
                  format(signif(1 - start/n, 5)), "\n"))
        cat(paste("Suggestion 2: Increase `start' above ", 
                  ceiling(length(data) * (1 - p)), "\n"))
    }
    k <- 1:n
    loggs <- logb(ordered)
    avesumlog <- cumsum(loggs)/(1:n)
    xihat <- c(NA, (avesumlog - loggs)[2:n])
    alphahat <- 1/xihat
    y <- switch(option, alpha = alphahat, xi = xihat, quantile = ordered * 
                    ((n * (1 - p))/k)^(-1/alphahat))
    ses <- y/sqrt(k)
    if (is.na(end)) 
        end <- n
    x <- trunc(seq(from = min(end, length(data)), to = start))
    y <- y[x]
    x2 <- if("xlim" %in% names(dots)) x[x>=dots$xlim[1] & x<=dots$xlim[2]] else x
    y2 <- y[which(x %in% x2)]
    
    ylabel <- option
    yrange <- range(y)
    if (ci && (option != "quantile")) {
        qq <- qnorm(1 - (1 - ci)/2)
        u <- y2 + ses[x2] * qq
        l <- y2 - ses[x2] * qq
        ylabel <- paste(ylabel, " (CI, p =", ci, ")", sep = "")
        yrange <- range(u, l)
    }
    if (option == "quantile") 
        ylabel <- paste("Quantile, p =", p)

    index <- x2
    if (reverse) 
        index <- -index
    if (auto.scale) 
        plot(index, y, ylim = yrange, type = "l", xlab = "", 
             ylab = "", axes = FALSE, ...)
    else plot(index, y2, type = "l", xlab = "", ylab = "", axes = FALSE, 
              ...)
    axis(1, at = at <- pretty(index,50), labels = paste(pretty(x2,50)))
    axis(2)
    threshold <- rev(findthresh(data, x2))
    
    at3 <- seq(from = min(index), to = length(data), length.out = 20)
    
    axis(3, at = at3, labels = paste(format(signif(threshold[at3], 
                                                     3))))

    box()
    if (ci && (option != "quantile")) {
        lines(index, u, lty = 2, col = 4)
        lines(index, l, lty = 2, col = 4)
    }
    if (labels) {
        title(xlab = "Order Statistics", ylab = ylabel)
        mtext("Threshold", side = 3, line = 3)
    }
    invisible(list(x = index, y = y))
}



pickand <- function(data,k){
    datastat <- sort(data,TRUE)
    1/log(2)*log((datastat[k]-datastat[2*k])/
                     (datastat[2*k] - datastat[4*k]))
}
pickand <- Vectorize(pickand,"k")

dedh <- function(data,k){
    datastat <- sort(data,TRUE)
    h1 <- mean(log(datastat[1:k])-log(datastat[k+1]))
    h2 <- 1/k * sum((log(datastat[1:k])-log(datastat[k+1]))^2)
    1 + h1 + 1/2 * (h1^2/h2 - 1)^(-1)
}
dedh <- Vectorize(dedh,"k")




## ---- script_manipulation ----

## Source only a partition of a .R script
source_work <- function(file,textend,...){
    ## textend an exact match of a codeline that is the last 
    ## to be read
    ## ... optional arguments to readLines
    text <- readLines(file,...)
    end <- which(text == ifelse(missing(textend),"## ENDREAD",
                                textend))
    con <- textConnection(paste(text[1:end],collapse = "\n"))
    source(con)
}



