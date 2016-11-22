
## ---- functions ----

plot2 <- function(...){
    plot(...,axes = FALSE,mgp = c(2.5,1,0))
    axis(1, tck = -.055, labels = NA)
    axis(2, tck = -.055, labels = NA)
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



