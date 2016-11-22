library(microbenchmark)


mymax1 <- function(a){
    sapply(1:length(a),function(i){
        max(a[1:i])
    })
}

mymax2 <- function(myseq){
    m <- myseq[1]
    for(i in 2:length(myseq)){
        if((newm <- myseq[i]) > (last <- tail(m,1))){
            m <- c(m,newm)
        } else {
            m <- c(m,last)
        }
    }
    m
}



n <- 1e3
a <- rcauchy(n)
b <- rexp(n)



amax1 <- mymax1(a)
amax2 <- mymax2(a)


microbenchmark(mymax1(a))
microbenchmark(mymax2(a))


test <- function(i){
    n <- 1e4
    a <- rcauchy(n)
    b <- rexp(n)
    
    amax <- mymax2(a)
    bmax <- mymax2(b)
    # par(mfrow = c(2,1))
    # plot(amax,type = "l")
    # plot(bmax,type = "l")
    
    out <- c(length(unique(amax)),length(unique(bmax)))
    out
}

out <- sapply(1:100,test)
mean(out[1,]);var(out[1,])
mean(out[2,]);var(out[2,])

