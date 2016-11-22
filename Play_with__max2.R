n <- 1e5


a <- rpareto(n,1,1)
hist(a)
b <- sapply(1:10, function(i) rpareto(n,1,1))
# c <- rpareto(n,1,1)

bc <- cbind(b)

bcmax <- sapply(1:nrow(bc),function(i) max(bc[i,]))
mean(a>bcmax)
