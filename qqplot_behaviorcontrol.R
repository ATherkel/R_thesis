


coords <- qqplot(qpareto(ppoints(500),fitpar$estimate[1],fitpar$estimate[2]),
                 udg)

x <- coords$x
y <- coords$y

which(x == max(x))
which(y == max(y))
max(x)
max(y)
str(x)
str(y)

plot(x,y)
plot(x[-500],y[-500])

max(x[-500])
max(a[-3])
which(2 == c(1,2,3,4,2,2))


plot(x,y)
plot(x[-500],y[-500])
plot(sort(rnorm(10)),1:10)
plot(sort(rnorm(100)),1:100)
plot(sort(rnorm(1000)),1:1000)
plot(sort(rexp(1000)),1:1000)
plot(sort(rexp(10000)),1:10000)
plot(sort(rpareto(10000,0.7,1)),1:10000)

## If both axis are monotone data, R will put first value in 
## bottomleft and last value in topright. This behavior 
## scales the axis appropriately and is thus both expected and 
## wanted behavior.




