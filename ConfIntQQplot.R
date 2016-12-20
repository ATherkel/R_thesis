qqplot5 <- function (x, distribution = "norm", ..., ylab = deparse(substitute(x)), 
          xlab = paste(distribution, "quantiles"), main = NULL, las = par("las"), 
          envelope = 0.95, col = palette()[1], col.lines = palette()[2], 
          lwd = 2, pch = 1, cex = par("cex"), line = c("quartiles", 
                                                       "robust", "none"), 
          labels = if (!is.null(names(x))) names(x) else seq(along = x), 
          id.method = "y", id.n = if (id.method[1] == "identify") Inf else 0, 
          id.cex = 1, id.col = palette()[1], id.location = "lr", grid = TRUE) 
{
    line <- match.arg(line)
    good <- !is.na(x)
    ord <- order(x[good])
    if (length(col) == length(x)) 
        col <- col[good][ord]
    if (length(pch) == length(x)) 
        pch <- pch[good][ord]
    if (length(cex) == length(x)) 
        cex <- cex[good][ord]
    ord.x <- x[good][ord]
    ord.lab <- labels[good][ord]
    q.function <- eval(parse(text = paste("q", distribution, 
                                          sep = "")))
    d.function <- eval(parse(text = paste("d", distribution, 
                                          sep = "")))
    n <- length(ord.x)
    P <- ppoints(n)
    z <- q.function(P, ...)
    plot(z, ord.x, type = "n", xlab = xlab, ylab = ylab, main = main, 
         las = las)
    if (grid) {
        grid(lty = 1, equilogs = FALSE)
        box()
    }
    points(z, ord.x, col = col, pch = pch, cex = cex)
    if (line == "quartiles" || line == "none") {
        Q.x <- quantile(ord.x, c(0.25, 0.75))
        Q.z <- q.function(c(0.25, 0.75), ...)
        b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
        a <- Q.x[1] - b * Q.z[1]
        if (line == "quartiles") 
            abline(a, b, col = col.lines, lwd = lwd)
    }
    if (line == "robust") {
        coef <- coef(rlm(ord.x ~ z))
        a <- coef[1]
        b <- coef[2]
        abline(a, b, col = col.lines, lwd = lwd)
    }
    conf <- if (envelope == FALSE) 
        0.95
    else envelope
    zz <- qnorm(1 - (1 - conf)/2)
    SE <- (b/d.function(z, ...)) * sqrt(P * (1 - P)/n)
    return(SE)
    fit.value <- a + b * z
    upper <- fit.value + zz * SE
    lower <- fit.value - zz * SE
    if (envelope != FALSE) {
        lines(z, upper, lty = 2, lwd = lwd, col = col.lines)
        lines(z, lower, lty = 2, lwd = lwd, col = col.lines)
    }
    showLabels(z, ord.x, labels = ord.lab, id.method = id.method, 
               id.n = id.n, id.cex = id.cex, id.col = id.col, id.location = id.location)
}