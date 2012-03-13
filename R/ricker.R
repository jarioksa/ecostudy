### Ricker function

`ricker` <-
    function(N, R, K, b)
{
    i <- 1
    a <- ((R-1)^(1/b))/K
    over <- 0
    tmp <- N
    repeat { 
        N <- (N*R)/(1+(a*N)^b)
        tmp <- c(tmp,N)
        i <- i+1
        if (N > 0.99*K) over <- over+1
        if (over > 10) break
    }
    out <- list(R = R, K = K, b = b, a = a, N = tmp,
                call = match.call())
    class(out) <- "ricker"
    out
}

## brief print

`print.ricker`  <-
    function(x, ...)
{
    cat("\n")
    cat("Ricker model with parameters\n")
    cat("R = ", x$R, ", K = ", x$K, ", b = ", x$b, sep = "", "\n")
    cat("and", length(x$N), "population values\n\n")
    invisible(x)
}

## extract ricker() results

`traj.ricker` <-
    function(x, ...)
{
    N <- x$N
    t <- seq_along(N) - 1
    out <- cbind(t, N)
    class(out) <- "traj"
    out
}

`plot.ricker`  <-
    function(x, which = c(1, 2),  ...)
{
    xy <- traj(x)
    if (any(which == 1)) {
        plot(xy, ...)
        abline(h = x$K, col = "gray")
    }
    if (any(which == 2)) {
        xy <- xy[,2]
        z <- seq(from=0, to=max(xy)+0.5, len=100)
        plot(z, (z*x$R)/(1+(x$a*z)^x$b), type="l", xlab=expression(N[t]),
             ylab=expression(N[t+1]))
        lines(z, z, type="l", lty=2)
        for (i in seq_len(length(xy)-1)) {
            segments(xy[i], xy[i], xy[i], xy[i+1], col = "red")
            segments(xy[i], xy[i+1], xy[i+1], xy[i+1], col = "red")
        }
    }
}

