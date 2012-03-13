### Random logistic growth following Pielou

`ranlogist` <-
  function(N0, b, d, K, bcontr=0.5)
{
    c <- (K-N0)/N0
    r <- b - d
    ## refuse to work with decreasing populations
    if (r <= 0)
        stop("'r' must be positive")
    t <- 0
    N <- N0
    ## Estimate time that is approximately needed to grow an
    ## error-free logistic model to 0.998*K
    if (r > 0)
        tend <- ceiling(-log(0.002/c)/r)
    time <- 0
    now <- 0
    simN <- N0
    size <- N0
    while (now < tend) {
        adj <- (b-d)*(size/K)
        bN <- b - bcontr*adj
        dN <- d + (1-bcontr)*adj
        events <- (bN+dN)*size
        step <- rexp(1,events)
        now <- now+step
        if (runif(1) <= bN/(bN+dN)) {
            size <- size+1
        }
        else {
            size <- size-1
        }
        simN <- c(simN,size)
        time <- c(time,now)
        if (size == 0) break
    }
    out <- list(N0 = N0, b = b, d = d, r = r, K = K,
                bcontr = bcontr, timend = tend, time = time,
                N = simN, call = match.call())
    class(out) <- "ranlogist"
    out
}

## Extract trajectories

`traj.ranlogist` <-
    function(x, ...)
{
    out <- cbind(time = x$time, N = x$N)
    class(out) <- "traj"
    out
}

## plot

`plot.ranlogist` <-
    function(x, ...)
{
    ## Plot the trajectory
    xy <- traj(x)
    xlim <- c(0, max(x$timend, xy[,1]))
    ylim <- c(0, max(x$K, xy[,2]))
    plot(xy, xlim = xlim, ylim = ylim, ...)
    ## Plot the expected curve
    z <- seq(0, x$timend, len=101)
    y <- poplogist(z, x$r, x$N0, x$K)
    lines(z, y, col=4, lwd=1, ...)
    abline(h = x$K, col = "gray")
}

## Only simulated curve

`lines.ranlogist` <-
    function(x, ...)
{
    xy <- traj(x)
    lines(xy, ...)
}

## print

`print.ranlogist` <-
    function(x, ...)
{
    cat("\nRandom Logistic Population Growth\n")
    cat("\nCall:", deparse(x$call), "\n\n")
    cat("End time", x$timend, "after", length(x$time), "events\n\n")
    invisible(x)
}
