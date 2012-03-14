### Lotka-Volterra Predator-Prey Model with Density Dependent Prey
### (Logistic Prey Growth).

### Set up the model

`lotkapreyK`  <-
    function(a, b, c, d, K)
{
    out <- list("a" = a, "b" = b, "c" = c, "d" = d, "K" = K,
                "Nslope" = -a/b/K, 
                "P0" = c/d/b, "N0" = a/b, call = match.call())
    out$sol <- out$N0 + out$Nslope * out$P0
    class(out) <- "lotkapreyK"
    out
}

`plot.lotkapreyK`  <-
    function(x, xlim, ylim, arrows = 5, ...)
{
    if (missing(ylim))
        ylim <- c(0, 2 * x$N0)
    if (missing(xlim))
        xlim <- c(0, 1.05 * x$K)
    plot(x$P0, x$N0, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
         xlab = "Prey Density", ylab = "Predator Density", type = "n", ...)
    abline(x$N0, x$Nslope, col = "blue", lwd = 2)
    abline(v = x$P0, col = "red", lwd = 2)
    ## arrows
    if (arrows) {
        ## prey
        py <- ppoints(arrows) * x$N0
        px <- (py - x$N0)/x$Nslope
        arrows(0, py, px, py, code = 2, length = 0.1,
               lwd = 0.5, col = "blue")
        ## predator
        px <- ppoints(arrows) * (xlim[2] - x$P0) + x$P0
        arrows(px, ylim[1], px, ylim[2],  code = 2, length = 0.1,
               lwd = 0.5, col = "red")
    }
}

### trajectory

`traj.lotkapreyK` <-
    function(x, N, P, time = 100, step = 1, ...)
{
    ## Function for ode
    odefun <- function(t, y, p) {
        N <- y[1]
        P <- y[2]
        with(as.list(p), {
            dN.dt = a * N * (1 - N/K) - b * N * P
            dP.dt = -c * P + d * b * N * P
            list(c(dN.dt, dP.dt))
        })
    }
    parms <- c(a = x$a, b = x$b, d = x$d, c = x$c, K = x$K)
    initialN <- c(N, P)
    time <- seq(from = 0, to = time, by = step)
    out <- ode(y = initialN, times = time, func = odefun,
               parms = parms)
    class(out) <- c("traj", class(out))
    out
}

## Add trajectory line to the phase plot

`lines.lotkapreyK` <-
    function(x, N, P, time = 100, step = 0.2, ...)
{
    out <- traj(x, N, P, time = time, step = step, ...)
    lines(out[,-1], ...)
}

`print.lotkapreyK` <-
    function(x, ...)
{
    cat("\n")
    cat("Lotka-Volterra Predator-Prey Model\n")
    cat("\nCall:", deparse(x$call), "\n\n")
    cat("Equilibrium at P =", x$P0, " N =" , x$sol, "\n\n")
    invisible(x)
}
