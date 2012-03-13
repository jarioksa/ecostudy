### Lotka-Volterra Predator-Prey Model

### Set up Lotka-Volterra Predator-Prey Model

`lotkaprey`  <-
    function(a, b, c, d)
{
    out <- list("a" = a, "b" = b, "c" = c, "d" = d, 
                "P0" = c/d/b, "N0" = a/b, call = match.call())
    class(out) <- "lotkaprey"
    out
}

`plot.lotkaprey`  <-
    function(x, xlim, ylim, arrows = 8, ...)
{
    if (missing(ylim))
        ylim <- c(0, 3 * x$N0)
    if (missing(xlim))
        xlim <- c(0, 3 * x$P0)
    plot(x$P0, x$N0, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
         xlab = "Prey Density", ylab = "Predator Density", type = "n", ...)
    abline(h = x$N0, col = "blue", lwd = 2)
    abline(v = x$P0, col = "red", lwd = 2)
    ## arrows
    if (arrows) {
        ## prey
        py <- ppoints(arrows) * ylim[2]
        inc <- py < x$N0
        dec <- py > x$N0
        arrows(0, py[inc], x$P0, py[inc], code = 2, length = 0.1,
               lwd = 0.5, col = "blue")
        arrows(x$P0, py[inc], xlim[2], py[inc], code = 2, 
               length = 0.1, lwd = 0.5, col = "blue")
        arrows(0, py[dec], x$P0, py[dec], code = 1, length = 0.1,
               lwd = 0.5, col = "blue")
        arrows(x$P0, py[dec], xlim[2], py[dec], code = 1, 
               length = 0.1, lwd = 0.5, col = "blue")
        ## predator
        yup <- ylim[2]
        px <- ppoints(arrows) * xlim[2]
        inc <- px > x$P0
        dec <- px < x$P0
        arrows(px[dec], yup, px[dec], x$N0,  code = 2, length = 0.1,
               lwd = 0.5, col = "red")
        arrows(px[dec], x$N0, px[dec], 0,  code = 2, length = 0.1,
               lwd = 0.5, col = "red")
        arrows(px[inc], yup, px[inc], x$N0,  code = 1, length = 0.1,
               lwd = 0.5, col = "red")
        arrows(px[inc], x$N0, px[inc], 0,  code = 1, length = 0.1,
               lwd = 0.5, col = "red") 
    }
}

### trajectory

`traj.lotkaprey` <-
    function(x, N, P, time = 100, step = 1, ...)
{
    ## Based on primer package: translate parameter names
    parms <- c(b = x$a, a = x$b, e = x$d, s = x$c)
    initialN <- c(N, P)
    time <- seq(from = 0, to = time, by = step)
    out <- ode(y = initialN, times = time, func = predpreyLV,
               parms = parms)
    class(out) <- c("traj", class(out))
    out
}

## Add trajectory line to the phase plot

`lines.lotkaprey` <-
    function(x, N, P, time = 100, step = 0.2, ...)
{
    out <- traj(x, N, P, time = time, step = step, ...)
    lines(out[,-1], ...)
}

`print.lotkaprey` <-
    function(x, ...)
{
    cat("\n")
    cat("Lotka-Volterra Predator-Prey Model\n")
    cat("\nCall:", deparse(x$call), "\n\n")
    cat("Predator isocline: N =", x$P0, "\n")
    cat("Prey isocline:     P =", x$N0, "\n\n")
    invisible(x)
}
