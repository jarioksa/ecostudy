### Lotka-Volterra competition models

### lotkacomp creates a Lotka-Volterra object

`lotkacomp` <-
    function(alpha, beta, K1, K2)
{
    if (missing(K1))
        stop("K1 needs a value")
    if (missing(K2))
        stop("K2 needs a value")
    out <- list(`alpha` = alpha, `beta` = beta, `K1` = K1, `K2` = K2,
                `sp1` = list("x" = K1, "y" = K1/alpha),
                `sp2` = list("x" = K2/beta, "y" = K2))
    class(out) <- "lotkacomp"
    out
}

## draw a Lotka-Volterra phase diagrama

`plot.lotkacomp` <-
    function(x, arrows = 9, ...)
{
    xlim = c(0, 1.1*max(x$sp1$x, x$sp2$x))
    ylim = c(0, 1.1*max(x$sp1$y, x$sp2$y))
    plot(c(0,0), xlim = xlim, ylim = ylim, axes = FALSE, type = "n",
         xlab = "Population size of species 1",
         ylab = "Population size of species 2", xaxs = "i", yaxs = "i", ...)
    box()
    segments(0, x$sp1$y, x$sp1$x, 0, col = 2, lwd = 3, ...)
    segments(0, x$sp2$y, x$sp2$x, 0, col = 4, lwd = 3, ...)
    axis(1, at = c(x$sp1$x, x$sp2$x),
         labels =c(expression(K[1]), expression(K[2]/beta)))
    axis(2, at = c(x$sp1$y, x$sp2$y),
         labels = c(expression(K[1]/alpha), expression(K[2])))
    if (arrows) {
        px <- ppoints(arrows) * x$sp2$x
        arrows(px, 0, px, x$sp2$y - x$beta * px, col="blue",
               length=0.1, lwd=0.5)
        py <- ppoints(arrows) * x$sp1$y
        arrows(0, py, x$alpha * (x$K1/x$alpha - py), py, col = "red",
               length=0.1, lwd=0.5)
    }
}


## trajectory uses 'primer::lvcomp2' and numerical integration

`traj.lotkacomp` <-
    function(x, N1 = 1, N2 = 1, r1 = 0.2, r2 = 0.2, time = 100, step = 1, ...)
{
    ## primer::lvcomp2 parametrization
    parms <- c(r1 = r1, r2 = r2, a11 = 1/x$K1, a22 = 1/x$K2,
               a12 = x$alpha/x$K1, a21 = x$beta/x$K2)
    initialN <- c(N1, N2)
    out <- ode(y = initialN, times = seq(from = 0, to = time, by = step),
               func = lvcomp2, parms = parms)
    class(out) <- c("traj", class(out))
    out
}


## add trajectory line to a plot

`lines.lotkacomp` <-
    function(x, N1 = 1, N2 = 1, r1 = 0.2, r2 = 0.2, ...)
{
    out <- traj(x, N1 = N1, N2 = N2, r1 = r1, r2 = r2, ...)
    out <- out[,-1]
    lines(out, ...)
}

