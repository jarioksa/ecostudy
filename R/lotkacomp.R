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
        arrows(px, 0, px, x$sp2$y - x$beta * px, col="blue")
        py <- ppoints(arrows) * x$sp1$y
        arrows(0, py, x$alpha * (x$K1/x$alpha - py), py, col = "red")
    }
}