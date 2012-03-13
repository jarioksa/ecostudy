### Rosenzweig-MacArthur model for predator-prey relationship
### following Hank Stevens's Primer

### Prey with density dependence, predator with type II functional
### response

`rmprey` <-
    function(b, K, w, D, e, s)
{
    alpha <- 1/K
    ## prey isocline as a 2nd degree polynomial b0+b1*x+b2*x^2
    b0 <- b/w*D
    b1 <- b/w*(1 - alpha*D)
    b2 <- -b/w*alpha
    preyfun <- function(x) b0 + (b1 + b2*x)*x
    ## solve preyfun(x) = y for x
    invpreyfun <- function(x) {
        b0 <- b0 - x
        disc <- b1*b1 - 4*b0*b2
        if (disc < 0)
            c(NA, NA)
        else
            c((-b1 + sqrt(disc))/2/b2, (-b1 - sqrt(disc))/2/b2)
    }
    ## prey shape parameters
    opt <- -b1/2/b2
    top <- preyfun(opt)
    lims <- pmax(c(0,0), invpreyfun(0))
    ## predator isocline
    prediso <- s*D/(e*w - s)
    ## out
    out <- list(preyfun = preyfun, invpreyfun = invpreyfun, preyopt = opt,
                preytop = top, preylimits = lims, prediso = prediso, b = b,
                K = K, D = D, w = w, e = e, s = s, call = match.call())
    class(out) <- "rmprey"
    out
}

## plot

`plot.rmprey` <-
    function(x, xlim, ylim, arrows = 5, ...)
{
    if (missing(xlim))
        xlim <- c(0, 1.2*max(x$preylimits, x$prediso))
    if (missing(ylim))
        ylim <- c(0, 2*x$preytop)
    plot(0,0, type = "n", xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
         xlab = "Prey Density", ylab = "Predator Density", ...)
    pfun <- x$preyfun
    curve(pfun, x$preylimits[1], x$preylimits[2], col = 4, lwd = 2, add = TRUE)
    abline(v = x$prediso, col = 2, lwd = 2)
    points(x$prediso, x$preyfun(x$prediso), pch = 16)
    if (arrows) {
        py <- ppoints(arrows) * x$preytop
        for (i in seq_len(arrows)) {
            pr <- x$invpreyfun(py[i])
            arrows(max(xlim[1], pr[1]), py[i], min(xlim[2], pr[2]), py[i],
                   col = 4, lwd = 0.5, len = 0.1)
        }
        px <- ppoints(arrows) * (xlim[2] - x$prediso) + x$prediso
        arrows(px, ylim[1], px, ylim[2], col = 2, lwd = 0.5, len=0.1)
    }
}

## population trajectories

traj.rmprey <-
    function(x, N, P, time = 100, step = 1, ...)
{
    parms <- c(b = x$b, e = x$e, s = x$s, w = x$w, D = x$D,
               alpha = 1/x$K)
    initialN <- c(N, P)
    time <- seq(from = 0, to = time, by = step)
    out <- ode(y = initialN, time = time, func = predpreyRM,
               parms = parms)
    class(out) <- c("traj", class(out))
    out
}

## lines

lines.rmprey <- function(x, N, P, time = 100, step = 0.2, ...)
{
    out <- traj(x, N, P, time = time, step = step, ...)
    lines(out[,-1], ...)
}
