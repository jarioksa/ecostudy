### Milne Stochastic Ceiling Model

### Population grows with Normally distributed rate N(r, sd) but is
### bound to upper ceiling K. Populations with abundance < 0.5 are
### regareded as extinct.

`ranmilne` <-
    function(N0, r, sd, K, time = 100, rep = 1)
{
    k <- log(K)
    n <- log(N0)
    dead <- log(0.5)
    out <- matrix(-Inf, nrow = time + 1, ncol = rep)
    out[1,] <- n
    ## loop over replicates and timesteps
    for (j in seq_len(rep)) {
        for (i in 1:time) {
            out[i+1, j] <- out[i, j] + rnorm(1, r, sd)
            if (out[i+1, j] < dead)
                break
            if (out[i+1, j] > k)
                out[i+1, j] <- k
        }
    }
    out <- exp(out)
    out <- cbind("Time" = seq_len(nrow(out)), out)
    class(out) <- "traj"
    obj <- list("traj" = out, "N0" = N0, "r" = r, "sd" = sd,
                "K" = K, "time" = time, "rep" = rep,
                "call" = match.call())
    class(obj) <- c("ranmilne")
    obj
}

## traj extracts simulation

`traj.ranmilne` <-
    function(x, ...)
{
    x$traj
}

## plot() and lines() use *.traj()

`plot.ranmilne`  <-
    function(x, lty = 1, lwd = 1, col = 4,  ...)
{
    plot(traj(x), lty = lty, lwd = lwd, col = col, ...)
}

`lines.ranmilne` <-
    function(x, ...)
{
    lines(traj(x), ...)
}

## print

`print.ranmilne` <-
    function(x, ...)
{
    cat("\nMilne Stochastic Ceilling Population Model\n")
    cat("\nCall:", deparse(x$call), "\n\n")
    cat("Model parameters: r =", x$r, " sd =", x$sd, " K =", x$K,"\n")
    cat("Initial Population: N0 =", x$N0, "\n")
    cat("Time steps", x$time, "  no. of replicate runs", x$rep, "\n\n")
    invisible(x)
}
