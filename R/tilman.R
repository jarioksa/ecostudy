### Tilman Resource Competition model


## parameters for number of species [S] and resources [R]
## N[S]: population sizes (starting sizes for populations: traj) 
## R[R]: amount of resources (staritng size for resource pool: traj)
## r[S]: species growth rate
## k[R,S]: half-saturation constants
## m[S]: mortality
## S[R]: resource supply point
## a[R]: resource mobilation rate
## c[R,S]: consumption rate

`tilman` <-
    function(S, a, r, m, k, c, Rnames, Spnames)
{
    ## Interpret number of resources and number of species from input
    ## parameter
    nres <- length(S)
    nsp <- length(r)
    ## Check or complement other model parameters
    if (length(a) < nres)
        a <- rep(a, nres)
    if (length(m) < nsp)
        m <- rep(m, nsp)
    ## Take care that c and k are matrices
    if (!is.matrix(k))
        k <- matrix(k, nrow = nres, ncol = nsp, byrow = TRUE)
    if (!is.matrix(c))
        c <- matrix(c, nrow = nres, ncol = nsp, byrow = TRUE)
    ## Names
    if (missing(Rnames))
        Rnames <- paste("R", seq_len(nres), sep = "")
    if (missing(Spnames))
        Spnames <- paste("Sp", seq_len(nsp), sep = "")
    names(S) <- names(a) <- rownames(k) <- rownames(c) <- Rnames
    names(r) <- names(m) <- colnames(k) <- colnames(c) <- Spnames 
    ## Model solutions
    Rstar <- sweep(sweep(k, 2, m, "*"), 2, r-m, "/")
    Nstar <- a * sweep(-Rstar, 1,  S, "+") / sweep(c, 2, m, "*")
    out <- list(S = S, a = a, r = r, m = m, k = k, c = c,
                Rstar = Rstar, Nstar = Nstar,
                nres = nres, nsp = nsp, call = match.call())
    class(out) <- "tilman"
    out
}

## print

`print.tilman`  <-
    function(x, ...)
{
    cat("\nTilman Competition Model for Essential Resources\n")
    cat("\nCall:", deparse(x$call), "\n\n")
    cat("Number of Resources:", x$nres, "\n")
    cat("Number of Species:  ", x$nsp, "\n\n")
    cat("Equilibrium Levels of Resources:\n")
    print(x$Rstar)
    cat("\nCorresponding Population Sizes:\n")
    print(x$Nstar)
    cat("\n")
    invisible(x)
}

## Differential equations

`diffTilman` <-
    function(t, y, p, ...)
{
    Ry <- y[1:p$nres]
    Ny <- y[-(1:p$nres)]
    up <- outer(Ry, p$r, "*")
    down <- sweep(p$k, 1, Ry, "+")
    dN.Ndt <- apply(sweep(up/down, 2, p$m, "-"), 2, min)
    dR.dt1 <- p$a * (p$S - Ry) 
    dR.dt2 <- rowSums(sweep(p$c, 2, Ny * (dN.Ndt + p$m), "*")) 
    dR.dt <- dR.dt1 - dR.dt2
    dN.dt <- Ny * dN.Ndt
    list(c(dR.dt, dN.dt))
}

## trajectories

## THINK THIS: should S (supply point) be configurable here instead of
## being a constant model parameter -- this will be important for
## multi-resource models where the supply point can define the end
## result.

`traj.tilman` <-
    function(x, R, N, time = 40, step = 1, ...)
{
    if (missing(R))
        R <- x$S
    if (is.null(names(R)))
        names(R) <- paste("R", seq_len(x$nres), sep = "")
    if (missing(N)) 
        N <- rep(1, x$nsp)
    names(N) <- names(x$r)
    initial <- c(R, N)
    time <- seq(from = 0, to = time, by = step)
    out <- ode(y = initial, time = time, func = diffTilman, parms = x)
    class(out) <- c("traj", class(out))
    out
}

## plot

`plot.tilman` <-
    function(x, R, N, kind = c("time", "resource"), time = 40, step = 0.2,
             lwd = 1, col, ...)
{
    kind <- match.arg(kind)
    if (kind == "resource")
        .NotYetUsed("kind", error = FALSE)
    ## don't know how to do "resource" plot for more than two
    ## resources
    if (x$nres > 2)
        kind <- "time"
    ## default col
    if (missing(col))
        col <- c(1,4,2,3,5:8)
    ## get traj
    tr <- traj(x, R = R, N = N, time = time, step = step, ...)
    ## Scale resource axes to the same range as population sizes
    rax <- seq_len(x$nres) + 1
    resran <- range(tr[, rax])
    popmax <- max(tr[, -c(1, rax)])
    resmul <- popmax/resran[2]
    tr[,rax] <- resmul * tr[,rax]
    plot(tr, col = col, lwd = lwd, ...)
    ## Add Resource legend to the right if there is space
    if ((mar4 <- par("mar")[4]) >= 2) {
        ## logarithmic y axis is trickier: this uses grDevices
        ## functions .axisPars and axisTicks
        if (par("ylog")) {
            usr <- log10(resran)
            p <- .axisPars(usr, log = TRUE)
            at <- axisTicks(usr, log = TRUE, c(p$axp, p$n))
            axis(4, at = resmul * at, labels = at)
        } else {
            at <- pretty(resran)
            axis(4, at = resmul*at, labels=at)
        }
        if (mar4 >= 3)
            mtext("Resources", side = 4, line = 2)
    }
    ## If there is only one resource, plot its isoclines
    if (x$nres == 1)
        abline(h = x$Rstar * resmul, col = col[-1], lty = 3)
    ## if there is only one speices, plot all resource isoclines
    else if (x$nsp == 1)
        abline(h = x$Rstar * resmul, col = col, lty=3)
    invisible(tr)
}
