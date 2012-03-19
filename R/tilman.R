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
        m <- rep(a, nsp)
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
                nres = nres, nsp = nsp)
    class(out) <- "tilman"
    out
}

## Differential equations

`diffTilman` <-
    function(t, y, p, ...)
{
    Ry <- y[1:p$nres]
    Ny <- y[-(1:p$nres)]
    up <- outer(Ry, p$r, "*")
    down <- sweep(p$k, 1, Ry, "+")
    dN.Ndt <- up/down - p$m
    dR.dt1 <- p$a * (p$S - Ry) 
    dR.dt2 <- rowSums(sweep(p$c * sweep(dN.Ndt, 2, p$m, "+"), 2, Ny, "*"))
    dR.dt <- dR.dt1 - dR.dt2
    dN.dt <- apply(Ny * dN.Ndt, 2, min)
    list(c(dR.dt, dN.dt))
}

## trajectories

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
