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
    function(S, a, r, m, k, c)
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
    ## Model solutions
    Rstar <- sweep(sweep(k, 2, m, "*"), 2, r-m, "/")
    Nstar <- a * sweep(-Rstar, 1,  S, "+") / sweep(c, 2, m, "*")
    list(Rstar = Rstar, Nstar = Nstar)
}

## Differential equations

odeTilman <-
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
