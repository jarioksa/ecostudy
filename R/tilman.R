### Tilman Resource Competition model

`tilman` <-
    function()
{
}

## parameters for number of species [S] and resources [R]
## N[S]: population sizes
## R[R]: amount of resources
## r[S]: species growth rate
## k[R,S]: half-saturation constants
## m[S]: mortality
## S[R]: resource supply point
## a[R]: resource mobilation rate
## c[R,S]: consumption rate

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
