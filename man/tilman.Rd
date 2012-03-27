\name{tilman}
\alias{tilman}
\alias{diffTilman}
\alias{traj.tilman}
\alias{plot.tilman}
\alias{tilmandiagr}
\alias{print.tilman}

\title{
  Tilman Resource Competition Model
}

\description{ 
  Tilman's competition model for essential resources. The
  model can in principle contain any number of resources or species.
}

\usage{
tilman(S, a, r, m, k, c, Rnames, Spnames)
\method{plot}{tilman}(x, R, N, S, time = 40, step = 0.2, lwd = 1, col, ...)
tilmandiagr(x, R, N, S, col = c(4,2,1), add = FALSE, xlim, ylim, 
   xlab, ylab, ...)
\method{traj}{tilman}(x, R, N, time = 40, step = 1, ...)
}

\arguments{
  \item{S, a}{Supply points and mobilation rates for resources.}
  
  \item{r, m}{Per capita rates of population increase and mortality for
    each species.}

  \item{k, c}{Half-saturation constants and consumption rates for each
    species and resource. For multi-resource models, these must be
    matrices with a row for each resource and column for each
    species. The matrix can also be given as a single vector elements
    ordered by species (column).}

  \item{Rnames, Spnames}{Names of resources and species. Defaults will
    be used if not given.}

  \item{x}{A \code{tilman} result object.}

  \item{R, N}{Initial resource levels and population sizes. Defaults
    to supply level \code{S} for resources, and 1 to population
    levels.}

  \item{add}{Add new resource trajectories to an existing diagram.}

  \item{time, step}{Total time and time step for trajectories}

  \item{lwd, col}{The line width and colours for trajectory lines.}

  \item{xlim, ylim}{Axis limits.}

  \item{xlab, ylab}{Axis labels.}

  \item{\dots}{Other parameters passed to functions.}
}

\details{

  The populations grow by consuming resources, and the population
  growths stop when resources are depleted to a critical level
  \eqn{R^*}{R*}. The differential equation for the change in the
  available amount of resource \eqn{i} and species \eqn{j} is given as

  \deqn{\frac{dR_i}{dt}=a_i(S_i-R_i) - \sum_j \left(c_{ij} \left(\frac{dN_j}{dt} + N_j m_j \right) \right)\;,}{dR[i]/dt = a[i]*(S[i]-R[i]) - SUM(j)(c[ij] * (dN[j]/dt + N[j]*m[j])),}

  where \eqn{R} and \eqn{S} are current resource level and its supply
  point, eqn{a} is the resource mobilization rate from \eqn{S} to
  \eqn{R}, \eqn{N} is the current population size, \eqn{m} its per
  capita mortality rate, and \eqn{c} the resource consumption
  rate. The corresponding change in the population size (that also
  appears in the equation above) is

  \deqn{\frac{dN_j}{dt} = \min_i N_j \left(\frac{r_j R_i}{R_i + k_{ij}} - m_j \right) \;, }{dN[j]/dt = min(i) N[j] * (r[j] * R[i] / (R[i] + k[ij]) - m[j]),}

  where \eqn{r} is the per capita growth rate of \eqn{N} and \eqn{k}
  is the half-saturation parameter of the Monod equation that
  describes the relationship between resource use and population size,
  and \eqn{min} selects the resource that gives the smallest response
  for the species.

  The population growth stops when the resource level for the critical
  resource is

  \deqn{R^*_i = \frac{k_{ij} m_j}{r_j - m_j}\;.}{R*[i] = k[ij]*m[j]/(r[j] - m[j]).}

  The corresponding maximum population size is

  \deqn{N^*_j = \frac{a_i (S_i - R^*_i)}{c_{ij} m_j}\;.}{N*[j] = a[i] * (S[i] - R*[i])/(c[ij] * m[j]).}

  With one resource, the winner will be the species with lowest
  \eqn{R^*}{R*} with population size \eqn{N^*}{N*} and all other
  species will have \eqn{N^* = 0}{N* = 0}. The outcome of multi-resource
  model is best inspected with simulations and graphics.

  Function \code{plot} shows resource levels and population
  trajectories against time. The critical resource levels
  (\eqn{R^*}{R*}) are also displayed if there is only one resource or
  one species. Function \code{tilmandiagr} makes the traditional
  Tilman graphics of resource trajectories in a two-resource plan
  together with \eqn{R^*}{R*} lines for each species.  The diagram
  will be only drawn for two-resource models, but the number of
  species is not limited. Function \code{traj} is used to calculate
  the resource and population trajectories.

  Krebs (2009) discusses the model on p. 168--169.  The current
  implementation is based on Tilman (1982).

}

\value{

  Function returns an object of class \code{"tilman"} with following
  items:
  \item{S, a, r, m, k, c}{Model input parameters.}
  \item{Rstar, Nstar}{Equilibrium resource levels and population sizes
    for each species and resource.}
  \item{nres, nsp}{Number of resources and species.}
  \item{call}{Function call.}

}

\references{
  
  Krebs, C. J. (2009) \emph{Ecology.} Benjaming Cummings. 6 ed., 655 p.

  Tilman, D. (1982) \emph{Resource competition and community
  structure.} Princeton UP, 296pp.

}
\author{
  Jari Oksanen
}

\examples{
## One resource
mod <- tilman(S = 20, a = 0.5, r = c(2,1), m = 0.5, k = c(5,0.5), c = 0.1)
mod
plot(mod)
## Two resources
mod2 <- update(mod, S=c(20,20), k=c(5,0.5, 1,2), c = c(0.15, 0.1, 0.1, 0.2))
tilmandiagr(mod2)
tilmandiagr(mod2, S = c(10,20), add = TRUE, time = 100)
tilmandiagr(mod2, S = c(20,10), add = TRUE, time = 100)
op <- par(mfrow=c(3,1), mar=c(4,4,0,1)+.1)
plot(mod2)
plot(mod2, S = c(10,20))
plot(mod2, S = c(20,10))
par(op)
}

\keyword{ models }

