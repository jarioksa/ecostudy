\name{tilman}
\alias{tilman}
\alias{diffTilman}
\alias{traj.tilman}

\title{
  Tilman Resource Competition Model
}

\description{ 
  Tilman's competition model for essential resources. The
  model can in principle contain any number of resources or species.
}

\usage{
tilman(S, a, r, m, k, c, Rnames, Spnames)
\method{traj}{tilman}(x, R, N, time = 40, step = 1, ...)
}

\arguments{
  \item{S, a}{Supply points and mobilation rates for resources.}
  
  \item{r, m}{Per capita rates of population increas and mortality for
    each species.}

  \item{k, c}{Half-saturation constants and consumption rates for each
    species and resource. For multi-resource models, these must be
    matrices with a row for each resource and column for each
    species. The matrix can also be given as a single vector by row.}

  \item{Rnames, Spnames}{Names of resources and species. Defaults will
    be used if not given.}

  \item{x}{A \code{tilman} result object.}

  \item{R, N}{Initial resource levels and population sizes. Defaults
    to supply level \code{S} for resources, and 1 to population
    levels.}

  \item{time, step}{Total time and time step for trajectories}

  \item{\dots}{Other parameters passed to functions.}
}

\details{
  Details to come.
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
  
  Tilman, D. (1982) \emph{Resource competition and community
   structure.} Princeton UP, 296pp.

}
\author{
  Jari Oksanen
}

\examples{
mod <- tilman(S = 20, a = 0.5, r = c(2,1), m = 0.5, k = c(5,0.5), c = 0.1)
str(mod)
tr <- traj(mod)
tr[nrow(tr),]
}

\keyword{ models }

