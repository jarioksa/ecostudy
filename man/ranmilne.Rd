\name{ranmilne}
\alias{ranmilne}
\alias{traj.ranmilne}
\alias{plot.ranmilne}
\alias{lines.ranmilne}

\title{
  Milne Stochastic Ceiling Model for Population Size
}
\description{
  
  Milne Ceiling model is a stochastic model where the population
  growth is in the range zero to carrying capacity.  It seems to be
  particularly suitable for small populations which are often close to
  their carrying capacity, and is much used in metapopulation research
  (Hanski 2005).

}
\usage{
ranmilne(N0, r, sd, K, time = 100, rep = 1)
\method{plot}{ranmilne}(x, lty = 1, lwd = 1, col = 4, ...)
\method{lines}{ranmilne}(x, ...)
\method{traj}{ranmilne}(x, ...)
}

\arguments{
  \item{N0}{
    Initial population size(s).
}
  \item{r, sd}{
    Normally distributed growth rate and its standard deviation
}
  \item{K}{
    Carrying capacity.
}
  \item{time}{
    Number of time steps.
}
  \item{rep}{
    Number of replicated simulations.
  }
  \item{x}{
    A \code{ranmilne} result object.
  }
  \item{lty, lwd, col}{
    Line type, width and colour.
  }
  \item{\dots}{
    Other arguments passed to functions.
  }
}

\details{

  The logarithmic population size \eqn{n} at time \eqn{t} is
  given as \eqn{n_{t-1} + N(r, s)} where \eqn{N()} is a random normal
  variate with expectation \eqn{r} and standard deviation \eqn{s}.
  Simulated responses \eqn{> K} are made to \eqn{K} and responses \eqn{<
  0.5} are made to \eqn{0} (extinct).

}

\value{
  The function returns an object of class \code{"ranmilne"} which is a
  list with following items:
  \item{traj}{Simulated results as a \code{\link{traj}} object. Function
    \code{traj} extracts the simulation results.}
  \item{N0, r, sd, K, time, rep}{Model parameters.}
  \item{call}{Function call.}
}

\author{
  Jari Oksanen
}
\note{
  The current version returns the simulated population size as a
  real number. It may be that the return value should be an integer, but
  this was not  yet checked.
}

\references{

  Hanski, I (2005) \emph{The Shrinking World: Ecological Consequences
  of Habitat Loss.} International Ecological Institute, 307 p.

}

\seealso{
  \code{\link{ranlogist}}.
}
\examples{
sim <- ranmilne(50, r = 0.5, sd = 1, K = 50, rep = 10)
sim
plot(sim)
}
\keyword{ models }

