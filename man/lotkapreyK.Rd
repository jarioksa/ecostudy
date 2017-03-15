\name{lotkapreyK}
\alias{lotkapreyK}
\alias{plot.lotkapreyK}
\alias{traj.lotkapreyK}
\alias{lines.lotkapreyK}

\title{
  Lotka-Volterra Predator-Prey Model with Density Dependent Prey
}
\description{
  Function implements Lotka-Volterra predator-prey model.
}
\usage{
lotkapreyK(a, b, c, d, K)
\method{plot}{lotkapreyK}(x, xlim, ylim, arrows = 5, ...)
\method{lines}{lotkapreyK}(x, N, P, time = 100, step = 0.2, ...)
\method{traj}{lotkapreyK}(x, N, P, time = 100, step = 1, ...)
}

\arguments{
  \item{a, b, c, d, K}{Lotka-Volterra model parameters: see Details.}

  \item{x}{\code{"lotkaprey"} result object.}

  \item{xlim, ylim}{Axis ranges.}

  \item{arrows}{Number of trend arrows displayed in the graph.}

  \item{N, P}{Initial prey and predator densities}
  
  \item{time, step}{Total time span and time step for which trajectories
    are found.}

  \item{\dots}{Other arguments passed to functions.}

}

\details{ The model is identical to \code{\link{lotkaprey}} except that
  the prey \eqn{N} is density-dependent following the logistic
  population growth model with carrying capacity \eqn{K}, prey growth
  rate \eqn{a}, attack rate \eqn{b} and predator density \eqn{P}:
  \deqn{\frac{dN}{dt} = aN(1 - N/K) - bNP}{dN/dt = a*N*(1-N/K) - b*N*P}
  Other parameters and the predator model are identical to the basic
  Lotka-Volterra model (\code{\link{lotkaprey}}). The
  Rosenzweig-MacArtur model uses the same prey model as this function,
  but applies type II functional response for the predator. The
  density-dependent Lotka-Volterra model is not commonly represented in
  ecology textbooks, but it is an educating step between basic
  Lotka-Volterra and Rosenzweig-MacArthur models. The predator isocline
  is identical to the basic model, but the prey isocline is a line
  reaching \eqn{P = 0} at carrying capacity \eqn{K}, and the intercept
  at \eqn{N} is equal to the prey isocline of the basic Lotka-Volterra
  model. The density dependent Lotka-Volterra model always has a stable
  equilibrium point.

  Function \code{lotkapreyK} set ups the model, \code{plot} displays the
  phase diagram and isoclines, and \code{lines} adds trajectories. The
  trajectories can be generated with \code{traj}.
}
\value{
  The function returns a \code{"lotkapreyK"} result object which contains
  the input parameters and isocline values. Function \code{traj} returns
  trajectories with given starting values.
}

\author{
  Jari Oksanen
}
\note{
  The parametrization vary: the one in this function is the same as used
  in my lectures for the basic Lotka-Volterra model.
}


\seealso{
  \code{\link{lotkacomp}}, \code{\link{lotkaprey}}, \code{\link{rmprey}}.
}
\examples{
mod <- lotkapreyK(0.5, 0.01, 0.2, 0.1, K = 1000)
mod
op <- par(no.readonly = TRUE)
layout(matrix(c(1,1, 1, 1, 2, 2), ncol=2, byrow = TRUE))
plot(mod)
lines(mod, 200, 20)
plot(traj(mod, 200, 20))
par(op)
}
\keyword{ models }

