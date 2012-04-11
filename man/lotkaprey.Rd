\name{lotkaprey}
\alias{lotkaprey}
\alias{plot.lotkaprey}
\alias{traj.lotkaprey}
\alias{lines.lotkaprey}

\title{
  Lotka-Volterra Predator-Prey Model
}
\description{
  Function implements Lotka-Volterra predator-prey model.
}
\usage{
lotkaprey(a, b, c, d)
\method{plot}{lotkaprey}(x, xlim, ylim, arrows = 8, ...)
\method{lines}{lotkaprey}(x, N, P, time = 100, step = 0.2, ...)
\method{traj}{lotkaprey}(x, N, P, time = 100, step = 1, ...)
}

\arguments{
  \item{a, b, c, d}{Lotka-Volterra model parameters: see Details.}

  \item{x}{\code{"lotkaprey"} result object.}

  \item{xlim, ylim}{Axis ranges.}

  \item{arrows}{Number of trend arrows displayed in the graph.}

  \item{N, P}{Initial prey and predator densities}
  
  \item{time, step}{Total time span and time step for which trajectories
    are found.}

  \item{\dots}{Other arguments passed to functions.}

}

\details{
  The function uses the following parametrization for the prey \eqn{N}
  \deqn{\frac{dN}{dt} = aN - bNP \, }{dN/dt = a*N - b*N*P, }
  and for predator \eqn{P}
  \deqn{\frac{dP}{dt} = -cP + dbNP}{dP/dt = -c*P + d*b*N*P .}

  Function \code{lotkaprey} set ups the model, \code{plot} displays the
  phase diagram and isoclines, and \code{lines} adds trajectories. The
  trajectories can be generated with \code{traj}.

  Krebs (2009, p. 191) hardly mentions these models, although they are
  the basis of more advanced models.  The implementation is based on
  Stevens (2009).
}

\value{
  The function returns a \code{"lotkaprey"} result object which contains
  the input parameters and isocline values. Function \code{traj} returns
  trajectories with given starting values.
}

\author{
  Jari Oksanen
}
\note{
  The parametrization vary: the one in this function is the same as used
  in my lectures.
}

\references{
Krebs, C. J. (2009) \emph{Ecology.} Benjamin Cummings. 6 ed., 655 p.

Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R.} Springer. 
}


\seealso{
  \code{\link{lotkacomp}}, \code{\link{lotkapreyK}},
  \code{\link{rmprey}}. 
}
\examples{
mod <- lotkaprey(0.5, 0.01, 0.2, 0.1)
mod
op <- par(no.readonly = TRUE)
layout(matrix(c(1,1,2, 1,1,3), nrow=2, byrow = TRUE))
plot(mod)
lines(mod, 200, 20)
lines(mod, 200, 30)
plot(traj(mod, 200, 20))
plot(traj(mod, 200, 30))
par(op)
}
\keyword{ models }

