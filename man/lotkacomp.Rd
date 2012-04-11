\name{lotkacomp}
\alias{lotkacomp}
\alias{plot.lotkacomp}
\alias{traj.lotkacomp}
\alias{lines.lotkacomp}
\alias{summary.lotkacomp}

\title{
  Lotka-Volterra Competition Model
}
\description{
  Lotka-Volterra competition model for two species.
}
\usage{
lotkacomp(alpha, beta, K1, K2)
\method{plot}{lotkacomp}(x, arrows = 9, ...)
\method{lines}{lotkacomp}(x, N1 = 1, N2 = 1, r1 = 0.2, r2 = 0.2, ...)
\method{traj}{lotkacomp}(x, N1 = 1, N2 = 1, r1 = 0.2, r2 = 0.2,
    time = 100, step = 1, ...)
\method{summary}{lotkacomp}(object, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{alpha, beta}{
    Competition coefficients.
}
  \item{K1, K2}{
    Carrying capacities.
}
  \item{x, object}{
    A \code{lotkacomp} result object.
}
  \item{arrows}{
    Number of arrays drawn to show the area where species increases. 
  }
  \item{N1, N2}{Initial population sizes.}
  \item{r1, r2}{Population growth rates.}
  \item{time, step}{Total time and length of single time step in the
    trajectories.}
  \item{digits}{Number of digits in the printed output.}
  \item{\dots}{
    Other parameters passed to functions.
  }
}
\details{
  Function \code{lotkacomp} sets a Lotka-Volterra competition model
  (Krebs 2009, p. 165--169). The implementation is based on Stevens
  (2009). Function \code{plot} draws a phase graph. Function
  \code{lines} adds species trajectory to the phase graph, and
  \code{traj} extracts the population trajectories. Use
  \code{\link{plot.traj}} to display the population sizes against
  time. The \code{summary} gives the winner of the competition and the
  outcome population sizes.
}

\value{
  Function \code{lotkacomp} returns a Lotka-Volterra competition object.
}

\author{
  Jari Oksanen
}

\references{
Krebs, C. J. (2009) \emph{Ecology.} Benjamin Cummings. 6 ed., 655 p.

Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R.} Springer. 
}


\examples{
mod <- lotkacomp(0.8, 0.6, 20, 20)
plot(mod)
## Add trajectory line and starting point
lines(mod, 1, 1, col="green", lwd=2)
points(1, 1, pch=16, col="green")
## Plot populations against time
plot(traj(mod, 1, 1))
## Textual output
mod
summary(mod)
}
\keyword{ models }

