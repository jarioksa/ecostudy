\name{lotkacomp}
\alias{lotkacomp}
\alias{plot.lotkacomp}
\alias{lines.lotkacomp}
\alias{lotkaplot}

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
lotkaplot(x, N1 = 1, N2 = 1, r1 = 0.2, r2 = 0.2, ...)
}

\arguments{
  \item{alpha, beta}{
    Competition coefficients.
}
  \item{K1, K2}{
    Carrying capacities.
}
  \item{x}{
    A \code{lotkacomp} result object.
}
  \item{arrows}{
    Number of trajectory arrays drawn.
  }
  \item{N1, N2}{Initial population sizes.}
  \item{r1, r2}{Population growth rates.}
  \item{\dots}{
    Other parameters passed to functions.
  }
}
\details{
  Function \code{lotkacomp} sets a Lotka-Volterra competition
  model. Function \code{plot} draws a phase graph. Function
  \code{lines} adds species trajectory to the phase graph, and
  \code{lotkaplot} shows population sizes against time.
}
\value{
  Function \code{lotkacomp} returns a Lotka-Volterra competition object.
}

\author{
  Jari Oksanen
}

\examples{
mod <- lotkacomp(0.8, 0.6, 20, 20)
plot(mod)
}
\keyword{ models }

