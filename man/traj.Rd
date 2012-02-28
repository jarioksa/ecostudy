\name{traj}
\alias{traj}
\alias{plot.traj}

\title{
  Extract Population Trajectories
}

\description{
  Function \code{traj} is a generic method that extracts time and
  population sizzes from population models.  The \code{plot} function
  presents population size(s) against time.
}

\usage{
traj(x, ...)
\method{plot}{traj}(x, xlab = "Time (t)", ylab = "Population size", lwd = 2, 
    lty = 1, col =c(2,4), ...)
}

\arguments{
  \item{x}{
  A population model or a \code{"traj"} object for \code{plot}.
}
  \item{xlab, ylab}{Axis labels.}
  \item{lwd, lty, col}{Line width, type and colours.}	
  \item{\dots}{
  Other arguments passed to functions; \code{plot} passes these to 
  \code{\link{matplot}}.
}
}
\details{
  The population growth models use \code{traj} to extract data for their
  plots.  The \code{plot} function only shows the population sizes against
  time.  For trajectories in models, see \code{lines} function with the
  appropriate model.
}
\value{
  Function returns a matrix with at least two columns. The first is time,
  and the later are population sizes.
}

\author{
  Jari Oksanen
}

\examples{
mod <- lotkacomp(1.4, 0.9, 30, 30)
plot(traj(mod, 10, 2))
}
\keyword{ graphics }
\keyword{ models }
