\name{lynxhare}
\alias{lynxhare}
\docType{data}
\title{
Hudson Bay Trading Compary Lynx and Hare Trapping Data
}
\description{
  The number of pelts (in thousands) traded with the Hudson Bay Trading
  Company in 1845 to 1935. 
}
\usage{data(lynxhare)}
\format{
  A data frame with 91 observations on the following 3 variables.
  \describe{
    \item{\code{Year}}{Trading year.}
    \item{\code{Hare}}{Number of snowshoe hare pelts (thousands).}
    \item{\code{Lynx}}{Number of lynx pelts (thousands).}
  }
}
\details{
  These data are regularly used to demonstrate predator--prey
  cycles. However, the species are a part of more complex ecosystem
  (including plant--hare) and collected from various areas, and are
  hardly a demonstration of a predator--prey cycle alone.

  The lynx data are also provided in the standard \R dataset
  \code{\link{lynx}} (with slightly different time range). The dataset
  is a data frame that also inherits from class \code{"traj"} so that
  \code{\link{plot.traj}} can be used to display its results.
}

\source{
  Stevens, M.H.H. (2009). \emph{A Primer of Ecology with R.} Springer.
}

\examples{
data(lynxhare)
## plot uses traj() method
plot(lynxhare)
}
\keyword{datasets}
