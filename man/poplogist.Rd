\name{poplogist}
\alias{poplogist}

\title{
  Logistic Population Growth
}
\description{
  Function returns the expected value at a given time for the logistic
  population growth.
}
\usage{
poplogist(t, r, N0, K)
}

\arguments{
  \item{t}{
    Time.
}
  \item{r}{
    Instantaneous population growth rate.
}
  \item{N0}{
    Population size at time \eqn{t = 0}.
}
  \item{K}{
    Carrying capacity.
}
}
\details{
  
  The logistic model is usually parametrized in ecology as:  

  \deqn{N = \frac{K}{1 + \frac{K-N_0}{N_0}\exp(-rt)}}{N = K/(1 + (K-N0)/N0 * exp(-r * t))}

  Function \code{\link{plogis}} defines a logistic model with slightly
  different parametrization. A similar model to \code{poplogist} can be
  defined as \code{K * plogis(-log(K/N0 - 1) + r * t)}. Function
  \code{\link{SSlogis}} can be used to fit logistic model in non-linear
  regression (\code{\link{nls}}, \code{\link[nlme]{lme}}) with yet
  another parametrization.  \code{\link{SSlogis}} uses terms
  \code{Asym}, \code{xmid} and \code{scal} which can be transformed to
  \code{poplogist} terms as \code{K = Asym}, \code{r = 1/scal} and
  \code{N0 = Asym/(exp(xmid/scal) + 1)}.

}
\value{
  Function returns a vector of populations sizes corresponding to input
  vector of \code{t}.
}

\author{
  Jari Oksanen
}

\examples{
x <- seq(0,100)
plot(x, poplogist(x, 0.1, 10, 100), xlab = "Time (t)",
  ylab="Population size (N)", type = "l")
abline(h=100, lty=2)
}

\keyword{ models }
