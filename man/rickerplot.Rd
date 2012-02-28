\encoding{UTF-8}
\name{rickerplot}
\alias{rickerplot}
 
\title{Draws a Ricker Plot}

\description{

  Function finds a deterministic discrete population growth
  model (Ricker model) and draws side by side population size against
  time and a phase diagram of population size and population increase. 

}

\usage{
rickerplot(N, R, K, b)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{N}{Initial populatin size. }
  \item{R}{Population growth rate. }
  \item{K}{Carrying capacity. }
  \item{b}{Time lag parameter. }
}

\details{

  The population size \eqn{N} at time \eqn{t+1} is given as
  \deqn{N_{t+1} = R N_t / (1 + (aN)^b) \,}{N[t+1] = r*N[t]/(1+(a*N)^b),}
  where
  \deqn{a = ((R-1)^{1/b})/K.}
}

\value{
  The function draws the plots and returns nothing.
}

\examples{
par(mfrow=c(2,2))
rickerplot(10, 2, 50, 1)
rickerplot(10, 2, 50, 4)
par(mfrow=c(1,1))
}
\keyword{ models }
