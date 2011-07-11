\name{mcmc.control}
\alias{mcmc.control}
\title{
  Control options for MCMC after ADMB fitting
}
\description{
  Determines the options (number of steps, save interval, etc.)
  for running MCMC based on the estimated mode (maximum likelihood
  estimate) and parameter variance-covariance matrix
}
\usage{
mcmc.control(mcmc = 1000, mcmc2=0, mcsave, mcnoscale = FALSE, mcgrope =
FALSE, mcmult = 1, mcmcpars=NULL)
}
\arguments{
  \item{mcmc}{
    Total number of MCMC steps
  }
  \item{mcmc2}{
    MCMC2 steps (see ADMB-RE manual)
  }
  \item{mcsave}{
    Thinning interval.  Default is \code{pmax(1,floor(mcmc/1000))},
    i.e. aim to save 1000 steps 
  }
  \item{mcnoscale}{
    
  }
  \item{mcgrope}{
    
}
  \item{mcmult}{
    
  }
  \item{mcmcpars}{

  }
}
\details{
  See the AD Model Builder reference manual
}
\value{
  Returns a list of options suitable for passing as the \code{mcmc.opts}
  argument to \code{\link{do_admb}}
}
\author{
  Ben Bolker
}
\note{
  Some options (\code{mcmc2}, etc.) that can be used in AD Model Builder
  and ADMB-RE may not be available
}
\examples{
mcmc.control(mcmc=2000)
}
\keyword{misc}
