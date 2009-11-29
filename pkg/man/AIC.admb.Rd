\name{AIC.admb}
\alias{AIC.admb}
\alias{vcov.admb}
\alias{logLik.admb}
\alias{coef.admb}
\alias{deviance.admb}
\title{
  Standard accessor functions for ADMB model fits
}
\description{
  Extract standard information such as log-likelihood, AIC,
  coefficients, etc. from ADMB model fits
}
\usage{
\method{AIC}{admb}(object, \dots, k = 2)
\method{logLik}{admb}(object, \dots)
\method{deviance}{admb}(object, \dots)
\method{vcov}{admb}(object, \dots)
\method{coef}{admb}(object,\dots)
}
\arguments{
  \item{object}{an ADMB model fit (of class "admb")}
  \item{k}{penalty value for AIC fits}
  \item{\dots}{other parameters (for S3 generic compatibility)}
}
\value{
  Extracts appropriate values: numeric (scalar) for AIC, type logLik for logLik,
  numeric vector of coefficients, numeric variance-covariance matrix
  of parameter estimates
}
\author{
  Ben Bolker
}
\examples{
  admbex <- system.file("doc","ADMB_runs.RData",package="R2admb")
  load(admbex)
  m1
  coef(m1)
  summary(m1)
  coef(summary(m1)) ## returns just z-table
  AIC(m1)
  vcov(m1)
  logLik(m1)
  deviance(m1)
}
\keyword{misc}
