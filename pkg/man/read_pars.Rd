\name{read_pars}
\alias{read_pars}
\title{
  Read in parameters from an AD Model Builder run
}
\description{
  Reads coefficients, standard errors, log-likelihoods,
  maximum gradients, correlation and variance-covariance matrices
  from AD Model Builder output files
}
\usage{
read_pars(fn)
}
\arguments{
  \item{fn}{
    (character) Base name of AD Model Builder file
}
}
\details{
  Given the output from an ADMB run on FOO.tpl, reads the files
  FOO.par (parameters, log-likelihood, max gradient) and FOO.std (standard
  deviations, correlations)
}
\value{
  \item{coefficients}{parameter estimates}
  \item{se}{estimated standard errors of coefficients}
  \item{loglik}{log-likelihood}
  \item{maxgrad}{maximum gradient of log-likelihood surface}
  \item{cor}{correlation matrix}
  \item{vcov}{variance-covariance matrix}
}
\author{
  Ben Bolker
}
\seealso{
  \code{\link{write_pin}}, \code{\link{write_dat}}
}
\keyword{misc}
