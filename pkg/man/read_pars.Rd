\name{read_pars}
\alias{read_psv}
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
read_pars(fn,drop_phase=TRUE)
read_psv(fn,names)
}
\arguments{
  \item{fn}{(character) Base name of AD Model Builder }
  \item{names}{(character) Names of variables}
  \item{drop_phase}{(logical) drop negative-phase (fixed) parameters
    from results?}
}
\details{
  Given the output from an ADMB run on FOO.tpl, \code{read_pars}
  reads the files
  FOO.par (parameters, log-likelihood, max gradient); FOO.std (standard
  deviations); and FOO.cor (correlations).
  \code{read_psv} reads the output of MCMC runs
}
\value{
  \item{coefficients}{parameter estimates}
  \item{coeflist}{parameter estimates in list format, with proper shape
    (vectors, matrices, etc.)}
  \item{se}{estimated standard errors of coefficients}
  \item{loglik}{log-likelihood}
  \item{maxgrad}{maximum gradient of log-likelihood surface}
  \item{cor}{correlation matrix}
  \item{vcov}{variance-covariance matrix}
  \item{npar}{number of parameters}
}
\section{Warning}{
The \code{coeflist} component is untested for data structures more
complicated than scalars, vectors or matrices (i.e. higher-dimensional
or ragged arrays)}
\author{
  Ben Bolker
}
\seealso{
  \code{\link{write_pin}}, \code{\link{write_dat}}
}
\keyword{misc}
