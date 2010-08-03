\name{do_admb}
\alias{do_admb}
\title{
  Compile and/or run an ADMB model, collect output
}
\description{
  Compile an ADMB model, run it, collect output
}
\usage{
do_admb(fn, data, params,
bounds=NULL,
phase=NULL,
re = FALSE,
re_vectors = NULL,
safe = TRUE,
mcmc = FALSE,
mcmc2 = FALSE,
mcmcsteps = 1000,
mcmcsave = round(mcmcsteps/1000),
mcmcpars,
impsamp = FALSE,
verbose=FALSE,
wd=getwd(),
  checkparam=c("stop","warn","write","ignore"),
  checkdata=c("stop","warn","write","ignore"),
objfunname="f",
clean=FALSE,
extra.args)
}
\arguments{
  \item{fn}{(character) base name of a TPL function,
    located in the working directory}
  \item{data}{
    a list of input data variables
    (order must match TPL file)}
  \item{params}{a list of starting parameter values
    (order must match TPL file)}
  \item{bounds}{named list of 2-element vectors of lower and upper
    bounds for specified parameters}
  \item{phase}{named numeric vector of phases (not implemented yet)}
  \item{re}{(logical) use random effects?}
  \item{re_vectors}{a named numeric vector of the identities and lengths
    of any random effects vectors used in the TPL file}
  \item{safe}{(logical) compile in safe mode?}
  \item{mcmc}{(logical) run MCMC around best fit?}
  \item{mcmc2}{(logical) run MCMC including random effects (RE only)?}
  \item{mcmcsteps}{integer) number of MCMC steps to run}
  \item{mcmcsave}{(integer) number of MCMC steps to save}
  \item{mcmcpars}{(character) vector of parameters for which to save
    MCMC output}
  \item{impsamp}{(logical) run importance sampling?}
  \item{verbose}{(logical) print details}
  \item{wd}{working directory}
  \item{checkparam}{result of parameter checks: "stop" to stop with
    an error when there is any mismatch, "warn" to warn,
    "write" to attach an appropriate section to the TPL file,
    and "ignore" to ignore mismatches}
  \item{checkdata}{result of data checks (as for checkparam)}
  \item{objfunname}{(character) name for objective function in TPL file
    (only relevant if \code{checkparam} is set to "write")}
  \item{clean}{remove output and working files when finished}
  \item{extra.args}{(character) extra argument string to pass to admb}
}
\value{
  An object of class \code{admb}.
}
\details{
  \code{do_admb} will attempt to do everything required to start from the
  model definition (TPL file) specified by \code{fn}, the data list, and
  the list of input parameters, compile and run (i.e. minimize the
  objective function of) the model in AD Model Builder, and read the
  results back into an object of class \code{admb} in R. If
  \code{checkparam} or \code{checkdata} are set to "write", it will
  attempt to construct a DATA section, and construct or
  (augment an existing) PARAMETER section (which may contain definitions
  of non-input parameters to be used in the model). It copies the input
  TPL file to a backup (.bak); on finishing, it restores the original TPL
  file and leaves the auto-generated TPL file in a file called [fn]_gen.tpl.
}
\note{1. Mixed-case file names are ignored by ADMB; this function makes
  a temporary copy with the file name translated to lower
  case. 2. Parameter names containing periods/full stops will not work,
  because this violates C syntax (currently not checked). 3. There are
  many, many, implicit restrictions and assumptions: for example,
  all vectors and matrices are assumed to be indexed starting from 1.
}
\author{
  Ben Bolker
}
\keyword{misc}
