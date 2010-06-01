\name{do_admb}
\alias{do_admb}
\title{
  Compile and/or run an ADMB model, collect output
}
\description{
  Compile an ADMB model, run it, collect output
}
\usage{
do_admb(fn, data_list, param_list,
re = FALSE,
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
clean=FALSE,
extra.args)
}
\arguments{
  \item{fn}{(character) base name of a TPL function,
    located in the working directory}
  \item{data_list}{
    a list of input data variables
    (order must match TPL file)}
  \item{param_list}{a list of starting parameter values
    (order must match TPL file)}
  \item{re}{(logical) use random effects?}
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
  \item{clean}{remove output and working files when finished}
  \item{extra.args}{(character) extra argument string to pass to admb}
}
\value{
  An object of class "admb"
}
\author{
  Ben Bolker
}
\keyword{misc}
