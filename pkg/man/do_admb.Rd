\name{do_admb}
\alias{do_admb}
\title{
  Compile and/or run an ADMB model, collect output
}
\description{
  Compile an ADMB model, run it, collect output
}
\usage{
do_admb(fn, input_list, param_list,
re = FALSE, safe = TRUE, mcmc = FALSE, mcmc2 = FALSE, mcmcsteps = 1000,
mcmcsave = round(mcmcsteps/1000), impsamp = FALSE, verbose=FALSE,
wd=getwd(), clean=FALSE)
}
\arguments{
  \item{fn}{(character) base name of a TPL function,
    located in the working directory}
  \item{input_list}{
    a list of input variables
    (order must match TPL file)}
  \item{param_list}{a list of starting parameter values
    (order must match TPL file)}
  \item{re}{(logical) use random effects?}
  \item{safe}{(logical) compile in safe mode?}
  \item{mcmc}{(logical) run MCMC around best fit?}
  \item{mcmc2}{(logical) run MCMC including random effects (RE only)?}
  \item{mcmcsteps}{integer) number of MCMC steps to run}
  \item{mcmcsave}{(integer) number of MCMC steps to save}
  \item{impsamp}{(logical) run importance sampling?}
  \item{verbose}{(logical) print details}
  \item{wd}{working directory}
  \item{clean}{remove output and working files when finished}
}
\value{
  An object of class "admb"
}
\author{
  Ben Bolker
}
\keyword{misc}
