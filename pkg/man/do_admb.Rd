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
re = NULL,
data_type = NULL,
safe = TRUE,
profile = FALSE,
profpars = NULL,
mcmc = FALSE,
mcmc.opts = mcmc.control(),
impsamp = FALSE,
verbose=FALSE,
run.opts=run.control(),
objfunname="f",
workdir=getwd(),
admb_errors=c("stop","warn","ignore"),
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
  \item{re}{a named list of the identities and dimensions
    of any random effects vectors or matrices used in the TPL file}
  \item{data_type}{a named vector specifying (optional)
    data types for parameters, in parname="storage mode" format
    (e.g. \code{c(x="integer",y="numeric")})}
  \item{safe}{(logical) compile in safe mode?}
  \item{profile}{(logical) generate likelihood profiles? (untested!)}
  \item{profpars}{(character) vector of names of parameters to profile}
  \item{mcmc}{(logical) run MCMC around best fit?}
  \item{mcmc.opts}{options for MCMC (see
    \code{\link{mcmc.control}} for details)}
  \item{impsamp}{(logical) run importance sampling?}
  \item{verbose}{(logical) print details}
  \item{run.opts}{options for ADMB run (see
    \code{\link{run.control}} for details)
  }
  \item{objfunname}{(character) name for objective function in TPL file
    (only relevant if \code{checkparam} is set to "write")}
  \item{workdir}{temporary working directory (dat/pin/tpl files will be copied)}
  \item{admb_errors}{how to treat ADMB errors (in either compilation or
    run): use at your own risk!}
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
\examples{
\dontrun{
setup_admb()
file.copy(system.file("tplfiles","ReedfrogSizepred0.tpl",package="R2admb"),"tadpole.tpl")
 tadpoledat <-
  data.frame(TBL = rep(c(9,12,21,25,37),each=3),
             Kill = c(0,2,1,3,4,5,0,0,0,0,1,0,0,0,0L),
             nexposed=rep(10,15))
m1 <- do_admb("tadpole",
              data=c(list(nobs=15),tadpoledat),
              params=list(c=0.45,d=13,g=1),
              bounds=list(c=c(0,1),d=c(0,50),g=c(-1,25)),
              run.opts=run.control(checkparam="write",
                checkdata="write",clean="all"))
unlink("tadpole.tpl")
}
}
\keyword{misc}
