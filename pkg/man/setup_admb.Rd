\name{setup_admb}
\alias{setup_admb}
\alias{clean_admb}
\title{
Set up AD Model Builder environment variables
}
\description{
Attempts to set environment variables so that AD Model Builder will "just work" when run from inside R
}
\usage{
setup_admb(admb_home)
clean_admb(fn, which=c("all","sys","output","none"),profpars=NULL)
}
\arguments{
  \item{admb_home}{(character) directory containing AD Model Builder
    binary files}
  \item{fn}{(character) base name of ADMB model files}
  \item{which}{what to remove (system/output files or input files)}
  \item{profpars}{(character) names of parameters that have been profiled (i.e.,
   \code{.plt} files to remove)}
}
\note{Only 'all' option is currently implemented in clean_admb()
  [i.e., remove all except TPL file]
  }
\details{
(1) If the environment variable ADMB_HOME is not already set and admb_home is not
specified, this function will try to set it sensibly. (I.e., on Unix systems, it will run a "locate" command (if one is available) to try
to find the binaries, and thereafter check if they are installed
in the default location (/usr/local/admb); on Windows it will assume they are installed in the default
location (C:/ADMB).) (2) If ADMB_HOME is set and admb_home is not specified,
it will leave the original setting alone. (3) If admb_home is specified, it
will set the environment variable ADMB_HOME to that value.

The function also appends
the admb_home value to the PATH variable.
}
\value{
A character vector containing the name of the current ADMB_HOME.
}
\author{
  Ben Bolker
}
\examples{
  orig <- Sys.getenv("ADMB_HOME")
  setup_admb("elsewhere")   ## this doesn't make sense but won't break anything until you actually try to run AD Model Builder
  Sys.setenv(ADMB_HOME="") ## erase environment variable
\dontrun{
  setup_admb()              ## auto-locate (fails if ADMB not found)
}
  Sys.setenv(ADMB_HOME=orig) ## restore sanity
}
\keyword{misc}
\keyword{environment}
