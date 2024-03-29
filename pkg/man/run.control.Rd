\name{run.control}
\alias{run.control}
\title{
  set run options for running ADMB via R
}
\description{
  A helper function 
}
\usage{
run.control(check_tpl = TRUE, write_files = TRUE, checkparam = c("stop", "warn", "write", "ignore"), checkdata = c("stop", "warn", "write", "ignore"), compile = TRUE, run = TRUE, read_files = TRUE, clean_files = "all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{check_tpl}{
    Check the specified TPL file for problems?
  }
  \item{write_files}{
    Write out data and initialization files?
  }
  \item{checkparam}{
    How to check \code{PARAMETERS} section of the TPL file: \code{stop}=stop if
    there are problems; \code{warn}=give a warning if there are problems,
    but try to proceed; \code{write}=modify TPL file, writing appropriate
    sections; \code{ignore}=assume TPL file is OK, proceed
  }
  \item{checkdata}{
    as with \code{checkparam}: how/whether to check/generate the
    \code{DATA} section of the TPL file
  }
  \item{compile}{
    compile the TPL file (via ADMB) into an executable?
  }
  \item{run}{
    run the executable file with the specified data/initial values?
  }
  \item{read_files}{
    read the results of an ADMB run into R?
  }
  \item{clean_files}{
    Delete working files after completion of the run?
    Options are \code{"all"}, \code{"sys"}, \code{"output"},
    \code{"none"}; \code{TRUE} is equivalent to \code{"all"}
    and \code{FALSE} is equivalent to \code{"none"}
      
  }
}
\value{
  A list with appropriate default values inserted
  for passing to \code{\link{do_admb}}
}
\author{
  Ben Bolker
}
\keyword{misc}

