setup_admb <- function(admb_home) {
  ## check whether already set up
  sys_home <- Sys.getenv("ADMB_HOME")
  if (missing(admb_home)) { ## not passed to function
    if (nchar(sys_home)>0) {
    ## previously defined
      admb_home <- sys_home
    } else {
      ## not defined & not passed: try to guess
      admb_home <- ""
      if (.Platform$OS.type=="unix") {
        ## unix.  Should (1) check that locate really exists,
        ## (2) check that it finds something, (3) try to
        ## use 'default' location??
        admb_home <- gsub("/bin/admb","",system("locate bin/admb | grep bin/admb$",intern=TRUE))
        ## n.b. extra slash at end of ADMB_HOME is **VERY BAD** **VERY CONFUSING**
        ##  provokes weird behavior where "bin/sedd..." turns into "binsedd..." ???
        if (length(admb_home)>1) {
          warning("'locate' found more than one instance of bin/admb: using last")
          ## FIXME: query user?
          admb_home <- admb_home[length(admb_home)]
        }
      }
      if (.Platform$OS.type=="windows") {
        ## default location from IDE setup
        admb_home <- "c:/admb"
      }
      if (nchar(admb_home)==0) stop("couldn't guess ADMB_HOME location,",
                 "you will have to configure it manually")
    }
  }
  Sys.setenv(ADMB_HOME=admb_home)
  path=Sys.getenv("PATH")
  Sys.setenv(PATH=paste(path,paste(admb_home,"bin",sep="/"),
               sep=":"))
  admb_home
}
