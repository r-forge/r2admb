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
        admb_home <- suppressWarnings(system("locate bin/admb | grep bin/admb$",intern=TRUE))
        if (length(admb_home)>0) {
          admb_home <- gsub("/bin/admb$","",admb_home)
        } else if (file.exists("/usr/local/admb")) {
          ## try default location
          admb_home <- "/usr/local/admb"
        } else {
          admb_home <- ""
        }
        ## n.b. extra slash at end of ADMB_HOME is **VERY BAD** **VERY CONFUSING**
        ##  provokes weird behavior where "bin/sedd..." turns into "binsedd..." ???
        if (length(admb_home)>1) {
          warning("'locate' found more than one instance of bin/admb: using last")
          ## FIXME: query user for which one to use?
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
  path <- Sys.getenv("PATH")
  if (.Platform$OS.type=="windows") {
      ## FIXME: don't know if this is general enough ?
      ##  are empty elements in path (;;) OK?
      ##  what happens if compiler is not found ... ?
      pathsepchr <-  ";"
      compiler_path <- "C:/MinGW/bin"
      pathstr <- paste(c(paste(admb_home,c("bin","utilities"),sep="/"),
                         if (file.exists(compiler_path)) compiler_path else "",
                         path),collapse=pathsepchr)
      Sys.setenv(PATH=pathstr)
  } else {
      ## assume that compiler etc. are already in PATH ...
      pathsepchr <-  ":"
      Sys.setenv(PATH=paste(path,paste(admb_home,"bin",sep="/"),
                 sep=pathsepchr))
  }
  admb_home
}
