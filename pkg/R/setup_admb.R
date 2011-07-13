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
        ## unix.
        ## (1) see if we can find a directory that contains 'bin/admb' *and* 'LICENSE'
        admb_home1 <- suppressWarnings(system(paste("locate bin/admb | grep bin/admb$ |",
                                                    "sed -e \"s:/bin/admb$::\""),intern=TRUE))
        admb_home2 <- suppressWarnings(system("locate admb/LICENSE | sed -e \"s:/LICENSE::\"",
                                              intern=TRUE))
        admb_home <- intersect(admb_home1,admb_home2)
        if (length(admb_home)>1) {
          warning("found more than one instance of bin/admb: using last")
          ## FIXME: query user for which one to use?
          admb_home <- tail(admb_home,1)
        } else {
          ## try default location
          if (file.exists("/usr/local/admb")) {        
            admb_home <- "/usr/local/admb"
          } else {
            admb_home <- ""
          }
        }
      } else {  ## not unix
        ## n.b. extra slash at end of ADMB_HOME is **VERY BAD** **VERY CONFUSING**
        ##  provokes weird behavior where "bin/sedd..." turns into "binsedd..." ???
        if (.Platform$OS.type=="windows") {
          ## default location from IDE setup
          admb_home <- "c:/admb"
          warning("admb_home setting is dicey for Windows,",
                  "suggest manual setting")
        }
        else stop("unknown platform")
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
