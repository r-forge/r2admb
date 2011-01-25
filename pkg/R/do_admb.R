do_admb <- function(fn,
                    data,
                    ## clean up argument list: bounds, phase as attributes on data?
                    params,
                    bounds=NULL,
                    phase=NULL,    ## FIXME: implement phases
                    ## clean up argument list: re if re_vectors is not NULL/missing?
                    re=FALSE,
                    re_vectors=NULL,
                    ## clean up argument list: kluge?
                    data_type=NULL,
                    safe=TRUE,
                    ## clean up argument list: profile if profpars is not NULL/missing?
                    profile=FALSE,
                    profpars=NULL,
                    ## clean up argument list: mcmcpars, mcmc.control?
                    mcmc=FALSE,
                    mcmc2=FALSE,
                    mcmcsteps=1000,
                    mcmcsave=max(1,round(mcmcsteps/1000)),
                    mcmcpars=NULL,
                    impsamp=FALSE,
                    verbose=FALSE,
                    wd=getwd(),
                    checkparam=c("stop","warn","write","ignore"),
                    checkdata=c("stop","warn","write","ignore"),
                    ## actions=c("compile","write","run","read"),
                    ## or should this be compile=TRUE, write=TRUE, run=TRUE, read=TRUE?
                    write_files=TRUE,
                    objfunname="f",
                    clean=TRUE,
                    extra.args) {
  ## TO DO: check to see if executables are found
  ## MODULARIZE (separate sub-function):
  ##  1. check or construct input&data files, optionally write TPL file
  ##  2. compile (tpl -> rem/cpp -> binary)
  ##  3. run
  ##  4. retrieve & package output
  checkparam <- match.arg(checkparam)
  checkdata <- match.arg(checkdata)
  if (mcmc && mcmc2) stop("only one of mcmc and mcmc2 can be specified")
  if (mcmc && missing(mcmcpars) && checkparam=="write") stop("must specify mcmcpars when checkparam=='write' and mcmc is TRUE")
  if (!mcmc && !missing(mcmcpars)) stop("mcmcpars specified but mcmc is FALSE")
  if (profile && missing(profpars) && checkparam=="write") stop("must specify profpars when checkparam=='write' and profile is TRUE")
  if (!profile && !missing(profpars)) stop("profpars specified but profile is FALSE")
  if (!re) {
    if (mcmc2) stop("mcmc2 only applies when re=TRUE")
    if (!is.null(re_vectors)) stop("re_vectors should only be specified when re=TRUE")
  } else {
    if (checkparam=="write" && is.null(re_vectors)) {
      stop("must specify random effects vectors")
    }
  }
  tplfile <- paste(fn,"tpl",sep=".")
  if (!file.exists(tplfile)) stop("can't find TPL file ",tplfile)
  tpldat <- read_tpl(fn)  ## extract info from TPL file
  tplinfo <- tpldat$info
  orig_tplfile <- tplfile
  ofn <- fn
  ## FIXME: need to make this work on MacOS/Windows, and safe
  if (!tolower(fn)==fn) {
    warning("Base name converted to lower case for ADMB compatibility: copying TPL file")
    tplfile <- tolower(tplfile)
    fn <- tolower(fn)
    ## if (file.exists(tplfile)) stop("refusing to write over existing (lowercase) TPL file")
    file.copy(orig_tplfile,tplfile,overwrite=TRUE)
  }
  
  ## require(glmmADMB)
  ## if (!re) {
  ## system(paste("makeadms",fn))
  ## } else {
  ##  system(paste("tpl2rem",fn,"; mygcss-re",fn))
  ## }
  if (!file.exists(tplfile))
    stop("could not find TPL file ",tplfile)
  ### check PARAMETER section
  if (!checkparam %in% c("write","ignore") && is.null(tplinfo$inits))
    stop("must specify PARAMETER section (or set 'checkparam' to 'write' or 'ignore')")
  if (checkparam!="ignore") {
    dmsg <- check_section(ofn,tpldat,"inits",params,
                          check=checkparam,
                          bounds=bounds,
                          data_type=data_type,
                          secname="PARAMETER",
                          objfunname=objfunname,
                          re_vectors=re_vectors,
                          mcmcpars=mcmcpars,
                          profpars=profpars)
  }
  if (!checkparam %in% c("write","ignore") && nchar(dmsg)>0) {
    if (checkparam=="stop") stop(dmsg)
    if (checkparam=="warn") warning(dmsg)
  } else if (checkparam=="write") {
    if (!is.null(tpldat$secs$PARAMETER)) {
      tpldat$secs$PARAMETER <- dmsg
    } else {
      ## insert immediately before PROCEDURE
      tpldat$secs <- append(tpldat$secs,list(PARAMETER=c(dmsg,"")),
                            after=which(names(tpldat$secs)=="PROCEDURE")-1)
    }
    ## modifications to PROCEDURE section:
    ## need to assign MCMC reporting variables
    if (mcmc) {
      
      mcmcparnames <- gsub("^ +sdreport_(number|vector) r_","",
                           gsub("\\(.*$","",
                                dmsg[grep("^ +sdreport",dmsg)]))
      tpldat$secs$PROCEDURE <- append(tpldat$secs$PROCEDURE,
                                      indent(paste("r_",mcmcparnames,"=",mcmcparnames,";",sep="")))
    }
    if (profile) {
      profparnames <- gsub("^ +likeprof_number p_","",
                           gsub("\\(.*$","",
                                dmsg[grep("^ +likeprof_",dmsg)]))
      tpldat$secs$PROCEDURE <- append(tpldat$secs$PROCEDURE,
                                      indent(paste("p_",profparnames,"=",profparnames,";",sep="")))
    }
  }
  ## check DATA section
  if (checkdata!="ignore") {
    dframes <- sapply(data,data.class)=="data.frame"
    if (any(dframes)) warning("attempted to convert data frame to matrix")
    data[dframes] <- lapply(data[dframes],as.matrix)
  }                              
  if (!checkdata %in% c("write","ignore") && is.null(tplinfo$data))
    stop("must specify DATA section (or set 'checkdata' to 'write' or 'ignore')")
  if (checkdata != "ignore") {
    dmsg <- check_section(ofn,tpldat,"data",data,
                          check=checkdata,
                          data_type=data_type,
                          secname="DATA")
    if (!checkdata %in% c("write","ignore") && nchar(dmsg)>0) {
      if (checkdata=="stop") stop(dmsg)
      if (checkdata=="warn") warning(dmsg)
      
    } else if (checkdata=="write") {
      if (!is.null(tpldat$secs$DATA)) {
        tpldat$secs$DATA <- dmsg
      } else {
        ## insert immediately before PARAMETER
        tpldat$secs <- append(tpldat$secs,list(DATA=c(dmsg,"")),
                              after=which(names(tpldat$secs)=="PARAMETER")-1)
      }
    }
  }
  ##
  if (checkdata=="write" || checkparam=="write") {
    parnames <- c(names(data),names(params))
    badnames <- grep("\\.",parnames)
    if (length(badnames)>0) {
      for (i in badnames) {
        old <- parnames[i]
        new <- gsub("\\.","_",parnames[i])
        tpldat$secs$PROCEDURE <- gsub(old,new,tpldat$secs$PROCEDURE)
      }
    }
    ## fn2 <- paste(fn,".tpl",sep="")
    fn <- paste(fn,"_gen",sep="")
    tplfile <- paste(fn,".tpl",sep="")
    ## fn2bak<- paste(fn,".tpl.bak",sep="")
    ## file.copy(fn2,fn2bak,overwrite=TRUE)
    ## FIXME: work with generated file, don't touch original!
    ## on exit, copy auto-generated file and restore original ...
    ## on.exit(file.copy(fn2,fn2gen,overwrite=TRUE),add=TRUE)
    ## on.exit(file.copy(fn2bak,fn2,overwrite=TRUE),add=TRUE)
    writeLines(do.call("c",tpldat$secs),con=tplfile)
    tpldat <- read_tpl(fn) ## get auto-generated info
  }
  if (write_files) {
    if (verbose) cat("writing data and parameter files ...\n")
    ## check order of data; length of vectors???
    dat_write(fn,data)
    ## check order of parameters ??
    ## add random effects to list of initial parameters
    if (re) {
      rv <- re_vectors[!names(re_vectors) %in% names(params)]
      params <- c(params,lapply(as.list(rv),rep,x=0))
    }
    pin_write(fn,params)
  }
  ## insert check(s) for failure at this point
  ## PART 2A: compile
  test <- try(system("admb",intern=TRUE))
  if (inherits(test,"try-error")) stop("base admb command failed: run setup_admb(), or check ADMB installation")
  args <- ""
  if (re) args <- "-r"
  if (safe) args <- paste(args,"-s")
  if (verbose) cat("compiling with args: '",args,"' ...\n")
  res0 <- system(paste("admb",args,fn," 2>",paste(fn,".cout",sep="")),
    intern=TRUE)
  coutfile <- readLines(paste(fn,".cout",sep=""))
  if (verbose) {
    cat("compile output:\n",res0,"\n")
    cat("compile log:\n")
    cat(coutfile,sep="\n")
  }
  ##  if (length(grep("error",coutfile)>0))
  ## HACK for ignorable ADMB-RE error messages
  admb_warning_index <- grep("warning:",coutfile)
  if (length(admb_warning_index)>0) {
    admb_warnings <- paste("from ADMB:",coutfile[admb_warning_index])
    sapply(admb_warnings,warning)
    coutfile <- coutfile[-admb_warning_index]
  }
  cred <- coutfile[!substr(coutfile,1,85) %in%
                   c("cat: xxalloc4.tmp: No such file or directory",
                     "cat: xxalloc5.tmp: No such file or directory",
                     "Error executing command cat xxglobal.tmp   xxhtop.tmp   header.tmp   xxalloc1.tmp   x")]
  if (length(cred)>0)
    stop("errors detected in compilation: run with verbose=TRUE to view")
  ## PART 2B: run executable file
  args <- ""
  if (mcmc) {
    args <- paste(args,"-mcmc",mcmcsteps)
    if (mcmcsave>0)
      args <- paste(args,"-mcsave",mcmcsave)
  }
  if (profile) args <- paste(args,"-lprof")
  if (!missing(extra.args)) {
    args <- paste(args,extra.args)
  }
  if (verbose) cat("running compiled executable with args: '",args,"'...\n")
  res <- system(paste("./",fn,args," 2>",fn,".out",sep=""),intern=TRUE)
  outfile <- readLines(paste(fn,".out",sep=""))
  ## replace empty res with <empty> ?
  if (verbose) {
    cat("Run output:\n",res,"\n",sep="\n")
    cat(outfile,"\n",sep="\n")
  }
  if (length(grep("^Error",outfile)>0))
    stop("errors detected in run: run with verbose=TRUE to view")
  ## PART 3
  if (verbose) cat("reading output ...\n")
  L <- c(list(fn=fn,txt=res),read_pars(fn))
  if (mcmc) {
    if (checkparam!="write") {
      warning("MCMC naming is probably wrong")
    }
    pnames <- get_names(c(mcmcpars,names(re_vectors)),
                        tpldat$info)
    L <- c(L,list(hist=read_hst(fn)))
    if (mcmcsave>0) {
      L$mcmc <- read_psv(fn,names=pnames)
      attr(L$mcmc,"mcpar") <- c(1,mcmcsteps,mcmcsave)
    }
    ## FIXME: if TPL file is user-written we need to recover
    ## the *order* of mcmc pars somehow?
  }
  if (profile) {
    L$prof <- lapply(paste("p_",profpars,sep=""),read_plt)
    names(L$prof) <- profpars
  }
  if (isTRUE(clean)) clean <- "all"
  ## need to do this **AFTER** auto-generated versions get swapped around
  if (is.character(clean)) {
    ## cover both cases
    on.exit(clean_admb(fn,clean),add=TRUE)
    on.exit(clean_admb(fn,clean,profpars),add=TRUE)
  }
  ## check for NA/NaN in logLik, errors in text?
  class(L) <- "admb"
  L
}
