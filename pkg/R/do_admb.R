run.control <- function(check_tpl=TRUE,
                        write_files=TRUE,
                        checkparam=c("stop","warn","write","ignore"),
                        checkdata=c("stop","warn","write","ignore"),
                        compile=TRUE,
                        run=TRUE,
                        read_files=TRUE,
                        clean_files="all") {
   checkparam <- match.arg(checkparam)
   checkdata <- match.arg(checkdata)
   if (is.logical(clean_files)) {
     if (clean_files) clean_files <- "all"
     else clean_files <- "none"
   }
   c(check_tpl=check_tpl,checkparam=checkparam,checkdata=checkdata,
     write_files=write_files,compile=compile,run=run,read_files=read_files,
     clean_files=clean_files)
 }

 mcmc.control <- function(mcmc=1000,
                          mcmc2=0,
                          mcsave,
                          mcnoscale=FALSE,
                          mcgrope=FALSE,
                          mcmult=1,
                          mcmcpars=NULL) {
   if (missing(mcsave)) mcsave <- pmax(1,floor(mcmc/1000))
   if (mcmc2>0) {
       if (missing(mcmc)) {
           mcmc <- 0
       }
       if (mcmc>0) stop("may not specify both mcmc and mcmc2>0")
   }
   r <- list(mcsave=mcsave,mcnoscale=mcnoscale,mcgrope=mcgrope,mcmult=mcmult,mcmcpars=mcmcpars)
   if (mcmc>0) c(list(mcmc=mcmc),r) else c(list(mcmc2=mcmc2),r)
 }

 mcmc.args <- function(L) {
   L[["mcmcpars"]] <- NULL ## don't want to include this
  argstr <- mapply(function(n,val) {
    if (is.numeric(val)) paste("-",n," ",val,sep="") else
    if (isTRUE(val)) paste("-",val,sep="")
  },names(L),L)
  paste(unlist(argstr),collapse=" ")
}

test_OScase <- function(dir=getwd()) {
  fn1 <- tempfile(tmpdir=dir)
  fn2 <- gsub("/([^/]+)$","/\\U\\1",fn1,perl=TRUE)
  if (!file.create(fn1)) stop("can't create temporary file")
  ## attempt to copy: will return FALSE if 'file already exists', i.e.
  ##  file system is case-insensitive.  Could also use file.rename(), but
  ##  would have to suppress warning msg
  res <- file.copy(fn1,fn2)
  unlink(fn1)
  unlink(fn2)
  res
}
 
do_admb <- function(fn,
                    data,
                    ## clean up argument list: bounds, phase as attributes on data?
                    params,
                    bounds=NULL,
                    phase=NULL,    ## FIXME: implement phases
                    re=NULL,
                    ## clean up argument list: kluge?
                    data_type=NULL,
                    safe=TRUE,
                    ## clean up argument list: profile if profpars is not NULL/missing?
                    profile=FALSE,
                    ## FIXME: replace profpars with profile.opts? or 
                    profpars=NULL,
                    mcmc=FALSE,
                    mcmc.opts=mcmc.control(),
                    impsamp=FALSE,
                    verbose=FALSE,
                    run.opts=run.control(),
                    objfunname="f",
                    workdir=getwd(),
                    admb_errors=c("stop","warn","ignore"),
                    extra.args) {
  ## TO DO: check to see if executables are found
  ## MODULARIZE (separate sub-function):
  ##  1. check or construct input&data files, optionally write TPL file
  ##  2. compile (tpl -> rem/cpp -> binary)  [DONE]
  ##  3. run
  ##  4. retrieve & package output
  admb_errors <- match.arg(admb_errors)
  if (!missing(workdir)) {
    file.copy(list.files(pattern=paste(fn,"\\.(dat|pin|tpl)",sep="")),
              workdir)
    cwd <- setwd(workdir)
    on.exit(setwd(cwd))
  }
  checkparam <- run.opts["checkparam"]
  checkdata <- run.opts["checkdata"]
  if (mcmc) {
    if (checkparam=="write" && !"mcmcpars" %in% names(mcmc.opts))
      stop("must specify mcmcpars when checkparam=='write' and mcmc is TRUE")
  }
  if (profile & !is.null(re)) {
    stop("profiling is not implemented for models with random effects")
  }
  if (profile && missing(profpars) && checkparam=="write")
    stop("must specify profpars when checkparam=='write' and profile is TRUE")
  if (!profile && !missing(profpars)) stop("profpars specified but profile is FALSE")
  if (is.null(re)) {
    if (mcmc && "mcmc2" %in% names(mcmc.opts)) stop("mcmc2 can only be used with random-effects models")
  }
  tplfile <- paste(fn,"tpl",sep=".")
  if (!file.exists(tplfile)) stop("can't find TPL file ",tplfile)
  if (run.opts["check_tpl"]) {
    tpldat <- read_tpl(fn)  ## extract info from TPL file
    tplinfo <- tpldat$info
    orig_tplfile <- tplfile
    ofn <- fn
    ## FIXME: need to make this work on MacOS/Windows, and safe
    if (test_OScase() && !tolower(fn)==fn) {
      warning("Base name converted to lower case for ADMB compatibility: copying TPL file")
      tplfile <- tolower(tplfile)
      fn <- tolower(fn)
      ## if (file.exists(tplfile)) stop("refusing to write over existing (lowercase) TPL file")
      file.copy(orig_tplfile,tplfile,overwrite=TRUE)
    }
    if (!file.exists(tplfile))
      stop("could not find TPL file ",tplfile)
    ## check PARAMETER section
    if (!checkparam %in% c("write","ignore") && is.null(tplinfo$inits))
      stop("must specify PARAMETER section (or set 'checkparam' to 'write' or 'ignore')")
    if (checkparam!="ignore") {
      dmsg <- check_section(ofn,tpldat,"inits",params,
                          check=checkparam,
                          bounds=bounds,
                          data_type=data_type,
                          secname="PARAMETER",
                          objfunname=objfunname,
                          re=re,
                          mcmcpars=mcmc.opts$mcmcpars,
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
      if (length(mcmcparnames)>0) {
        tpldat$secs$PROCEDURE <- append(tpldat$secs$PROCEDURE,
                                        indent(paste("r_",mcmcparnames,"=",mcmcparnames,";",sep="")))
      }
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
  }
  if (run.opts["write_files"]) {
    if (verbose) cat("writing data and parameter files ...\n")
    ## check order of data; length of vectors???
    dat_write(fn,data)
    ## check order of parameters ??
    ## add random effects to list of initial parameters
    if (!is.null(re)) {
      rv <- re[!names(re) %in% names(params)]
      params <- c(params,lapply(as.list(rv),rep,x=0))
    }
    pin_write(fn,params)
  }
  ## insert check(s) for failure at this point
  ## PART 2A: compile
  if (run.opts["compile"]) {
    compile_admb(fn,safe,re=!is.null(re),verbose)
  }
  ## PART 2B: run executable file
  if (run.opts["run"]) {
    res <- run_admb(fn,verbose,mcmc,mcmc.opts,profile,extra.args,admb_errors)
  } else res <- NULL
  if (run.opts["read_files"]) {
    ## pnames <- get_names(c(mcmc.opts[["mcmcpars"]],names(re)),
    ## tpldat$info)
    ## PART 3
    ## FIXME: we should be able to recover info about prof and MCMC par names from
    ##  somewhere ... re-read TPL files etc.?
    L <- read_admb(fn,verbose,
                   profile,
                   mcmc,
                   admbOut=res)
  }
  ## need to do this **AFTER** auto-generated versions get swapped around
  ## cover both cases
  ## FIXME: check -- why is this here? it's "on.exit" -- can't it go earlier?
  on.exit(clean_admb(fn,run.opts["clean_files"]),add=TRUE)
  ## check for NA/NaN in logLik, errors in text?
  if (run.opts["read_files"]) L else NULL
}
