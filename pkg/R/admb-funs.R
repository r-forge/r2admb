indent <- function(str,n=2) {
  paste(paste(rep(" ",n),collapse=""),str,sep="")
}

## format numbers equal width with leading zeros if necessary
numfmt <- function(x,len=length(x)) {
  paste(x,
        formatC(seq(len),width=format.info(seq(len)),flag="0"),
        sep="")
}

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

check_section <- function(fn,
                          tpldat,
                          tplsec,
                          R_list,
                          check,
                          bounds,
                          phase,
                          re_vectors,
                          mcmcpars,
                          profpars,
                          secname,
                          objfunname="f",
                          intcheck=c("strict","sloppy")) {
  intcheck <- match.arg(intcheck)
  Rnames  <- names(R_list)
  msg <- ""
  if (check!="write") {
    info <- tpldat$info[[tplsec]]
    tplnames <- info$vname
    if (length(setdiff(tplnames,Rnames))>0) {
      msg <- paste("missing values in list:",
                   paste(setdiff(tplnames,Rnames),sep=","))
    } else if (length(setdiff(Rnames,tplnames))==0 && !all(tplnames==Rnames)) {
      msg <- "all values present, but order doesn't match"
    } else {
      msg <- ""
      if (length(grep("\\.",tplnames)>1)) {
        msg <- paste(msg,"dots in parameter/variable names")
      }
      attach(R_list,warn.conflicts=FALSE)
      on.exit(detach(R_list),add=TRUE)
      ## now need to check dimensions etc...
      for (i in 1:nrow(info)) {
        v <- info[i,]
        ## x <- get(v$vname)
        ##  i.e. search only in data, not everywhere ...
        x <- R_list[[v$vname]]
        v$type <- gsub("init_","",v$type)
        if (v$type %in% c("int","ivector","imatrix")) {
          if (any(trunc(x)!=x)) msg <- paste(msg,v$vname,
                     "non-integer;")
        }
        if (v$type %in% c("int","number")) {
          if (length(x)>1) msg <- paste(msg,"length(",v$vname,
                                        ")>1;")
        }
        if (v$type %in% c("ivector","vector")) {
          if (any(is.na(c(v$X1,v$X2)))) {
            msg <- paste(msg,"NAs in dimensions in ",v$vname)
          } else {
            tpllen <- eval(parse(text=paste(v$X2,"-",v$X1)))+1
            if (length(x)!=tpllen)
              msg <- paste(msg,"length mismatch in ",v$vname,
                           ": ",length(x)," (r) != ",tpllen," (tpl)",
                           sep="")
          }
        }
        if (v$type %in% c("imatrix","matrix")) { 
          tpldim <- with(v,c(
            eval(parse(text=paste(v$X2,"-",v$X1)))+1,
            eval(parse(text=paste(v$X4,"-",v$X3)))+1))
          rdim <- dim(x)
          if (is.null(rdim)) {
            msg <- paste(msg,v$vname,"not a matrix;")
          } else {
            if (any(rdim!=tpldim))
            msg <- paste(msg,"dimension mismatch in ",v$vname,
                         ": (",rdim[1],",",rdim[2],"), (r) != ",
                         " (",tpldim[1],",",tpldim[2],") (tpl)",
                         sep="")
          } ## dimensions
        } ## if matrix
        if (length(grep(v$type,"array"))>0) {
          arraydim <- as.numeric(substr(v$type,1,1))
          tpldim <- numeric(arraydim)
          for (j in 1:arraydim) {
            tpldim[j] <-
              eval(parse(text=paste(v[2*j+2],"-",v[2*j+1])))+1
          }
          rdim <- dim(x)
          if (is.null(rdim)) {
            msg <- paste(msg,v$vname,"not an array;")
          } else {
            if (any(rdim!=tpldim))
              msg <- paste(msg,"dimension mismatch in ",v$vname,
                         ": (",
                           paste(rdim,sep=","),"), (r) != ",
                         " (",
                           paste(tpldim,sep=","),") (tpl)",
                         sep="")
          } ## !is.null(rdim)
        } ## array
      } ## loop over variables
    } ## checking
    return(msg)
  } else { ## check=="write"
    ## WRITING THE TPL SECTION
    ## FIXME: handle bounds, phases, and 1D matrices??
    ## bounds: list OR (named) 2-col matrix or 2xn matrix
    ## phases: list OR (named) vector or n-vector
    ## arrays must be 1-based: if you want zero-based arrays
    ##     then write the TPL section yourself
    if (!missing(phase))
      stop("phase support not yet implemented")
    ## no existing section: need title line
    sectitle <- paste(secname,"_SECTION",sep="")
    nvals <- length(R_list)
    if (length(grep("\\.",names(R_list)>1))) {
      warning("dots changed to underscores in parameter/variable names")
      names(R_list) <- gsub("\\.","_",names(R_list))
    }
    objstr <- NULL
    pars <- (secname=="PARAMETER")  ## only check bounds for parameters
    if (pars) {
      if (is.null(objfunname)) stop("must specify a name for the objective function")
      objstr <- indent(paste("objective_function_value",objfunname))
    }
    parstr <- character(nvals)
    ## FIXME: parameter tables should be more standardized --
    ## room for all possible dimensions, bounds, phase?
    partab <- data.frame(type=character(nvals),
                         vname=names(R_list),
                         dim1=rep(NA,nvals),
                         dim2=rep(NA,nvals),
                         dim3=rep(NA,nvals),
                         dim4=rep(NA,nvals),
                         lower=rep(NA,nvals),
                         upper=rep(NA,nvals),
                         phase=rep(NA,nvals),
                         stringsAsFactors=FALSE)
    for (i in 1:length(R_list)) {
      x <- R_list[[i]]
      n <- names(R_list)[i]
      ## attempt to coerce (if not all numeric, will end up as character and stop ...)
      if (is.data.frame(x)) x <- as.matrix(x)
      is.int <- function(x) {
        ((intcheck=="strict" && storage.mode(x)=="integer") |
         (intcheck=="trunc" && all(trunc(x)==x)))}
      if (is.int(x) && !is.null(dim(x)) &&
          length(dim(x))>2) {
        stop("can't handle integer arrays of dimension>2")
      }
      if (is.int(x)) {
        if (length(x)==1 && is.null(dim(x))) {
          parstr[i] <- paste("init_int",n)
          partab$type[i] <- "int"
        } else if (length(x)>1 && is.null(dim(x))) {
          parstr[i] <- paste("init_ivector ",n," (1,",length(x),")",sep="")
          partab$type[i] <- "ivector"
          partab$dim1 <- length(x)
        } else if (!is.null(dim(x)) && length(dim(x))==2) {
          parstr[i] <- paste("init_imatrix",n,
                  " (1,",dim(x)[1],",1,",dim(x)[2],")")
          partab$dim1 <- dim(x)[1]
          partab$dim2 <- dim(x)[2]
        }
      } else if (storage.mode(x) %in% c("numeric","double")) {
        if (length(x)==1 && is.null(dim(x))) {
          if (pars && !is.null(bounds) && n %in% names(bounds)) {
            parstr[i] <- paste("init_bounded_number ",n,
                               "(",bounds[[n]][1],",",bounds[[n]][2],")",sep="")
          } else {
            parstr[i] <- paste("init_number",n)           
          }
          partab$type[i] <- "number"
        } else if (length(x)>1 && is.null(dim(x))) {
          if (pars && !is.null(bounds) && n %in% names(bounds)) {
            parstr[i] <- paste("init_bounded_vector ",n,
                    "(1,",length(x),bounds[[n]][1],",",bounds[[n]][2],")",sep="")
          } else {
            parstr[i] <- paste("init_vector ",n,"(1,",length(x),")",sep="")
          }
          partab$type[i] <- "vector"
          partab$dim1[i] <- length(x)
        } else if (!is.null(dim(x)) && length(dim(x))==2) {
          parstr[i] <- paste("init_matrix ",n,
                  "(1,",dim(x)[1],",1,",dim(x)[2],")",sep="")
          partab$type[i] <- "matrix"
          partab$dim1 <- dim(x)[1]
          partab$dim2 <- dim(x)[2]
        } else if (!is.null(dim(x)) && length(dim(x))>2) {
          ndim <- length(dim(x))
          if (ndim>7) stop("can't handle arrays of dim>7")
          parstr[i] <- paste("init_",ndim,"array",
                  n," (",
                  paste(c(rbind(rep(1,ndim),dim(x))),
                        collapse=","),")",sep="")
          partab$type[i] <- "array"
          ## FIXME: store array dimensions?
          if (!is.null(mcmcpars)) stop("arrays currently incompatible with MCMC")
        } ## multi-dim array
      } else stop("can only handle numeric values")
    } ## loop over R list
    parstr <- indent(parstr)
    cursec <- tpldat$secs[[secname]]
    if (!is.null(cursec)) {
      cursec <- cursec[-1] ## drop title
      cursec <- grep("^ *$",cursec,invert=TRUE,value=TRUE)  ## drop blank lines
    }
    restr <- mcmcstr <- profstr <- NULL
    if (pars) {
      ## deal with random effects vectors
      if (!is.null(re_vectors)) {
        nre <- length(re_vectors)
        restr <- character(nre)
        if (is.list(re_vectors)) re_vectors <- unlist(re_vectors)
        restr <- indent(paste("random_effects_vector ",names(re_vectors),
                "(1,",re_vectors,")",sep=""))
      }
      ## FIXME: uuuuuugly! need a better, more consistent way
      ## of handling parameter attributes ...
      make_names <- function(z,pref1="sdreport_",pref2="r_") {
        if (z %in% partab$vname) { ## in newly specified parameters
          i <- match(z,partab$vname)
          type <- partab$type[i]
          name <- partab$vname[i]
          dimvals <- na.omit(unlist(partab[i,c("dim1","dim2","dim3","dim4")]))
          dimvals <- c(rbind(rep(1,length(dimvals)),dimvals))
        } else if (z %in% tt$vname) { ## in existing parameters
          i <- match(z,tt$vname)
          type <- tt$type[i]
          name <- tt$vname[i]
          dimvals <- na.omit(unlist(tt[i,3:7]))
        }
        indent(paste(pref1,type," ",pref2,name,
                     if (length(dimvals)==0) "" else
                     paste("(",paste(dimvals,collapse=","),")",sep=""),
                     sep=""))
      }
      if(!is.null(mcmcpars) || !is.null(profpars)) {
        tt <- tpldat$info$other
        allnames <- c(partab$vname,tt$vname)
        if (!is.null(mcmcpars)) {
          bad <- which(!mcmcpars %in% allnames)
          if (length(bad)>0) {
            stop("some mcmcpars not found in parameter table:",
                 paste(mcmcpars[bad],collapse=", "))
          }
          mcmcstr <- sapply(mcmcpars,make_names,pref1="sdreport_",pref2="r_")
        }
        if (!is.null(profpars)) {
          bad <- which(!profpars %in% allnames)
          if (length(bad)>0) {
            stop("some profpars not found in parameter table:",
                 paste(profpars[bad],collapse=", "))
          }
          profstr <- sapply(profpars,make_names,pref1="likeprof_",pref2="p_")
        }
      }
    }
    return(c(sectitle,"",objstr,parstr,cursec,restr,mcmcstr,profstr))
  }
}

read_pars <- function (fn) {
  rt <- function(f,ext,...) {
    if (file.exists(f)) read.table(paste(f,ext,sep="."),...) else NA
  }
  rs <- function(f,ext,comment.char="#",...) {
    if (file.exists(f)) scan(paste(f,ext,sep="."),
                             comment.char=comment.char,quiet=TRUE,...) else NA
  }
  par_dat <- rs(fn,"par", skip = 1)
  npar <- length(par_dat)
  sd_dat <- rt(fn,"std", skip = 1,as.is=TRUE)
  if (length(sd_dat)==1 && is.na(sd_dat)) return(NA)
  ## need col.names hack so read.table knows how many
  ##  columns to read: ?read.table, "Details"
  ncorpar <- length(readLines(paste(fn,"cor",sep=".")))-2
  cor_dat <- rt(fn,"cor", skip = 2, fill=TRUE, 
                as.is=TRUE,col.names=paste("X",1:(4+ncorpar),sep=""))
  ## drop cors that are not parameters
  ## (have dropped mc parameters)
  cormat <- as.matrix(cor_dat[1:npar,4+(1:npar)])
  cormat[upper.tri(cormat)] <- t(cormat)[upper.tri(cormat)]
  est <- unlist(par_dat)
  parnames <- sd_dat[1:npar, 2]
  if (any(duplicated(parnames))) {
    parnames <- unlist(lapply(split(parnames,factor(parnames)),
                              function(x) {
                                if (length(x)==1) x else numfmt(x)
                              }))
  }
  std <- sd_dat[1:npar, 4]
  tmp <- rs(fn, "par", what = "", comment.char="")
  loglik <- as.numeric(tmp[11])
  grad <- as.numeric(tmp[16])
  vcov <- outer(std,std) * cormat 
  names(est) <- names(std) <- rownames(vcov) <- rownames(cormat) <-
    colnames(vcov) <- colnames(cormat) <- parnames
  if (!is.finite(loglik)) warning("bad log-likelihood: fitting problem in ADMB?")
  list(coefficients=est, se=std, loglik=-loglik, grad=-grad, cor=cormat, vcov=vcov)
}

do_admb <- function(fn,
                    data,
                    params,
                    bounds=NULL, ## maybe specify some other way, e.g. as attributes on data?
                    phase=NULL,
                    ## FIXME: allow phases?
                    re=FALSE,
                    re_vectors=NULL,
                    safe=TRUE,
                    profile=FALSE,
                    profpars=NULL,
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
                    objfunname="f",
                    clean=TRUE,
                    extra.args) {
  ## TO DO: check to see if executables are found
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
  tpldat <- read_tpl(fn)  ## extract info from TPL file
  tplinfo <- tpldat$info
  ofn <- fn
  if (!tolower(fn)==fn) {
    warning("base name converted to lower case for ADMB compatibility")
    fn <- tolower(fn)
    file.copy(ofn,fn,overwrite=TRUE)
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
  dmsg <- check_section(ofn,tpldat,"inits",params,
                        check=checkparam,
                        bounds=bounds,
                        secname="PARAMETER",
                        objfunname=objfunname,
                        re_vectors=re_vectors,
                        mcmcpars=mcmcpars,
                        profpars=profpars)
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
  dmsg <- check_section(ofn,tpldat,"data",data,
                        check=checkdata,
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
    fn2 <- paste(fn,".tpl",sep="")
    fn2gen <- paste(fn,"_gen.tpl",sep="")
    fn2bak<- paste(fn,".tpl.bak",sep="")
    file.copy(fn2,fn2bak,overwrite=TRUE)
    ## on exit, copy auto-generated file and restore original ...
    on.exit(file.copy(fn2,fn2gen,overwrite=TRUE),add=TRUE)
    on.exit(file.copy(fn2bak,fn2,overwrite=TRUE),add=TRUE)
    writeLines(do.call("c",tpldat$secs),con=fn2)
    tpldat <- read_tpl(fn) ## get auto-generated info
  }
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
                     "cat: xxalloc5.tmp: No such file or directory",                               "Error executing command cat xxglobal.tmp   xxhtop.tmp   header.tmp   xxalloc1.tmp   x")]
  if (length(cred)>0)
    stop("errors detected in compilation: run with verbose=TRUE to view")
  ## insert check(s) for failure at this point
  if (verbose) cat("writing data and parameter files ...\n")
  ## check order of data; length of vectors???
  dat_write(fn,data)
  ## check order of parameters ??
  ## add random effects to list of initial parameters
  if (re) {
    rv <- re_vectors[!names(re_vectors) %in% names(params)]
    params <- c(params,lapply(as.list(re_vectors),rep,x=0))
  }
  pin_write(fn,params)
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
      L$mcmc <- read_psv(ofn,names=pnames)
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
    on.exit(clean_admb(ofn,clean),add=TRUE)
    on.exit(clean_admb(fn,clean,profpars),add=TRUE)
  }
  ## check for NA/NaN in logLik, errors in text?
  class(L) <- "admb"
  L
}

str_contains <- function(x,y) {
  length(grep(x,y)>1)
}

get_names <- function(pars,info) {
  unlist(sapply(pars,
         function(p) {
           if (p %in% info$inits$vname) {
             i <- match(p,info$inits$vname)
             tt <- info$inits$type[i]
             if (str_contains("number$",tt)) {
               p
             } else if (str_contains("vector$",tt)) {
               numfmt(p,as.numeric(info$inits$X2[i]))
               ## FIXME: may fail if this value needs to be parsed?
             } else stop("can't handle matrix names yet")
           } else if (p %in% info$raneff$vname) {
             i <- match(p,info$raneff$vname)
             numfmt(p,as.numeric(info$raneff$X2[i]))
           }
         }))
}

print.admb <- function(x, verbose=FALSE, ...) {
  cat("Model file:",x$fn,"\n")
  if (is.null(x$loglik)) {
    cat("No fit\n")
    return(invisible(NULL))
  }
  cat("Negative log-likelihood:",-x$loglik,"\n")
  cat("Coefficients:\n")
  print(unlist(x$coefficients))
  if (!is.null(x$mcmc)) {
    cat("MCMC: ",nrow(x$mcmc)," steps\n")
  }
  if (verbose) cat(x$txt,sep="\n")
}      

summary.admb <- function(object, correlation=FALSE, symbolic.cor = FALSE, ...) {
  p1 <- 1:length(object$coefficients)
  coef.p <- unlist(object$coefficients)
  covmat <- object$vcov
  var.cf <- diag(covmat)
  s.err <- sqrt(var.cf)
  tvalue <- coef.p/s.err
  dn <- c("Estimate", "Std. Error")
  pvalue <- 2 * pnorm(-abs(tvalue))
  coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(names(coef.p), c(dn, 
                "z value", "Pr(>|z|)"))
  ans <- list(coefficients=coef.table,loglik=object$loglik,fn=object$fn)
  class(ans) <- "summary.admb"
  ans
}

print.summary.admb <- function(x,
                         digits = max(3, getOption("digits") - 3),
                         symbolic.cor = x$symbolic.cor, 
                         signif.stars = getOption("show.signif.stars"), ...) {
  coefs <- x$coefficients
  cat("Model file:",x$fn,"\n")
  cat("Negative log-likelihood:",-x$loglik,"\n")
  cat("Coefficients:\n")
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
}

coef.admb <- function(object,...) object$coefficients
logLik.admb <- function(object,...) object$loglik
vcov.admb <- function(object,...) object$vcov
deviance.admb <- function(object,...) -2*object$loglik
AIC.admb <- function(object,...,k=2) {
  if (length(list(...))>0) stop("multi-object AIC not yet implemented")
  deviance(object)+k*length(coef(object))
}

## summary() method ...
##  save model file with object???

clean_admb <- function(fn,which=c("all","sys","output"),profpars=NULL) {
  tplfile <- paste(fn,"tpl",sep=".")
  ## need logic here!  enumeration of all files ???
  sys.ext <- c("htp","cpp","o","rep","rhes","bar","eva",
               "bgs","ecm","luu","mc2","mcm","tpl.bak")
  input.ext <- c("pin","dat")
  output.ext <- c("log","cor","std", "par","psv","hst","prf")
  sys2.ext <- c("out","cout")
  other <- c("eigv.rpt","fmin.log","variance","sims",
             "hesscheck","hessian.bin","dgs2","diags",
             paste("admodel",c("dep","hes","cov"),sep="."))
  which <- match.arg(which)
  if (which=="all") {
    ## erase only targeted extensions -- NOT everything with the basename except the tpl file
    delfiles <-  paste(fn,c(sys.ext,input.ext,output.ext,sys2.ext),sep=".")
    ## list.files(pattern=paste("^",fn,"\\..*",sep=""))
    ##    delfiles <- setdiff(delfiles,tplfile)
    delfiles <- c(delfiles,other)
    if (!is.null(profpars)) {
      delfiles <- c(delfiles,paste(profpars,".plt",sep=""))
    }
  } else {
    stop("only 'all' option is currently implemented")
  } 
  unlink(delfiles)
}

## @rm -vf $(PROGRAM_NAME){.htp,.cpp,.std,.rep,.b??,.p??,.r??,.cor,.eva,.log,.rhes,.luu,}
## @rm -vf admodel{.cov,.dep,.hes} eigv.rpt fmin.log variance hess*

read_chunk <- function(fn,sep="^#",maxlines=1000) {
  end <- FALSE
  ans <- character(maxlines)
  i <- 1
  has_sep <- function(x) length(grep(sep,x))>0
  while (!end) {
    tmp <- readLines(fn,n=1)
    if (i>1 && has_sep(tmp)) {
      end=TRUE
      pushBack(tmp,fn)
    } else if (length(tmp)==0) {
      end=TRUE
    } else {
      ans[i] <- tmp
      i <- i+1
    }
  }
  ans[1:(i-1)]
}

read_hst <- function(fn) {
  fn <- paste(fn,"hst",sep=".")
  if (!file.exists(fn)) {
    warning("file ",fn," not found: returning NULL")
    return(NULL)
  }
  f <- file(fn,open="r")
  r <- list()
  repeat {
    chunk <- read_chunk(f)
    ## cat(length(chunk),":",chunk[1],"\n") 
    if (length(chunk)<2) break
    r <- c(r,list(chunk))
  }
  labs <- sapply(r,"[[",1)
  ## single values
  w <- c(1:2,8,10)
  r[w] <- lapply(r[w],
                 function(x) as.numeric(x[2]))
  names(r)[w] <- c("sampsize","stepsize_scale","npars","rseed")
  ## vectors of value
  w <- c(3:7,9)
  r[w] <- lapply(r[w],function(x) 
                 as.numeric(strsplit(gsub("^ +","",x[2])," ")[[1]]))
  names(r)[w] <- c("stepsizes","means","sdevs","lower","upper","mcmcparms")
  ## r$npars is NOT RELIABLE! use length(stepsizes instead)
  ## parameter matrices
  w <- 11:(10+length(r$stepsizes))
  r[w] <- lapply(r[w],
                 function(z) {
                   do.call(rbind,
                           lapply(z[-c(1,length(z))],
                                  function(x) {as.numeric(strsplit(x," ")[[1]])}))})
  names(r)[w] <- gsub("^#","",
                      gsub("\\[([0-9]+)\\]",".\\1",
                           gsub("; *//.*","",labs[w])))
  ans <- c(r[1:10],hists=list(r[w]))
  class(ans) <- "admb_hist"
  ans
}


## don't know how to structure this properly.
##  I would like to have plot.admb_hist plot
##  a graph (as a side effect) and invisibly return a
##  data frame restructured for plotting,
##  but that causes problems when incorporating lattice and ggplot
##  plots into Sweave documents: one would normally say
##  print(plotfun(...)), but that doesn't work when plotfun()
##  doesn't actually return the plot (I think)
## I broke the restructuring function out so it can be
## used separately, but at the moment am still 
## returning the restructured data.  A possible workaround for
## ggplot() plots:
##  plotfun(...); print(last_plot())
rhist <- function(x,pars) {
  nbars <-sapply(x,nrow)
  xx <- data.frame(do.call(rbind,x),
                   param=rep(names(x),nbars))
  if (!missing(pars)) {
    if (is.numeric(pars))
      pars <- levels(xx$param)[pars]
      xx <- subset(xx,as.character(xx$param) %in% pars)
  }
  xx
}

plot.admb_hist <- function(x,type=c("lattice","ggplot"),
                           dtype=c("hist","density"),
                           pars,
                           ...) { ## dots for generic compat
  type <- match.arg(type)
  dtype <- match.arg(dtype)
  if (dtype=="hist") {
    xx <- rhist(x$hists,pars)
    if (type=="ggplot") {
      ## if (!require(ggplot2)) stop("must install ggplot2 package")
      X1 <- ""; X2 <- ""  ## hack to circumvent NOTE in R CMD check
      vplot <- ggplot2::ggplot(xx,aes(x=X1,y=X2))+
        geom_step()+
          ## geom_bar(stat="identity",fill="darkgray")+
          facet_wrap(~param,scales="free")+
            labs(y="Frequency",x="")
    } else if (type=="lattice") {
      require(lattice)
      ##barchart(X2~X1|param,data=xx,horiz=FALSE,
      vplot <- xyplot(X2~X1|param,type="s",data=xx,
        scales=list(x=list(relation="free",tick.number=3),
          y=list(relation="free")),
        as.table=TRUE)
    }
  }
  vplot
  ## invisible(xx)
}


## from Steve Martell
reptoRlist <- function(fn) {
  ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
  idx=sapply(as.double(ifile),is.na)
  vnam=ifile[idx] #list names
  nv=length(vnam) #number of objects
  A=list()
  r=0
  for(i in 1:nv) {
    ir=match(vnam[i],ifile)
    if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
    dum=NA
    if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
    if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
    if(is.numeric(dum))#Logical test to ensure dealing with numbers
      {
        A[[ vnam[i ]]]=dum
      }
  }
  return(A)
}

if (FALSE) {
  readBin("mccoypred6.b01","double",n=11)
  scan("mccoypred6.p01",comment="#")
  zz <- matrix(readBin("mccoypred6.mcm","double",n=6000),
    byrow=TRUE,ncol=6)
  zz <- as.data.frame(zz)
  names(zz) <- c("c","d","h","g","sigma_c","xx")
  zz <- data.frame(it=1:1000,zz)
  library(reshape)
  mz <- melt(zz,id=1)
  library(lattice)
  xyplot(value~it|variable,data=mz,type="l",as.table=TRUE,
         scales=list(y=list(relation="free")))
  ##
  zz <- as.data.frame(read_admbbin("mccoypred6.psv"))
  names(zz) <- c("c","d","h","g","sigma_c",paste("c",1:6,sep=""))
  zz <- data.frame(it=1:1000,zz)
  mz <- melt(zz[1:6],id=1)
  xyplot(value~it|variable,data=mz,type="l",as.table=TRUE,
         scales=list(y=list(relation="free")))
  ggplot(mz,aes(x=it,y=value))+geom_line()+
    facet_wrap(~variable,scale="free_y")
  ggplot(mz,aes(x=value))+geom_density()+
    facet_wrap(~variable,scale="free")
}


strip_comments <- function(s) {
  ## strip comments (and terminal whitespace)
  gsub("[ \\\t]*//.*$","",s)
}

proc_var <- function(s,drop.first=TRUE,maxlen) {
  if (drop.first) s <- s[-1]
  ## strip comments & whitespace
  s2 <- gsub("^ *","",gsub("[;]*[ \\\t]*$","",strip_comments(s)))
  s2 <- s2[nchar(s2)>0] 
  words <- strsplit(s2," ")
  type <- sapply(words,"[[",1)
  rest <- sapply(words,"[[",2)
  rest2 <- strsplit(gsub("[(),]"," ",rest)," ")
  vname <- sapply(rest2,"[[",1)
  maxlen0 <- max(sapply(rest2,length))
  if (missing(maxlen)) maxlen <- maxlen0
  else maxlen <- pmax(maxlen,maxlen0)
  opts <- t(sapply(rest2,
           function(w) {
             ## as.numeric()?
             c(w[-1],rep(NA,maxlen+1-length(w)))
           }))
  data.frame(type,vname,opts,stringsAsFactors=FALSE)
}

drop_calcs <- function(s) {
  startcalc <- grep("^ *LOC_CALCS",s)
  endcalc <- grep("^ *END_CALCS",s)
  if (length(endcalc)==0) endcalc <- length(s)
  if (length(startcalc)>0) {
    s <- s[-(startcalc:endcalc)]
  }
  commcalc <- grep("^ +!!",s)
  if (length(commcalc)>0) s <- s[-commcalc]
  s
}

read_tpl <- function(f) {
  r <- readLines(paste(f,"tpl",sep="."))
  secStart <- which(substr(r,1,1) %in% LETTERS)
  if (length(secStart)==0) stop("tpl file must contain at least one section (capitalized header)")
  if (secStart[1]!=1) { ## add first (comments etc.) section
    secStart <- c(1,secStart)
  }
  nsec <- length(secStart)
  ## length (in lines) of each chunk
  L <- c(secStart[-1],length(r)+1)-secStart
  sec <- rep(1:nsec,L)
  splsec <- split(r,sec)
  ## not QUITE right: we get some stuff in here that is not
  ##  a "SECTION" but is SEPARABLE_FUNCTION or TOP_OF_MAIN_CALCS
  ##  or something ...
  splnames <- sapply(splsec,"[",1)
  names(splsec) <- gsub("_.+","",splnames)
  splsec_proc <- lapply(splsec,drop_calcs)
  L1 <- L2 <- NULL
  pp <- splsec_proc$PARAMETER
  if (!is.null(pp)) {
      pp <- proc_var(pp,maxlen=7)
      type <- 1 ## kluge for R CMD check warnings; will be masked
      L1 <- with(pp,
             list(inits=pp[grep("^init",type),],
                  raneff=pp[grep("^random",type),],
                  sdnums=pp[grep("^sdreport_number",type),],
                  sdvecs=pp[grep("^sdreport_vector",type),],
                  other=pp[grep("^init|random|sdreport",type,invert=TRUE),],
                  ## FIXME: don't know what I needed this for
                  ## sdvecdims <- gsub("^ +sdreport_vector[ a-zA-Z]+","",
                  ## gsub("[()]","",
                  ## grep( "^ +sdreport_vector",splsec$PARAMETER,
                  ## value=TRUE)))
                  profparms=pp[grep("^likeprof",type),]))
    }
  pp <- splsec_proc$DATA
  if (!is.null(pp)) {
    pp <- proc_var(pp,maxlen=7)
    L2 <- with(pp,
               list(data=pp[grep("^init",type),]))
  }
  L <- c(L1,L2)
  L <- L[!sapply(L,is.null)]
  list(secs=splsec,info=L[sapply(L,nrow)>0])
}

read_psv <- function(f,names=NULL) {
  f <- tolower(f) ## arghv
  fn <- paste(f,"psv",sep=".")
  if (!file.exists(fn)) stop("no PSV file found")
  ans <- read_admbbin(fn)
  colnames(ans) <- names
  ans <- as.data.frame(ans)
  ans
}

read_plt <- function(varname) {
  fn <- paste(varname,"plt",sep=".")
  r <- readLines(fn)
  cisecline <- grep("Minimum width confidence limits",r)
  normline <- grep("Normal approximation$",r)
  prof1 <- matrix(scan(textConnection(r[3:(cisecline[1]-1)]),quiet=TRUE),ncol=2,
                  byrow=TRUE,
                  dimnames=list(NULL,c("value","logliK")))
  ci1 <- matrix(scan(textConnection(r[cisecline[1]+(2:4)]),quiet=TRUE),ncol=3,
                  byrow=TRUE,
                  dimnames=list(NULL,c("sig","lower","upper")))
  profnorm <- matrix(scan(textConnection(r[(normline+1):(cisecline[2]-1)]),quiet=TRUE),ncol=2,
                  byrow=TRUE,
                  dimnames=list(NULL,c("value","logliK")))
  cinorm <- matrix(scan(textConnection(r[cisecline[2]+(2:4)]),quiet=TRUE),ncol=3,
                  byrow=TRUE,
                  dimnames=list(NULL,c("sig","lower","upper")))
  list(prof=prof1,ci=ci1,prof_norm=profnorm,ci_norm=cinorm)
}
## read a "standard" ADMB format binary file into R:
##  standard format is: 1 integer describing number
##  of (double) values per vector
read_admbbin <- function(fn) {
  f <- file(fn,open="rb")
  nv <- readBin(f,"int")
  fs <- file.info(fn)$size
  isize <- 4; dsize <- 8
  m <- matrix(readBin(f,"double",n=(fs-isize)/dsize),byrow=TRUE,
         ncol=nv)
  m
}
         
## from glmmADMB, by Hans Skaug
"dat_write" <-
function (name, L) 
{
    n <- nchar(name)
    if (substring(name, n - 3, n) == ".dat") {
      file_name <- name
    } else file_name <- paste(name, "dat", sep = ".")
    cat("# \"", file_name,"\" produced by dat_write() from R2admb ", 
        date(), "\n", file = file_name, sep = "")
    for (i in 1:length(L)) {
        x <- L[[i]]
        if (data.class(x) == "numeric") 
            cat("#", names(L)[i], "\n", L[[i]], "\n\n", file = file_name, 
                append = TRUE)
        if (data.class(x) == "matrix") {
            cat("#", names(L)[i], "\n", file = file_name, append = TRUE)
            write.table(L[[i]], , col = FALSE, row = FALSE, quote = FALSE, 
                file = file_name, append = TRUE)
            cat("\n", file = file_name, append = TRUE)
        }
    }
  }

## from glmmADMB, by Hans Skaug
"pin_write" <-
function (name, L) 
{
    n <- nchar(name)
    if (substring(name, n - 3, n) == ".pin") 
        file_name <- name
    else file_name <- paste(name, ".pin", sep = "")
    cat("# \"", name, ".pin\" produced by pin_write() from R2admb ", 
        date(), "\n", file = file_name, sep = "")
    for (i in 1:length(L)) {
        x <- L[[i]]
        if (data.class(x) == "numeric") 
            cat("#", names(L)[i], "\n", L[[i]], "\n\n", file = file_name, 
                append = TRUE)
        if (data.class(x) == "matrix") {
            cat("#", names(L)[i], "\n", file = file_name, append = TRUE)
            write.table(L[[i]], , col = FALSE, row = FALSE, quote = FALSE, 
                file = file_name, append = TRUE)
            cat("\n", file = file_name, append = TRUE)
        }
    }
  }


if (FALSE) {
  ## test: can we read all ADMB examples without crashing?
  dir <- "/usr/local/src/admb/examples/admb/"
  dir <- "/usr/local/src/admb/examples/admb-re/"
  setwd(dir)
  ## omit files with '.' (happen to be non-directories)
  L <- list.files(pattern="^[a-zA-Z_]+$")
  source("/home/ben/lib/R/pkgs/r2admb/pkg/R/admb-funs.R")
  for (i in seq_along(L)) {
    setwd(file.path(dir,L[i]))
    tpls <- gsub(".tpl","",list.files(pattern=".tpl"))
    for (j in seq_along(tpls)) {
      cat(L[i],tpls[j],"\n")
      invisible(read_tpl(tpls[j])$info)
    }
  }
}

