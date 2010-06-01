## TO DO:
##   finish bounds/phase support
##    bounds checking?
##   integrate writing/reading/checking
##   add objective_function_value default
##   check for random effects vectors (don't redefine)
##   integrate sdreport stuff for MCMC
##   parameter/data order checking?
##   more checks/flags of compiling step -- stop if error



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
        admb_home <- gsub("bin/admb","",system("locate bin/admb",intern=TRUE))
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
                          tplsec,
                          R_list,
                          check,
                          bounds,
                          phase,
                          secname,
                          intcheck=c("strict","sloppy")) {
  intcheck <- match.arg(intcheck)
  Rnames  <- names(R_list)
  msg <- ""
  if (check!="write") {
    info <- read_tpl(fn)[[tplsec]]
    tplnames <- info$vname
    if (length(setdiff(tplnames,Rnames))>0) {
      msg <- paste("missing values in list:",
                   paste(setdiff(tplnames,Rnames),sep=","))
    } else if (length(setdiff(Rnames,tplnames))==0 && !all(tplnames==Rnames)) {
      msg <- "all values present, but order doesn't match"
    } else {
      msg <- ""
      attach(R_list,warn.conflicts=FALSE)
      on.exit(detach(R_list))
      ## now need to check dimensions etc...
      for (i in 1:nrow(info)) {
        v = info[i,]
        x = get(v$vname)
        v$type = gsub("init_","",v$type)
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
            tpllen = eval(parse(text=paste(v$X2,"-",v$X1)))+1
            if (length(x)!=tpllen)
              msg <- paste(msg,"length mismatch in ",v$vname,
                           ": ",length(x)," (r) != ",tpllen," (tpl)",
                           sep="")
          }
        }
        if (v$type %in% c("imatrix","matrix")) {
          tpldim = with(v,c(
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
          arraydim = as.numeric(substr(v$type,1,1))
          tpldim = numeric(arraydim)
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
    if (!missing(bounds) || !missing(phase))
      warning("bounds and phase support incomplete")
    secstr <- paste(secname,"_SECTION",sep="")
    for (i in 1:length(R_list)) {
      x <- R_list[[i]]
      n <- names(R_list)[i]
      is.int <- function(x) {
        ((intcheck=="strict" && storage.mode(x)=="integer") |
         (intcheck=="trunc" && all(trunc(x)==x)))}
      if (is.int(x) && !is.null(dim(x)) &&
          length(dim(x))>2) {
        stop("can't handle integer arrays of dimension>2")
      }
      if (is.int(x)) {
        if (length(x)==1 && is.null(dim(x))) {
          secstr[i+1] = paste("init_int",n)
        } else if (length(x)>1 && is.null(dim(x))) {
          secstr[i+1] = paste("init_ivector ",n," (1,",length(x),")",sep="")
        } else if (!is.null(dim(x)) && length(dim(x))==2) {
          secstr[i+1] = paste("init_imatrix",n,
                  " (1,",dim(x)[1],",1,",dim(x)[2],")")
        }
      } else if (storage.mode(x) %in% c("numeric","double")) {
        if (length(x)==1 && is.null(dim(x))) {
          if (n %in% names(bounds)) {
            secstr[i+1] = paste("init_bounded_number ",n,
                    "(",c(bounds[[n]][1],",",bounds[[n]][2]),")",sep="")
          } else {
            secstr[i+1] = paste("init_number",n)
          }
        } else if (length(x)>1 && is.null(dim(x))) {
          secstr[i+1] = paste("init_vector ",n,"(1,",length(x),")",sep="")
        } else if (!is.null(dim(x)) && length(dim(x))==2) {
          secstr[i+1] = paste("init_matrix ",n,
                  " (1,",dim(x)[1],",1,",dim(x)[2],")",sep="")
        } else if (!is.null(dim(x)) && length(dim(x))>2) {
          ndim <- length(dim(x))
          if (ndim>7) stop("can't handle arrays of dim>7")
          secstr[i+1] = paste("init_",ndim,"array",
                  n," (",
                  paste(c(rbind(rep(1,ndim),dim(x))),
                        collapse=","),")",sep="")
        } ## multi-dim array
      } else stop("can only handle numeric values")
    } ## loop over R list
    fn2 <- paste(fn,".tpl",sep="")
    file.copy(fn2,
              paste(fn,".tpl.bak",sep=""))
    ff <- readLines(fn2)
    ff <- c(secstr,"",ff)
    writeLines(ff,con=fn2)
    return(msg)
  }
}

read_pars <- function (fn) {
  rt <- function(f,...) {
    if (file.exists(f)) read.table(f,...) else NA
  }
  par_dat = rt(paste(fn,"par",sep="."), skip = 1)
  npar = nrow(par_dat)
  sd_dat = rt(paste(fn,"std",sep="."), skip = 1,as.is=TRUE)
  if (length(sd_dat)==1 && is.na(sd_dat)) return(NA)
  ## need col.names hack so read.table knows how many
  ##  columns to read: ?read.table, "Details"
  ncorpar = length(readLines(paste(fn,"cor",sep=".")))-2
  cor_dat = rt(paste(fn,"cor",sep="."), skip = 2, fill=TRUE, 
    as.is=TRUE,col.names=paste("X",1:(4+ncorpar),sep=""))
  ## drop cors that are not parameters
  ## (have dropped mc parameters)
  cormat = as.matrix(cor_dat[1:npar,4+(1:npar)])
  cormat[upper.tri(cormat)] = t(cormat)[upper.tri(cormat)]
  est = unlist(par_dat)
  parnames = sd_dat[1:npar, 2]
  std = sd_dat[1:npar, 4]
  tmp = scan(paste(fn, ".par", sep = ""), what = "", quiet = TRUE)
  loglik = as.numeric(tmp[11])
  grad = as.numeric(tmp[16])
  vcov = outer(std,std) * cormat 
  names(est) = names(std) = rownames(vcov) = rownames(cormat) =
    colnames(vcov) = colnames(cormat) = parnames
  if (!is.finite(loglik)) warning("bad log-likelihood: fitting problem in ADMB?")
  list(coefficients=est, se=std, loglik=-loglik, grad=-grad, cor=cormat, vcov=vcov)
}

do_admb = function(fn,
  data_list,param_list,
  re=FALSE,
  safe=TRUE,
  mcmc=FALSE,
  mcmc2=FALSE,
  mcmcsteps=1000,
  mcmcsave=round(mcmcsteps/1000),
  mcmcpars,
  impsamp=FALSE,
  verbose=FALSE,
  wd=getwd(),
  checkparam=c("stop","warn","write","ignore"),
  checkdata=c("stop","warn","write","ignore"),
  clean=FALSE,
  extra.args) {
  ## TO DO: check to see if executables are found
  checkparam <- match.arg(checkparam)
  checkdata <- match.arg(checkdata)
  if (mcmc && mcmc2) stop("only one of mcmc and mcmc2 can be specified")
  if (!re && mcmc2) stop("mcmc2 only applies when re=TRUE")
  tplfile = paste(fn,"tpl",sep=".")
  tplinfo <- read_tpl(fn)
  ofn <- fn
  if (!tolower(fn)==fn) {
    warning("base name converted to lower case for ADMB compatibility")
    fn <- tolower(fn)
  }
  ## require(glmmADMB)
  ## if (!re) {
  ## system(paste("makeadms",fn))
  ## } else {
  ##  system(paste("tpl2rem",fn,"; mygcss-re",fn))
  ## }
  if (!file.exists(tplfile))
    stop("could not find TPL file ",tplfile)
  dmsg <- check_section(ofn,"data",data_list,
                        check=checkdata,
                        ## lower,upper,
                        secname="PARAMETER")
  if (nchar(dmsg)>0) {
    if (checkdata=="stop") stop(dmsg)
    if (checkdata=="warn") warning(dmsg)
  }
  dmsg <- check_section(ofn,"inits",param_list,
                        check=checkparam,
                        ## lower,upper,
                        secname="")
  if (nchar(dmsg)>0) {
    if (checkparam=="stop") stop(dmsg)
    if (checkparam=="warn") warning(dmsg)
  }
  args = ""
  if (re) args = "-r"
  if (safe) args = paste(args,"-s")
  if (verbose) cat("compiling with args: '",args,"' ...\n")
  res0 = system(paste("admb",args,ofn," 2>",paste(fn,".cout",sep="")),
    intern=TRUE)
  coutfile <- readLines(paste(fn,".cout",sep=""))
  if (verbose) {
    cat("compile output:\n",res0,"\n")
    cat("compile log:\n")
    cat(coutfile,sep="\n")
  }
  ##  if (length(grep("error",coutfile)>0))
  ## HACK for ignorable ADMB-RE error messages
  cred <- coutfile[!substr(coutfile,1,85) %in%
                   c("cat: xxalloc4.tmp: No such file or directory",
                     "cat: xxalloc5.tmp: No such file or directory",                               "Error executing command cat xxglobal.tmp   xxhtop.tmp   header.tmp   xxalloc1.tmp   x")]
  if (length(cred)>0)
    stop("errors detected in compilation: run with verbose=TRUE to view")
  ## insert check(s) for failure at this point
  if (verbose) cat("writing data and parameter files ...\n")
  ## check order of data; length of vectors???
  dat_write(fn,data_list)
  ## check order of parameters ??
  pin_write(fn,param_list)
  args <- ""
  if (mcmc) {
    args = paste(args,"-mcmc",mcmcsteps)
    if (mcmcsave>0)
      args = paste(args,"-mcsave",mcmcsave)
  }
  if (!missing(extra.args)) {
    args <- paste(args,extra.args)
  }
  if (verbose) cat("running compiled executable with args: '",args,"'...\n")
  res = system(paste("./",ofn,args," 2>",fn,".out",sep=""),intern=TRUE)
  outfile <- readLines(paste(fn,".out",sep=""))
  ## replace empty res with <empty> ?
  if (verbose) {
    cat("Run output:\n",res,"\n",sep="\n")
    cat(outfile,"\n",sep="\n")
  }
  if (length(grep("^Error",outfile)>0))
    stop("errors detected in run: run with verbose=TRUE to view")
  if (verbose) cat("reading output ...\n")
  L = c(list(fn=fn,txt=res),read_pars(fn))
  if (mcmc) {
    L = c(L,list(hist=read_hst(fn),
      mcmc=read_psv(ofn)))
  }
  ## check for NA/NaN in logLik, errors in text?
  class(L) = "admb"
  if (isTRUE(clean)) clean <- "all"
  if (is.character(clean)) {
    ## cover both cases
    clean_admb(ofn,clean)
    clean_admb(fn,clean)
  }
  L
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

clean_admb <- function(fn,which=c("all","sys","output")) {
  tplfile <- paste(fn,"tpl",sep=".")
  ## need logic here!  enumeration of all files ???
  sys.ext <- c("htp","cpp","o","rep","rhes","bar","eva")
  input.ext <- c("pin","dat")
  output.ext <- c("log","cor","std", "par")
  sys2.ext <- c("out","cout")
  other <- c("eigv.rpt","fmin.log","variance","sims",
             paste("admodel",c("dep","hes","cov"),sep="."))
  which <- match.arg(which)
  if (which=="all") {
    delfiles <- list.files(pattern=paste("^",fn,"\\..*",sep=""))
    delfiles <- setdiff(delfiles,tplfile)
    delfiles <- c(delfiles,other)
  } else {
    stop("only 'all' option is currently implemented")
  } 
  unlink(delfiles)
}

## @rm -vf $(PROGRAM_NAME){.htp,.cpp,.std,.rep,.b??,.p??,.r??,.cor,.eva,.log,.rhes,.luu,}
## @rm -vf admodel{.cov,.dep,.hes} eigv.rpt fmin.log variance hess*

read_chunk <- function(fn,sep="^#",maxlines=1000) {
  end = FALSE
  ans = character(maxlines)
  i = 1
  has_sep <- function(x) length(grep(sep,x))>0
  while (!end) {
    tmp = readLines(fn,n=1)
    if (i>1 && has_sep(tmp)) {
      end=TRUE
      pushBack(tmp,fn)
    } else if (length(tmp)==0) {
      end=TRUE
    } else {
      ans[i] = tmp
      i = i+1
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
  ## parameter matrices
  w <- 11:(10+r$npars)
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
reptoRlist = function(fn) {
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
    if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=T,what=""))
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
  zz = matrix(readBin("mccoypred6.mcm","double",n=6000),
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

proc_var <- function(s,drop.first=TRUE) {
  if (drop.first) s <- s[-1]
  ## strip comments & whitespace
  s2 <- gsub("^ *","",gsub("[;]*[ \\\t]*$","",strip_comments(s)))
  s2 <- s2[nchar(s2)>0] 
  words <- strsplit(s2," ")
  type <- sapply(words,"[[",1)
  rest <- sapply(words,"[[",2)
  rest2 <- strsplit(gsub("[(),]"," ",rest)," ")
  vname <- sapply(rest2,"[[",1)
  maxlen <- max(sapply(rest2,length))
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
  names(splsec) <- gsub("_.+","",sapply(splsec,"[",1))
  splsec_proc <- lapply(splsec,drop_calcs)
  splsec_proc <- lapply(splsec_proc[c("PARAMETER","DATA")],proc_var)
  pp <- splsec_proc$PARAMETER
  type <- 1 ## kluge for R CMD check warnings; will be masked
  L1 <- with(pp,
             list(inits=pp[grep("^init",type),],
                  raneff=pp[grep("^random",type),],
                  sdnums=pp[grep("^sdreport_number",type),],
                  sdvecs=pp[grep("^sdreport_vector",type),],
                  ## FIXME: don't know what I needed this for
                  ## sdvecdims <- gsub("^ +sdreport_vector[ a-zA-Z]+","",
                  ## gsub("[()]","",
                  ## grep( "^ +sdreport_vector",splsec$PARAMETER,
                  ## value=TRUE)))
                  profparms=pp[grep("^likeprof",type),]))
  pp <- splsec_proc$DATA
  L2 <- with(pp,
             list(data=pp[grep("^init",type),]))
  L <- c(L1,L2)
  L[sapply(L,nrow)>0]
}

read_psv <- function(f) {
  tpl <- read_tpl(f)
  f <- tolower(f) ## argh
  fn <- paste(f,"psv",sep=".")
  if (!file.exists(fn)) stop("no PSV file found")
  ans <- read_admbbin(fn)
  nnums <- nrow(tpl$sdnums)
  colnames(ans) <- c(tpl$sdnums$vname,rep("",ncol(ans)-nnums))
  ## assume parameters come before random effects?
  ans <- as.data.frame(ans)
  ans
}

## read a "standard" ADMB format binary file into R:
##  standard format is: 1 integer describing number
##  of (double) values per vector
read_admbbin <- function(fn) {
  f = file(fn,open="rb")
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
    n = nchar(name)
    if (substring(name, n - 3, n) == ".dat") {
      file_name <- name
    } else file_name <- paste(name, "dat", sep = ".")
    cat("# \"", file_name,"\" produced by dat_write() from R2admb ", 
        date(), "\n", file = file_name, sep = "")
    for (i in 1:length(L)) {
        x = L[[i]]
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
    n = nchar(name)
    if (substring(name, n - 3, n) == ".pin") 
        file_name = name
    else file_name = paste(name, ".pin", sep = "")
    cat("# \"", name, ".pin\" produced by pin_write() from R2admb ", 
        date(), "\n", file = file_name, sep = "")
    for (i in 1:length(L)) {
        x = L[[i]]
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
  dir = "/usr/local/src/admb/examples/admb/"
  dir = "/usr/local/src/admb/examples/admb-re/"
  setwd(dir)
  ## omit files with '.' (happen to be non-directories)
  L = list.files(pattern="^[a-zA-Z_]+$")
  source("/home/ben/lib/R/pkgs/r2admb/pkg/R/admb-funs.R")
  for (i in seq_along(L)) {
    setwd(file.path(dir,L[i]))
    tpls <- gsub(".tpl","",list.files(pattern=".tpl"))
    for (j in seq_along(tpls)) {
      cat(L[i],tpls[j],"\n")
      invisible(read_tpl(tpls[j]))
    }
  }
}

