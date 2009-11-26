## TO DO:
##   bounds checking?
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
  list(coef=est, se=std, loglik=-loglik, grad=-grad, cor=cormat, vcov=vcov)
}

do_admb = function(fn,input_list,param_list,
  re=FALSE,
  safe=TRUE,
  mcmc=FALSE,
  mcmc2=FALSE,
  mcmcsteps=1000,
  mcmcsave=round(mcmcsteps/1000),
  impsamp=FALSE,
  verbose=FALSE,
  wd=getwd(),
  clean=FALSE) {
  ## TO DO: check to see if executables are found
  if (mcmc && mcmc2) stop("only one of mcmc and mcmc2 can be specified")
  if (!re && mcmc2) stop("mcmc2 only applies when re=TRUE")
  tplfile = paste(fn,"tpl",sep=".")
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
  if (length(coutfile)>0)
    stop("errors detected in compilation: run with verbose=TRUE to view")
  ## insert check(s) for failure at this point
  if (verbose) cat("writing data and input files ...\n")
  ## check order of data; length of vectors???
  dat_write(fn,input_list)
  ## check order of parameters ??
  pin_write(fn,param_list)
  args <- ""
  if (mcmc) {
    args = paste(args,"-mcmc",mcmcsteps)
    if (mcmcsave>0)
      args = paste(args,"-mcsave",mcmcsave)
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
  print(unlist(x$coef))
  if (verbose) cat(x$txt,sep="\n")
}      

summary.admb <- function(object, correlation=FALSE, symbolic.cor = FALSE, ...) {
  p1 <- 1:length(object$coef)
  coef.p <- unlist(object$coef)
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

coef.admb <- function(object,...) object$coef
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

read_tpl <- function(f) {
  r <- readLines(paste(f,"tpl",sep="."))
  secStart <- which(substr(r,1,1) %in% LETTERS)
  nsec <- length(secStart)
  L <- c(secStart[-1],length(r)+1)-secStart
  sec <- rep(1:nsec,L)
  splsec <- split(r,sec)
  names(splsec) <- gsub("_[A-Z]+","",sapply(splsec,"[",1))
  inits <- gsub("^ +init[a-z_]+ ","",
                gsub("\\([)a-zA-Z0-9,. ]+$","",
                     grep( "^ +init",splsec$PARAMETER,
                          value=TRUE)))
  raneff <- gsub("^ +random[a-z_]+ ","",
                 gsub("\\([)a-zA-Z0-9,. ]+$","",
                      grep( "^ +random",splsec$PARAMETER,
                           value=TRUE)))
  sdlines <- grep( "^ +sdreport_number",splsec$PARAMETER,
                           value=TRUE)
  sdnums <- gsub("^ +sdreport_number +","",
                 gsub("//.*$","",sdlines))
  sdnums <- strsplit(gsub(" *","",
                          paste(sdnums,collapse="")),";")[[1]]
  sdvecs <- gsub("^ +sdreport_vector +","",
                 ## gsub("\\([)a-zA-Z0-9,.;/ ]+$","",
                 gsub(";.*$","",
                           grep( "^ +sdreport_vector",splsec$PARAMETER,
                                value=TRUE)))
  sdvecdims <- gsub("^ +sdreport_vector[ a-zA-Z]+","",
                 gsub("[()]","",
                      grep( "^ +sdreport_vector",splsec$PARAMETER,
                           value=TRUE)))

  profparms <- gsub("^ +likeprof[a-z_]+ ","",
                 gsub("\\([)a-zA-Z0-9,. ]+$","",
                      grep( "^ +likeprof",splsec$PARAMETER,
                           value=TRUE)))

  list(inits=inits,raneff=raneff,
       sdnums=sdnums,sdvecs=sdvecs,profparms=profparms)
}

read_psv <- function(f) {
  tpl <- read_tpl(f)
  f <- tolower(f) ## argh
  fn <- paste(f,"psv",sep=".")
  if (!file.exists(fn)) stop("no PSV file found")
  ans <- read_admbbin(fn)
  nnums <- length(tpl$sdnums)
  colnames(ans) <- c(tpl$sdnums,rep("",ncol(ans)-nnums))
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
    cat("# \"", file_name,"\" produced by dat_write() from ADMButils; ", 
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
    cat("# \"", name, ".pin\" produced by pin_write() from ADMButils; ", 
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
