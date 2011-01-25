indent <- function(str,n=2) {
  paste(paste(rep(" ",n),collapse=""),str,sep="")
}

## format numbers equal width with leading zeros if necessary
numfmt <- function(x,len=length(x)) {
  paste(x,
        formatC(seq(len),width=format.info(seq(len)),flag="0"),
        sep="")
}

read_pars <- function (fn) {
  ## see
  ##  http://admb-project.org/community/admb-meeting-march-29-31/InterfacingADMBwithR.pdf
  ## for an alternate file reader -- does this have equivalent functionality?
  rt <- function(f,ext,...) {
    fn <- paste(f,ext,sep=".")
    if (file.exists(fn)) read.table(fn,...) else NA
  }
  rs <- function(f,ext,comment.char="#",...) {
    fn <- paste(f,ext,sep=".")
    if (file.exists(fn)) scan(fn,
                              comment.char=comment.char,
                              quiet=TRUE,...) else NA
  }
  ## parameter estimates
  par_dat <- rs(fn,"par", skip = 1)
  npar <- length(par_dat)
  tmp <- rs(fn, "par", what = "", comment.char="")
  ## COULD get parnames out of par file, but a big nuisance
  ##  for vectors etc.
  ## parnames <- gsub(":$","",tmp[seq(18,by=3,length=npar)])
  loglik <- as.numeric(tmp[11])
  maxgrad <- as.numeric(tmp[16])
  ## second pass to extract names from par file (ugh)
  tmp2 <- readLines(paste(fn,".par",sep=""))
  parlines <- grep("^#",tmp2)[-1]
  parlen <- count.fields(paste(fn,".par",sep=""))
  parnames <- gsub("^# +","",gsub(":$","",tmp2[parlines]))
  parnames <- unname(unlist(mapply(function(x,len) {
    if (len==1) x else paste(x,1:len,sep=".")
  },
                                   parnames,parlen)))
  est <- unlist(par_dat)
  names(est) <- parnames
  if (!is.finite(loglik)) warning("bad log-likelihood: fitting problem in ADMB?")
  ## if non-pos-def hessian, cor and std files will be missing ... but
  ##   we should still be able to retrieve some info
  sd_dat <- rt(fn,"std", skip = 1,as.is=TRUE)
  if (length(sd_dat)==1 && is.na(sd_dat)) {
    warning("std file missing: some problem with fit, but retrieving parameter estimates anyway")
    cormat <- vcov <- matrix(NA,nrow=npar,ncol=npar)
    std <- rep(NA,npar)
  } else {
    ## need col.names hack so read.table knows how many
    ##  columns to read: ?read.table, "Details"
    ncorpar <- length(readLines(paste(fn,"cor",sep=".")))-2
    cor_dat <- rt(fn,"cor", skip = 2, fill=TRUE, 
                  as.is=TRUE,col.names=paste("X",1:(4+ncorpar),sep=""))
    ## drop cors that are not parameters
    ## (have dropped mc parameters)
    cormat <- as.matrix(cor_dat[1:npar,4+(1:npar)])
    cormat[upper.tri(cormat)] <- t(cormat)[upper.tri(cormat)]
    parnames <- sd_dat[1:npar, 2]  ## FIXME: check with parnames above
    if (any(duplicated(parnames))) {
      parnames <- unlist(lapply(split(parnames,factor(parnames)),
                                function(x) {
                                  if (length(x)==1) x else numfmt(x)
                                }))
    }
    std <- sd_dat[1:npar, 4]
    vcov <- outer(std,std) * cormat
  }
  names(std) <- rownames(vcov) <- rownames(cormat) <-
    colnames(vcov) <- colnames(cormat) <- parnames
  list(coefficients=est, se=std, loglik=-loglik, maxgrad=-maxgrad, cor=cormat, vcov=vcov)
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
    mcpar <- attr(x$mcmc,"mcpar")
    cat("MCMC: start=",mcpar[1],", end=",mcpar[2],", thin=",mcpar[3],"\n")
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

stdEr <- function(x, ...) {
  UseMethod("stdEr")
}

coef.admb <- function(object,...) object$coefficients
logLik.admb <- function(object,...) object$loglik
vcov.admb <- function(object,...) object$vcov
stdEr.admb <- function(object,...) sqrt(diag(object$vcov))
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

## processing variables
proc_var <- function(s,drop.first=TRUE,maxlen) {
  if (drop.first) s <- s[-1]
  ## drop LOCAL CALCS sections
  calclocs <- grep("_CALCS *$",s)
  if (length(calclocs)>0) {
    droplines <- unlist(apply(matrix(-calclocs,nrow=2,byrow=TRUE),1,function(x) x[1]:x[2]))
    s <- s[droplines]
  }
  ## strip comments & whitespace
  s2 <- gsub("^[ \\\t]*","",gsub("[;]*[ \\\t]*$","",strip_comments(s)))
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
  t1 <- textConnection(r[3:(cisecline[1]-1)])
  prof1 <- matrix(scan(t1,quiet=TRUE),ncol=2,
                  byrow=TRUE,
                  dimnames=list(NULL,c("value","logLik")))
  close(t1)
  t2 <- textConnection(r[cisecline[1]+(2:4)])
  ci1 <- matrix(scan(t2,quiet=TRUE),ncol=3,
                  byrow=TRUE,
                  dimnames=list(NULL,c("sig","lower","upper")))
  close(t2)
  t3 <- textConnection(r[(normline+1):(cisecline[2]-1)])
  profnorm <- matrix(scan(t3,quiet=TRUE),ncol=2,
                  byrow=TRUE,
                  dimnames=list(NULL,c("value","logliK")))
  close(t3)
  t4 <- textConnection(r[cisecline[2]+(2:4)])
  cinorm <- matrix(scan(t4,quiet=TRUE),ncol=3,
                  byrow=TRUE,
                  dimnames=list(NULL,c("sig","lower","upper")))
  close(t4)
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

## confint.default works
confint.admb <- function(object, parm, level=0.95, method="default", ...) {
  if (method %in% c("default","quad")) {
    confint.default(object)
  } else if (method=="profile") {
    vals <- object[["prof"]]
    if (is.null(vals)) stop("model not fitted with profile=TRUE")
    if (!level %in% c(0.9,0.95,975)) stop("arbitrary levels not yet implemented:",
                                          "level must be in (0.9,0.95,0.975)")
    tab <- t(sapply(vals,function(x) {
      x$ci[x$ci[,"sig"]==level,c("lower","upper")]
    }))
    colnames(tab) <- paste(c((1-level)/2,(1+level)/2)*100,"%")
    tab
  } else if (method %in% c("quantile","HPDinterval")) {
    vals <- object[["mcmc"]]
    if (is.null(vals)) stop("model not fitted with mcmc=TRUE")
    if (method=="quantile") {
      tab <- t(apply(vals,2,quantile,c((1-level)/2,(1+level)/2)))
    } else {
      require(coda)
      tab <- HPDinterval(as.mcmc(vals))
      colnames(tab) <- paste(c((1-level)/2,(1+level)/2)*100,"%")
    }
  }
}
