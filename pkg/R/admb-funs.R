indent <- function(str,n=2) {
  paste(paste(rep(" ",n),collapse=""),str,sep="")
}

## format numbers equal width with leading zeros if necessary
numfmt <- function(x,len=length(x),sep=".") {
  paste(x,
        formatC(seq(len),width=format.info(seq(len)),flag="0"),
        sep=sep)
}

## format matrix elements equal width
numfmt2 <- function(x,xdim,sep=".",sep2=".") {
  ff1 <- format.info(seq(xdim[1]))
  ff2 <- format.info(seq(xdim[2]))
  paste(x,
        outer(seq(xdim[1]),seq(xdim[2]),
              function(i,j) paste(formatC(i,width=ff1,flag="0"),
                                  formatC(j,width=ff2,flag="0"),sep=sep2)),
        sep=sep)
}

numfmt3 <- function(x,lens,sep=".") {
  ff1 <- format.info(seq(length(lens)))
  ff2 <- format.info(unlist(lens))
  paste(x,
        formatC(rep(seq(length(lens)),lens),width=ff1,flag="0"),
        formatC(unlist(sapply(lens,seq)),   width=ff2,flag="0"),
        sep=sep)
}
  

rep_pars <- function(parnames) {
  parnames <- as.character(parnames)
  parnames <- unlist(lapply(split(parnames,factor(parnames,levels=unique(parnames))),
                            function(x) {
                              if (length(x)==1) x else numfmt(x)
                            }))
  parnames
}

read_pars <- function (fn,drop_phase=TRUE) {
    ## read .par files from ADMB
    ## see
    ##  http://admb-project.org/community/admb-meeting-march-29-31/InterfacingADMBwithR.pdf
    ## for an alternate file reader -- does this have equivalent functionality?
    ## FIXME: get hessian.bin ?
    # safe read.table/scan functions
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
    par_dat <- rs(fn,"par", skip = 1)
    tmp <- rs(fn, "par", what = "", comment.char="")
    ## COULD get parnames out of par file, but this is a big nuisance for vectors etc.
    ## parnames <- gsub(":$","",tmp[seq(18,by=3,length=npar)])
    loglik <- as.numeric(tmp[11])
    maxgrad <- as.numeric(tmp[16])
    ## second pass to extract names from par file (ugh)
    tmp2 <- readLines(paste(fn,".par",sep=""))
    parlines <- grep("^#",tmp2)[-1]
    npar <- length(par_dat)      ## TOTAL number of numeric values
    npar2 <- length(parlines)  ## number of distinct parameters
    parlen <- count.fields(paste(fn,".par",sep=""))
    parlen2 <- count.fields(paste(fn,".par",sep=""),comment.char="")
    parnames0 <- parnames <- gsub("^# +","",gsub(":$","",tmp2[parlines]))
    parlist <- vector("list",npar2)
    parnameslist <- vector("list",npar2)
    names(parlist) <- parnames
    cumpar <- 1
    cumline <- 1
    ## browser()
    pp <- c(parlines,length(tmp2)+1)
    ## reshape parameters properly
    parid <- numeric(npar2)
    for (i in seq(npar2)) {
        nrows <- diff(pp)[i]-1
        curlines <- cumline:(cumline+nrows-1)
        curlen <- sum(parlen[curlines])
        parvals <- par_dat[cumpar:(cumpar+curlen-1)]
        if (nrows==1) {
            if (curlen==1) {
                parnameslist[[i]] <- parnames[i]
            } else {
                parnameslist[[i]] <- numfmt(parnames[i],curlen)
            }
            parlist[[i]] <- parvals
        } else {
            if (length(unique(parlen[curlines]))==1) {
                ## matrix
                parlist[[i]] <- matrix(parvals,nrow=nrows,byrow=TRUE)
                parnameslist[[i]] <- numfmt2(parnames[i],dim(parlist[[i]]))
            } else {
                ## ragged array
                parlist[[i]] <- split(parvals,rep(seq(curlines),parlen[curlines]))
                parnameslist[[i]] <- numfmt3(parnames[i],parlen[curlines])
            }
        }
        ## cat(parnames[i],cumline,cumline+nrows-1,cumpar,cumpar+curlen-1,parnameslist[[i]],"\n")
        ## print(parlist[[i]])
        cumline <- cumline + nrows
        cumpar <- cumpar + curlen
    }
    ## FIXME: watch out, 'short param names' has now been overwritten by 'long param names'
    parnames <- unlist(parnameslist)
    
    ## parnames <- unname(unlist(mapply(function(x,len) {
    ## if (len==1) x else numfmt(x,len)  ## paste(x,1:len,sep=".")
    ##},
    ## parnames,parlen)))
    est <- unlist(par_dat)
    names(est) <- parnames[1:length(est)]
    if (!is.finite(loglik)) warning("bad log-likelihood: fitting problem in ADMB?")
    ## if non-pos-def hessian, cor and std files will be missing ... but
    ##   we should still be able to retrieve some info
    sdmiss <- FALSE
    sd_dat <- rt(fn,"std", skip = 1,as.is=TRUE)
    if (length(sd_dat)==1 && is.na(sd_dat)) {
        sdmiss <- TRUE
        warning("std file missing: some problem with fit, but retrieving parameter estimates anyway")
        cormat <- vcov <- matrix(NA,nrow=npar,ncol=npar)
        std <- rep(NA,npar)
        sdrptvals <- numeric(0)
    } else {
        nsdpar <- nrow(sd_dat)
        ## need col.names hack so read.table knows how many
        ##  columns to read: ?read.table, "Details"
        ## FIXME: go gracefully if .cor missing?
        ncorpar <- length(readLines(paste(fn,"cor",sep=".")))-2
        cor_dat <- rt(fn,"cor", skip = 2, fill=TRUE, 
                      as.is=TRUE,col.names=paste("X",1:(4+ncorpar),sep=""))
        ## drop cors that are not parameters
        ## (have dropped mc parameters)
        cormat <- as.matrix(cor_dat[1:nsdpar,4+(1:nsdpar)])
        cormat[upper.tri(cormat)] <- t(cormat)[upper.tri(cormat)]
        ## be careful here -- need to adjust for phase<0 parameters,
        ##  which will be in parameter vector but not in
        ##  sd
        if (!sdmiss) {
            sdparnames <- sd_dat[, 2]
            misspars <- setdiff(parnames0,sdparnames)
            ## only names of positive-phase parameters
            parnames2 <- unlist(parnameslist[!parnames0 %in% misspars])
            sdparnames <- c(parnames2,sdparnames[-seq_along(parnames2)])
            ## parnames <- c(parnames,sd_dat[-seq_along(parnames),2])
            if (any(duplicated(sdparnames))) {
                sdparnames <- rep_pars(sdparnames)
            }
            npar3 <- length(parnames2) ## positive-phase only
            if (drop_phase) {
                parlist <- parlist[!parnames0 %in% misspars]
                est <- unlist(parlist)
                names(est) <- parnames2
                npar <- npar3
            }
        }
        std <- sd_dat[, 4]
        sdrptvals <- sd_dat[-(1:npar3),3]
        vcov <- outer(std,std) * cormat
        
        ## hes <- read_admbbin("admodel.hes")
        ## FIXME: can read this, but I don't know what it means!
        ##  it doesn't seem to be the raw Hessian ...
        names(std) <- rownames(vcov) <- rownames(cormat) <-
            colnames(vcov) <- colnames(cormat) <- sdparnames
    }
    list(coefficients=c(est,sdrptvals),
         coeflist=parlist,
         se=std, loglik=-loglik, maxgrad=-maxgrad, cor=cormat, vcov=vcov,
         npar=npar)
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
    print(coef(x))
    ## FIXME: indicate extra parameters?
    if (!is.null(x$mcmc)) {
        mcpar <- attr(x$mcmc,"mcpar")
        cat("MCMC parameters: start=",mcpar[1],", end=",mcpar[2],", thin=",mcpar[3],"\n",sep="")
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

stdEr <- function(object, ...) {
    UseMethod("stdEr")
}

coef.admb <- function(object,type=c("par","extra","all"),...) {
    type <- match.arg(type)
    x <- object$coefficients
    n <- object$npar
    if (is.null(n)) n <- length(x)
    switch(type,par=x[1:n],
           extra=x[-(1:n)],
           all=x)
}
logLik.admb <- function(object,...) object$loglik
vcov.admb <- function(object,type=c("par","extra","all"),...) {
    type <- match.arg(type)
    v <- object$vcov
    n <- object$npar
    if (is.null(n)) n <- ncol(v)
    switch(type,par=v[1:n,1:n,drop=FALSE],
           extra=v[-(1:n),-(1:n),drop=FALSE],
           all=v)
}
stdEr.admb <- function(object,type=c("par","extra","all"),...) {
    type <- match.arg(type)  
    s <- sqrt(diag(object$vcov))
    n <- object$npar
    switch(type,par=s[1:n],
           extra=s[-(1:n)],
           all=s)
}

deviance.admb <- function(object,...) -2*object$loglik
AIC.admb <- function(object,...,k=2) {
    if (length(list(...))>0) stop("multi-object AIC not yet implemented")
    deviance(object)+k*length(coef(object))
}

## summary() method ...
##  save model file with object???

clean_admb <- function(fn,which=c("sys","output")) {
    if (length(which)==1) {
        if (which=="none") return()
        if (which=="all") which <- c("sys","input","output","gen")
    }
    
    sys.ext <- c("bar","bgs","cpp","ecm","eva","htp","luu","mc2","mcm","o","rep","rhes",
                 "luu","mc2","mcm","tpl.bak","out","cout","shess")
    sys.files <- paste(fn,sys.ext,sep=".")
    gen.files <- list.files(pattern="_gen(\\.tpl)*")
    sys.other <- c("eigv.rpt","fmin.log","variance","sims",
                   "hesscheck","hessian.bin","dgs2","diags",
                   paste("admodel",c("dep","hes","cov"),sep="."),
                   list.files(pattern="xx.*.tmp"),
                   list.files(pattern=".*f1b2list.*"),
                   list.files(pattern=paste(fn,"\\.[bpr][0-9]+",sep="")))
    ## FIXME: clean up abandoned 'buffer' files too
    ## f1b2list etc.
    input.ext <- c("pin","dat")
    input.files <- paste(fn,input.ext,sep=".")
    output.ext <- c("log","cor","std", "par","psv","hst","prf","mcinfo")
    output.files <- paste(fn,output.ext,sep=".")
    output.files <- c(output.files,list.files(pattern="\\.plt$"))
    if ("sys" %in% which) unlink(c(sys.files,sys.other))
    if ("input" %in% which) unlink(input.files)
    if ("output" %in% which) unlink(output.files)
    if ("gen" %in% which) unlink(gen.files)
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
            end <- TRUE
            pushBack(tmp,fn)
        } else if (length(tmp)==0) {
            end <- TRUE
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
            stop("ggplot disabled to avoid dependency")
            ## if (!require(ggplot2)) stop("must install ggplot2 package")
            ## X1 <- ""; X2 <- ""  ## hack to circumvent NOTE in R CMD check
            ## vplot <- ggplot2::ggplot(xx,aes(x=X1,y=X2))+
            ##   geom_step()+
            ##     ## geom_bar(stat="identity",fill="darkgray")+
            ##     facet_wrap(~param,scales="free")+
            ##       labs(y="Frequency",x="")
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
        if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrows=irr-ir-1,fill=TRUE))
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
    ## ggplot(mz,aes(x=it,y=value))+geom_line()+
    ##   facet_wrap(~variable,scale="free_y")
    ## ggplot(mz,aes(x=value))+geom_density()+
    ##   facet_wrap(~variable,scale="free")
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
        droplines <- unlist(apply(matrix(-calclocs,
                                         ncol=2,byrow=TRUE),1,function(x) seq(x[1],x[2])))
        s <- s[droplines]
    }
    ## strip comments & whitespace
    s2 <- gsub("^[ \\\t]*","",gsub("[;]*[ \\\t]*$","",strip_comments(s)))
    s2 <- s2[nchar(s2)>0]
    s2 <- s2[!grepl("+[ \\\t]*!!",s2)] ## strip !! lines
    words <- strsplit(s2," ")
    words <- lapply(words,function(x) x[nzchar(x)]) ## fix from Jeff Laake
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
    names(splsec) <- ifelse(grepl("SECTION$",splnames),
                            gsub("_.+","",splnames),
                            splnames)
    ## FIXME: disambiguate?
    splsec_proc <- lapply(splsec,drop_calcs)
    L1 <- L2 <- NULL
    pp <- splsec_proc$PARAMETER
    ## EXPERIMENTAL:
    pp <- pp[!grepl("^ *!!",pp)]
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

read_psv <- function(fn,names=NULL) {
    fn <- tolower(fn) ## arghv
    fn <- paste(fn,"psv",sep=".")
    if (!file.exists(fn)) stop("no PSV file found")
    ans <- read_admbbin(fn)
    if (is.null(names)) names <- paste("V",seq(ncol(ans)),sep="")
    if (length(names)!=ncol(ans)) {
        warning("mismatch between number of columns and number of names")
        names <- c(names,paste("V",seq(length(names)+1,ncol(ans)),sep=""))
    }
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
## FIXME: check bin sizes 
read_admbbin <- function(fn) {
    f <- file(fn,open="rb")
    nv <- readBin(f,"int")
    fs <- file.info(fn)$size
    isize <- 4; dsize <- 8
    m <- matrix(readBin(f,"double",n=(fs-isize)/dsize),byrow=TRUE,
                ncol=nv)
    close(f)
    m
}

## from glmmADMB, by Hans Skaug
write_dat <- "dat_write" <-
    function (name, L, append=FALSE) 
{
    n <- nchar(name)
    file_name <- if (tools::file_ext(name) == ".dat") {
        name
    } else paste(name, "dat", sep = ".")
    cat("# \"", file_name,"\" produced by dat_write() from R2admb ", 
        date(), "\n", file = file_name, sep = "", append=append)
    for (i in 1:length(L)) {
        x <- L[[i]]
        dc <- data.class(x)
        if (dc=="numeric") {
            cat("#", names(L)[i], "\n", L[[i]], "\n\n", file = file_name, 
                append = TRUE)
        } else {
            if (dc == "matrix") {
                cat("#", names(L)[i], "\n", file = file_name, append = TRUE)
                write.table(L[[i]], , col.names = FALSE, row.names = FALSE, quote = FALSE, 
                            file = file_name, append = TRUE)
                cat("\n", file = file_name, append = TRUE)
            } else {
                stop(paste("can't handle data type '",dc,"' (variable ",names(L)[i],")",sep=""))
            }
        }
    }
}

## from glmmADMB, by Hans Skaug
write_pin <- "pin_write" <-
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
            write.table(L[[i]], col.names = FALSE, row.names = FALSE, quote = FALSE, 
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
        tab <- confint.default(object)
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
    tab
}

compile_admb <- function(fn,safe=FALSE,re=FALSE,verbose=FALSE,
                         admb_errors=c("stop","warn","ignore")) {
    admb_errors <- match.arg(admb_errors)
    if (!file.exists(paste(fn,"tpl",sep="."))) stop("can't find TPL file")
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
    ## sorting out the lines that come BEFORE the warnings
    admb_warnerr_index <- grepl("warning|error",coutfile)
    csplit <- split(coutfile,head(c(0,cumsum(admb_warnerr_index)),-1))
    wchunks <- which(sapply(lapply(csplit,grep,pattern="warning"),length)>0)
    echunks <- which(sapply(lapply(csplit,grep,pattern="error"),length)>0)
    if (length(wchunks)>0) {
        if (!verbose) {
            ## figure we don't need these warnings
            ## if we are spitting them out above anyway
            admb_warnings <- paste("from ADMB:",unlist(csplit[wchunks]))
            sapply(admb_warnings,warning)
        }
        csplit <- csplit[-wchunks]
    }
    Sys.chmod(fn,mode="0755")
    if (length(echunks)>0) {
        comperrmsg <- "errors detected in compilation: run with verbose=TRUE to view"
        if (admb_errors=="stop") stop(comperrmsg) else if (admb_errors=="warn") warning(comperrmsg)
    }
}

run_admb <- function(fn,verbose=FALSE,mcmc=FALSE,mcmc.opts=mcmc.control(),
                     profile=FALSE,extra.args="",
                     admb_errors=c("stop","warn","ignore")) {
    admb_errors <- match.arg(admb_errors)
    args <- ""
    if (mcmc) {
        if (is.null(mcmc.opts$mcmcpars)) stop("you must specify at least one parameter in 'mcmc.opts$mcmcpars' (see ?mcmc.control)")
        args <- paste(args,mcmc.args(mcmc.opts))
    }
    if (profile) args <- paste(args,"-lprof")
    if (!missing(extra.args)) {
        args <- paste(args,extra.args)
    }
    if (verbose) cat("running compiled executable with args: '",args,"'...\n")

    outfn <- paste(fn,"out",sep=".")

    if (.Platform$OS.type=="windows") {
        cmdname <- paste(fn,".exe",sep="")
        shellcmd <- shell
    } else {
        cmdname <- paste("./",fn,sep="")
        shellcmd <- system
    }
    if (!file.exists(cmdname)) stop("executable ",cmdname," not found: did you forget to compile it?")
    res <- shellcmd(paste(cmdname,args,">",outfn),intern=TRUE)
    
    outfile <- readLines(paste(fn,".out",sep=""))
    ## replace empty res with <empty> ?
    if (mcmc) {
        ## write MC info to disk so it will be retrievable ...
        mcinfofile <- file(paste(fn,"mcinfo",sep="."),"w")
        mctab <- unlist(mapply(function(x,y) {
            c(paste("# ",x),if (is.null(y))  "" else paste(y,collapse=" "))
        },names(mcmc.opts),mcmc.opts))
        writeLines(mctab,paste(fn,"mcinfo",sep="."))
    }
    if (verbose) {
        cat("Run output:\n",res,"\n",sep="\n")
        cat(outfile,"\n",sep="\n")
    }
    if (length(grep("^Error",outfile)>0)) {
        runerrmsg <- "errors detected in ADMB run: run with verbose=TRUE to view"
        if (admb_errors=="stop") stop(runerrmsg) else if (admb_errors=="warn") warning(runerrmsg)
    }
    invisible(res)
}

read_admb <- function(fn,verbose=FALSE,
                      profile=FALSE,
                      mcmc=FALSE,
                      mcmc.opts=NULL,
                      admbOut=NULL,
                      checkterm=TRUE) {
    tpldat <- read_tpl(fn)  ## extract info from TPL file
    if (verbose) cat("reading output ...\n")
    parfn <- paste(fn,"par",sep=".")
    if (!file.exists(parfn)) stop("couldn't find parameter file ",parfn)
    L <- c(list(fn=fn,txt=admbOut),read_pars(fn))
    if (mcmc) {
        ## if (checkparam!="write") {
        ## warning("MCMC naming is probably wrong")
        ## }
        ## FIXME: get MCMC names -- how?
        if (file.exists(paste(fn,"hst",sep="."))) {
            L <- c(L,list(hist=read_hst(fn)))
        }
        if (is.null(mcmc.opts)) {
            ## try to retrieve mc info from file
            mcinfofile <- paste(fn,"mcinfo",sep=".")
            if (file.exists(mcinfofile)) {
                w <- readLines(mcinfofile)
                wnames <- gsub("^# +","",w[seq(1,length(w),by=2)])
                wvals <- as.list(w[seq(2,length(w),by=2)])
                wvals[c(1,2,5)] <- as.numeric(wvals[c(1,2,5)])
                wvals[3:4] <- as.logical(wvals[3:4])
                wvals[[6]] <- strsplit(wvals[[6]]," ")[[1]]
                names(wvals) <- wnames
                mcmc.opts <- wvals
            } else warning("having difficulty retrieving MCMC info, will try to continue anyway")
        }
        ## if (is.null(mcmc.opts) || is.null(mcmc.opts$mcmcpars) || nchar(mcmc.opts$mcmcpars)==0) {
        ## pnames <- gsub("r_","",tpldat$info$sdnums$vname)
        ## } else pnames <- mcmc.opts$mcmcpars
        sdinfo <- read.table(paste(fn,"std",sep="."),skip=1)
        pnames <- rep_pars(sdinfo[,2])
        ## pnames <- grep("^r_.*",pnames,value=TRUE,invert=TRUE)
        sdreportvars <- as.character(tpldat$info$sdnums$vname)
        pnames <- setdiff(pnames,sdreportvars)
        if (is.null(mcmc.opts) || mcmc.opts[["mcsave"]]>0) {
            L$mcmc <- read_psv(fn,names=pnames)
            ## FIXME: account for mcmc2 if appropriate
            attr(L$mcmc,"mcpar") <- c(1,mcmc.opts[["mcmc"]],mcmc.opts[["mcsave"]])
        }
    }
    if (profile) {
        if (!is.null(tpldat$info$raneff)) {
            stop("something's wrong -- profiling is not implemented for random effects models")
        }
        profpars <- tpldat$info$profparms$vname
        L$prof <- lapply(profpars,read_plt)
        names(L$prof) <- gsub("p_","",profpars)  ## FIXME: maybe dangerous?
    }
    ## compute termination criteria
    ##  can we retrieve hessian directly???
    if (checkterm) {
        v <- with(L,vcov[seq(npar),seq(npar)])
        ev <- try(eigen(solve(v))$value,silent=TRUE)
        L$eratio <- if (inherits(ev,"try-error")) NA else min(ev)/max(ev)
    }
    class(L) <- "admb"
    L
}


