ReedfrogSizepred <- 
  data.frame(TBL = rep(c(9,12,21,25,37),each=3),
             Kill = c(0,2,1,3,4,5,0,0,0,0,1,0,0,0,0L))

library(R2admb)
setup_admb()
m1 <- do_admb("ReedfrogSizepred0",
              data=c(list(nobs=nrow(ReedfrogSizepred),
                nexposed=rep(10,nrow(ReedfrogSizepred))),
                ReedfrogSizepred),
              params=list(c=0.45,d=13,g=1),
              run.opts=run.control(checkparam="write",
                checkdata="write"),
              verbose=TRUE)
m1P <- do_admb("ReedfrogSizepred0",
              data=c(list(nobs=nrow(ReedfrogSizepred),
                nexposed=rep(10,nrow(ReedfrogSizepred))),
                ReedfrogSizepred),
               params=list(c=0.45,d=13,g=1),
               profile=TRUE,
               profpars=c("c","d","g"),
               run.opts=run.control(checkparam="write",
                 checkdata="write"),
               admb_errors="warn",  ## because of Hessian problem
               verbose=TRUE)
m1MC <- do_admb("ReedfrogSizepred0",
              data=c(list(nobs=nrow(ReedfrogSizepred),
                nexposed=rep(10,nrow(ReedfrogSizepred))),
                ReedfrogSizepred),
                params=list(c=0.45,d=13,g=1),
                run.opts=run.control(checkparam="write",
                  checkdata="write"),
                mcmc=TRUE,
                mcmc.opts=mcmc.control(mcmcpars=c("c","d","g")),
                verbose=TRUE)
unlink("reedfrogsizepred0")
save("m1","m1MC","m1P",file="Reedfrog_runs.RData")

