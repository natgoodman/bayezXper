#################################################################################
##
## Author:  Nat Goodman
## Created: 19-08-04
##          from doc_confi.R created 19-07-16
##          from confi.R created 19-07-04
##
## Copyright (C) 2019 Nat Goodman.
## 
## Specialized data generation functions for confi document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## compute ci stats (lo,hi,coverage) from sim
doci_confi=function(sim,n,d0,
                    ci.method=cq(std,sim,simq,d2tpost,normpost,meta,bayes),meta,prior,
                    conf.level=param(conf.level)) {
  if (missing(ci.method)) ci.method='std'
  else ci.method=match.arg(ci.method,several.ok=TRUE);
  if (('meta' %in% ci.method)&&missing(meta))
    stop("ci.method contains 'meta' but meta argument missing");
  if (('bayes' %in% ci.method)&&missing(prior))
    stop("ci.method contains 'bayes' but prior argument missing");
  sim.byn=split(sim,sim$n);
  ci=do.call(rbind,lapply(n,function(n) {
    sim=sim.byn[[as.character(n)]];
    d.pop=sim$d.pop;
    m=length(d.pop);
    init_bayes(n=n,d0=d0,prior=prior);
    ## ci=matrix(nrow=length(conf.level),ncol=0);
    ci=do.call(cbind,lapply(ci.method,function(method) {
      switch(method,
             std=doci1(ci_d2t,list(n=n,d=d0),d.pop,conf.level),
             sim=doci1(ci_sim,list(data=d.pop,d0=d0),d.pop,conf.level),
             simq=doci1(ci_simq,list(data=d.pop),d.pop,conf.level),
             d2tpost=doci1(ci_d2tpost,list(n=n,d0=d0),d.pop,conf.level),
             normpost=doci1(ci_norm,list(mean=d0,sd=sqrt(2/n)),d.pop,conf.level),
             meta=doci1(ci_norm,list(mean=meta$mean,sd=meta$sd),d.pop,conf.level),
             bayes=doci1(ci_bayes,list(),d.pop,conf.level),
             stop(paste0('unknown method: ',method,'. should have been caught earlier')));
    }));
    colnames(ci)=paste(sep='.',rep(ci.method,each=3),cq(lo,hi,cvr));
    data.frame(n=n,conf.level,ci);
  }));
  ci;
}
doci1=function(ci_fun,ci.args,d.pop,conf.level) {
  ci.args=c(ci.args,list(simplify=F,conf.level=conf.level));
  ci=do.call(ci_fun,ci.args);
  lo=ci[1,]; hi=ci[2,];
  cvr=do.call(rbind,lapply(d.pop,function(d.pop) between(d.pop,lo,hi)));
  cvr=colSums(cvr)/length(d.pop);
  cbind(lo,hi,cvr);
}
