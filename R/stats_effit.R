#################################################################################
##
## Author:  Nat Goodman
## Created: 19-08-04
##          from doc_confi.R created 19-07-16 & bayes_confi.R created 19-07-21
##            with additional content from stats.R
##          doc_confi.R created from confi.R created 19-07-04
##
## Copyright (C) 2019 Nat Goodman.
## 
## Specialized stats functions for confi document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## create functions for bayesian posterior probability of d.pop given d0 (d.sdz)
##   and assigns into calling environment
##   also assigns bayes list containing these functions
## like R's built-in distributions, imports functions for
##   density, cumulative probability, quantiles, random with prefixes 'd', 'p', 'q', 'r'
## uses suffixes '_exact', -_interp' for exact and interpolated versions, and
##   '_bayes' for whichever is selected as default
## n, d0 are sample size and given effect size (d.sdz)
## prior is prior density of d.pop - must be function
## bayes.function controls which functions are  assigned to default '_bayes' notation
## use.interp controls whether interp functions used in creating next tier of exact functions
## dlim, dinc are range of d values used for interpolation
init_bayes=
  function(n,d0,prior,bayes.function=cq(interp,exact),use.interp=T,
           dlim=c(d0-3,d0+3),dinc=.005,qlim=c(0.005,1-0.005),qinc=.005) {
    if (!is.function(prior)) stop("'prior' must be function");
    bayes.function=match.arg(bayes.function);
    dx=seq(min(dlim),max(dlim),by=dinc);
    qx=seq(min(qlim),max(qlim),by=qinc);
    d_exact=post_bayes(n,d0,prior);
    d_interp=dbayes_interp(d_exact,dx);
    p_exact=pbayes_exact(if(use.interp) d_interp else d_exact);
    p_interp=pbayes_interp(p_exact,dx);
    q_exact=qbayes_exact(if(use.interp) p_interp else p_exact);
    q_interp=qbayes_interp(q_exact,qx);
    r_exact=rbayes_exact(if(use.interp) q_interp else q_exact);
    ## interp makes no sense for 'r'. use exact so r_bayes set correctly
    r_interp=r_exact;
    ## assign function to bayes list and parent environment
    bayes=do.call(c,lapply(cq(d,p,q,r),function(letter) {
      do.call(c,lapply(cq(exact,interp),function(what) {
        bayes=list();
        name=paste(sep='_',letter,what);
        value=get(name);
        bayes[[name]]=value;
        if (bayes.function==what) bayes[[paste(sep='_',letter,'bayes')]]=value;
        bayes;
      }));
    }));
    ## assign functions to global
    list2env(bayes,envir=.GlobalEnv);
    ## assign bayes to parent
    assign('bayes',bayes,envir=parent.frame(n=1));
    invisible(bayes);
  }
## return function for posterior probability of d.pop given d0 (d.sdz)
## prior is function of d.pop, eg, function(d.pop) dnorm(d.pop,mean=mean.prior,sd=sd.prior)
post_bayes=function(n,d0,prior) {
  ## prior probability of d.pop
  pA=prior;
  ## likelihood - probability of d0 (d.sdz) given d.pop
  pBgivenA=function(d.pop) d_d2t(n=n,d0=d.pop,d=d0);
  ## joint prob density
  pjoint=Vectorize(function(d.pop) pA(d.pop)*pBgivenA(d.pop));
  ## prob of d.pop given d0 (d.sdz)
  pAgivenB=pjoint;
  ## average liklihood - denominator in Bayes formula
  pB=integrate(function(d.pop) pAgivenB(d.pop),-Inf,Inf)$value;
  ## final function
  ## postAgivenB=Vectorize(function(d.pop) pAgivenB(d.pop)/pB);
  postAgivenB=function(d.pop) pAgivenB(d.pop)/pB;
}
## return probability functions for d.pop.
## unif
prior_unif=dprior_unif=function(d0,u) {
  span=u/2;
  lim=d0+c(-span,span);
  ## function(d.pop) ifelse((d.pop>=lim[1])&(d.pop<=lim[2]),1/u,0);
  function(d.pop) dunif(d.pop,lim[1],lim[2]);
}
pprior_unif=function(d0,u) {
  span=u/2;
  lim=d0+c(-span,span);
  ## function(d.pop) ifelse((d.pop>=lim[1])&(d.pop<=lim[2]),1/u,0);
  function(d.pop) punif(d.pop,lim[1],lim[2]);
}
qprior_unif=function(d0,u) {
  span=u/2;
  lim=d0+c(-span,span);
  function(p) qunif(p,lim[1],lim[2]);
}
rprior_unif=function(d0,u) {
  span=u/2;
  lim=d0+c(-span,span);
  function(m) runif(m,lim[1],lim[2]);
}
## norm
prior_norm=dprior_norm=
  function(mean.prior,sd.prior) function(d.pop) dnorm(d.pop,mean=mean.prior,sd=sd.prior);
pprior_norm=
  function(mean.prior,sd.prior) function(d.pop) pnorm(d.pop,mean=mean.prior,sd=sd.prior);
qprior_norm=
  function(mean.prior,sd.prior) function(p) qnorm(p,mean=mean.prior,sd=sd.prior);
rprior_norm=
  function(mean.prior,sd.prior) function(m) rnorm(m,mean=mean.prior,sd=sd.prior);

## mix - via nor1mix
prior_mix=dprior_mix=function(mix) function(d.pop) dnorMix(d.pop,mix)
pprior_mix=function(mix) function(d.pop) pnorMix(d.pop,mix)
qprior_mix=function(mix) function(p) qnorMix(p,mix)
rprior_mix=function(mix) function(m) rnorMix(m,mix)
median_mix=function(mix) function() qprior_mix(mix)(0.5)

## return probability functions for bayes
## NG 19-08-31: stop.on.error=F needed for 'mix' prior. dunno why...
pbayes_exact=function(d.bayes) Vectorize(function(x)
  integrate(d.bayes,-10,x,stop.on.error=FALSE)$value);
qbayes_exact=function(p.bayes) Vectorize(function(p) {
  if (p>0&p<1) uniroot(function(x) p.bayes(x)-p,interval=c(-10,10))$root
  else if (p<0|p>1) NaN else if (p==0) -Inf else Inf;
});
rbayes_exact=function(q.bayes) function(m) q.bayes(runif(m));

dbayes_interp=function(d_exact,dx) approxfun(dx,d_exact(dx),yleft=0,yright=0);
pbayes_interp=function(p_exact,dx) approxfun(dx,p_exact(dx),yleft=0,yright=1);
qbayes_interp=function(q_exact,qx) approxfun(qx,q_exact(qx),rule=2);

## these assume bayes already initialized
mean_bayes=function() integrate(function(x) x*d_bayes(x),-Inf,Inf)$value;
median_bayes=function() q_bayes(0.5);
  
var_bayes=function() {
  mu=mean_bayes();
  integrate(function(x) ((x-mu*x)^2)*d_bayes(x),-Inf,Inf)$value;
}
sd_bayes=function() sqrt(var_bayes());
ci_bayes=function(simplify=T,conf.level=param(conf.level)) {
  ci=ci_bayes_(conf.level);
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}
ci_bayes_=Vectorize(function(conf.level) {
  p0=(1-conf.level)/2; p1=1-p0;
  lo=suppressWarnings(
    uniroot(function(d0) p_bayes(d0)-p0,interval=c(-10,10))$root);
  hi=suppressWarnings(
    uniroot(function(d0) p_bayes(d0)-p1,interval=c(-10,10))$root);
  c(lo,hi);
},vectorize.args='conf.level')
                              
########################################
## some of these functions may get moved to stats.R when suitably generalized
## empirical distribution from sim
## NG 19-07-22: VERY INCOMPLETE! just realized I needed this. sigh
##   quick web search suggests this is more complicated than it looks...
## data is vector 
d_sim=function(data) stop('d_sim not yet implemented')
p_sim=function(data,d,...) ecdf(data,...)(d);
q_sim=function(data,p,...) quantile(data,probs=p,...)
r_sim=function(m,data) as.vector(q_sim(data,runif(m,0,1)))

## this one is equivalent to the one below
## ci_sim_=function(data,d,conf.level) {
##   p0=(1-conf.level)/2; p1=1-p0;
##   p_sim=ecdf(data);
##   lo=suppressWarnings(uniroot(function(d0) p_sim(d0)-p0,interval=c(-10,10))$root);
##   hi=suppressWarnings(uniroot(function(d0) p_sim(d0)-p1,interval=c(-10,10))$root);
##   c(lo,hi);
## }
## ci_sim=Vectorize(ci_sim_,vectorize.args=c('d','conf.level'));

## this one adapted from ci_norm. suffix 'q' because uses q_sim
ci_simq=function(data,simplify=T,conf.level=param(conf.level)) {
  p0=(1-conf.level)/2; p1=1-p0;
  lo=quantile(data,probs=p0,type=1);
  hi=quantile(data,probs=p1,type=1);
  ci=rbind(lo,hi);
  colnames(ci)=conf.level;
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}
## this one is from first principles
ci_sim=function(data,d0,simplify=T,conf.level=param(conf.level)) {
  p0=1-conf.level; p1=1-p0;
  data.lt=data[data<=d0];
  data.gt=data[data>=d0];
  lo=quantile(data.lt,probs=p0,type=1);
  hi=quantile(data.gt,probs=p1,type=1);
  ci=rbind(lo,hi);
  colnames(ci)=conf.level;
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}
## this one is equivalent to the one below
## ci_norm_=function(mean,sd,d,conf.level) {
##   p0=(1-conf.level)/2; p1=1-p0;
##   p_norm=function(d0) pnorm(d0,mean,sd);
##   lo=suppressWarnings(uniroot(function(d0) p_norm(d0)-p0,interval=c(-10,10))$root);
##   hi=suppressWarnings(uniroot(function(d0) p_norm(d0)-p1,interval=c(-10,10))$root);
##   c(lo,hi);
## }
## ci_norm=Vectorize(ci_norm_,vectorize.args=c('d','conf.level'));

## normal
## adapted from https://www.cyclismo.org/tutorial/R/confidence.html
ci_norm=function(mean,sd,simplify=T,conf.level=param(conf.level)) {
  p0=(1-conf.level)/2; p1=1-p0;
  lo=qnorm(p0,mean=mean,sd=sd);
  hi=qnorm(p1,mean=mean,sd=sd);
  ci=rbind(lo,hi);
  ## colnames(ci)=conf.level;
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}
## d2t posterior probability of d.pop - I think this is what unif is
## adapted from ci_norm
ci_d2tpost=function(n,d0,simplify=T,conf.level=param(conf.level)) {
  ci=ci_d2tpost_(n,d0,conf.level);
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}
ci_d2tpost_=Vectorize(function(n,d0,conf.level) {
  p0=(1-conf.level)/2; p1=1-p0;
  lo=q_d2t(n,p0,d0);
  hi=q_d2t(n,p1,d0);
  c(lo,hi);
})
#################### STILL USED!! ####################
## return list of meta mean,sd. for normal prior
meta_normd2t=function(mean.prior,sd.prior,d0,n) 
    meta_genl(mean.prior,sd.prior,d0,sd_d2t(n=n,d0=d0))
meta_genl=function(m1,sd1,m2,sd2) {
  w1=1/(sd1^2); w2=1/(sd2^2);
  ww=w1+w2;
  m12=(w1*m1+w2*m2)/(w1+w2);
  var12=1/ww;
  sd12=sqrt(var12);
  list(mean=m12,sd=sd12);
}

