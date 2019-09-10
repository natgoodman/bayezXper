#################################################################################
##
## Author:  Nat Goodman
## Created: 19-07-21
##
## Copyright (C) 2019 Nat Goodman.
## 
## Bayeian analysis for confi. Extracted from shell transcript
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## first block NOT USED in running code. Just a record of what I did
##
## prior probability of d.pop
pA=function(d.pop) dnorm(d.pop,mean=0.2,sd=0.2)
## likelihood - probability pf d.sdz given d.pop
pBgivenA=function(d.pop,d.sdz) d_d2t(n=N,d0=d.pop,d=d.sdz)
## joint prob density for fixed mean.prior, sd.prior, n
pjoint=function(d.pop,d.sdz) pA(d.pop)*pBgivenA(d.pop,d.sdz)
## another way to express prob of d.pop given d.sdz. returns function
pAgivenB=function(d.sdz) {f=function(d.pop) pjoint(d.pop,d.sdz); Vectorize(f)}
## joint prob specialized for d.sdz=0.5 & another way of saying the same thing
pAgiven0.5=function(d.pop) pjoint(d.pop,d.sdz=0.5)
pAgiven0.5=pAgivenB(0.5)
## average liklihood is denominator in Bayes formula, aka pB
## for reasons I don't understand, it's not a probability even though I get it by
##   integrating what I thought were probability density functions
pB=function(d.sdz) integrate(pAgivenB(d.sdz),-Inf,Inf)$value
## final function
postAgivenB=function(d.sdz) {
  avelik=pB(d.sdz);
  function(d.pop) pAgivenB(d.sdz)(d.pop)/avelik;
}

##########
## quick try for 50:50 mixture of two normals
pA=function(d.pop) 0.5*(dnorm(d.pop,mean=0.0,sd=0.2)+dnorm(d.pop,mean=0.5,sd=0.2))
## rest are as above

## test
chk_norm2=function() {
  sim0=load_sim_all(n=N,m,d0,id='norm_0.0_0.2',tol,prune=T);
  sim1=load_sim_all(n=N,m,d0,id='norm_0.5_0.2',tol,prune=T);
  sim=rbind(sim0,sim1)
  d.pop=sim$d.pop
  hist.obj=hist(d.pop,plot=FALSE);
  ylim=c(0,2.5)                           # quick guess after seeing plot
  xlim=range(hist.obj$breaks)
  dev.new(); 
  plot(hist.obj,freq=F,col='grey90',border='grey80',ylim=ylim)
  d.pop=seq(xlim[1],xlim[2],by=.01)
  lines(d.pop,postAgivenB(0.5)(d.pop),col='green')
}

##########
## version that could be used in running code
## return function for posterior probability of d.pop given d.sdz
##
## this version direct adaptation of first code block
meta_bayes_bak=function(n,d.sdz,mean.prior,sd.prior) {
  ## prior probability of d.pop
  pA=function(d.pop) dnorm(d.pop,mean=mean.prior,sd=sd.prior);
  ## likelihood - probability of d.sdz given d.pop
  pBgivenA=function(d.pop,d.sdz) d_d2t(n=n,d0=d.pop,d=d.sdz);
  ## joint prob density
  pjoint=function(d.pop,d.sdz) pA(d.pop)*pBgivenA(d.pop,d.sdz);
  ## prob of d.pop given d.sdz. 
  pAgivenB=function(d.sdz) {f=function(d.pop) pjoint(d.pop,d.sdz); Vectorize(f)}
  ## average liklihood - denominator in Bayes formula
  pB=function(d.sdz) {
    f=pAgivenB(d.sdz);
    integrate(function(d.pop) f(d.pop),-Inf,Inf)$value;
  }
  ## final function
  postAgivenB=function(d.sdz) {
    avelik=pB(d.sdz);
    function(d.pop) pAgivenB(d.sdz)(d.pop)/avelik;
  }
  postAgivenB;
}
## this version simplified - used in running code
meta_bayes=function(n,d0,mean.prior,sd.prior) {
  ## prior probability of d.pop
  pA=function(d.pop) dnorm(d.pop,mean=mean.prior,sd=sd.prior);
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
  postAgivenB;
}

ci_bayes=
  function(n,d,mean.prior,sd.prior,dlim=c(-2,2),dinc=.005,simplify=T,
           conf.level=param(conf.level)) {
    dx=seq(min(dlim),max(dlim),by=dinc);
    d_bayes=meta_bayes(n,d,mean.prior,sd.prior);
    d_interp=approxfun(dx,d_bayes(dx),yleft=0,yright=0);
    p_bayes=Vectorize(function(d0) integrate(d_interp,-10,d0)$value);
    ci=ci_bayes_(p_bayes,conf.level);
    if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}

  
ci_bayes_=function(n,d,mean.prior,sd.prior,conf.level) {
  d_bayes=meta_bayes(n,d,mean.prior,sd.prior);
  p_bayes=function(d0) integrate(d_bayes,-10,d0)$value;
  p0=(1-conf.level)/2; p1=1-p0;
  lo=suppressWarnings(
    uniroot(function(d0) p_bayes(d0)-p0,interval=c(-10,10))$root);
  hi=suppressWarnings(
    uniroot(function(d0) p_bayes(d0)-p1,interval=c(-10,10))$root);
  c(lo,hi);
}

##########
## derive posterior for unif
## u is width (ie, max-min). so prior is 1/u
## return function for posterior probability of d.pop given d0 (d.sdz)
meta_bayes_unif=function(n,d0,u) {
  ## prior probability of d.pop
  pA=function(d.pop) 1/u
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
  postAgivenB;
}
## quantile function of meta_bayes_unif
q_bayes_unif=function(p,n,d0,u,dlim=c(d0-(u/2),d0=(u/2)),dinc=.005) {
  dx=seq(min(dlim),max(dlim),by=dinc);
  d.bayes=meta_bayes_unif(n,d0,u);
  d.interp=approxfun(dx,d.bayes(dx),yleft=0,yright=0);
  p.bayes=Vectorize(function(x) integrate(d.interp,dlim[1],x)$value);
  sapply(p,function(p) uniroot(function(x) p.bayes(x)-p,interval=c(-10,10))$root);
}
## quick hack to compute ci for bayes_unif posterior probability of d.pop
## adapted from ci_xxx
ci_uuu=function(n,d0,u,simplify=T,conf.level=param(conf.level),
                dlim=c(d0-(u/2),d0=(u/2)),dinc=.005) {
  dx=seq(min(dlim),max(dlim),by=dinc);
  ci=ci_uuu_(n,d0,u,conf.level,dx);
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}
ci_uuu_=Vectorize(function(n,d0,u,conf.level,dx) {
  p0=(1-conf.level)/2; p1=1-p0;
  d.bayes=meta_bayes_unif(n,d0=d0,u);
  d.interp=approxfun(dx,d.bayes(dx),yleft=0,yright=0);
  p.bayes=Vectorize(function(x) integrate(d.interp,dlim[1],x)$value);
  lo=uniroot(function(x) p.bayes(x)-p0,interval=dlim)$root;
  hi=uniroot(function(x) p.bayes(x)-p1,interval=dlim)$root;
  c(lo,hi);
},vectorize.args='conf.level')

## mean and sd of meta_bayes_unif
mean_bayes_unif=function(n,d0,u) {
  d.bayes=meta_bayes_unif(n,d0,u);
  integrate(function(x) x*d.bayes(x),-Inf,Inf)$value;
}
var_bayes_unif=function(n,d0,u) {
  d.bayes=meta_bayes_unif(n,d0,u);
  mu=integrate(function(x) x*d.bayes(x),-Inf,Inf)$value;
  integrate(function(x) ((x-mu*x)^2)*d.bayes(x),-Inf,Inf)$value;
}
sd_bayes_unif=function(n,d0,u) sqrt(var_bayes_unif(n,d0,u))


