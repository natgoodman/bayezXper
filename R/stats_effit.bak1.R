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
## return list of meta mean,sd
meta_normd2t=function(mean.prior,sd.prior,d0,n) 
    meta(mean.prior,sd.prior,d0,sd_d2t(n=n,d0=d0))
meta=function(m1,sd1,m2,sd2) {
  w1=1/(sd1^2); w2=1/(sd2^2);
  ww=w1+w2;
  m12=(w1*m1+w2*m2)/(w1+w2);
  var12=1/ww;
  sd12=sqrt(var12);
  list(mean=m12,sd=sd12);
}
## return function for posterior probability of d.pop given d0 (d.sdz)
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
    d.bayes=meta_bayes(n,d,mean.prior,sd.prior);
    d.interp=approxfun(dx,d.bayes(dx),yleft=0,yright=0);
    p.bayes=Vectorize(function(d0) integrate(d.interp,-10,d0)$value);
    ci=ci_bayes_(p.bayes,conf.level);
    if (simplify&ncol(ci)==1) ci=as.vector(ci);
    ci;
  }
ci_bayes_=Vectorize(function(p.bayes,conf.level) {
  p0=(1-conf.level)/2; p1=1-p0;
  lo=suppressWarnings(
    uniroot(function(d0) p.bayes(d0)-p0,interval=c(-10,10))$root);
  hi=suppressWarnings(
    uniroot(function(d0) p.bayes(d0)-p1,interval=c(-10,10))$root);
  c(lo,hi);
},vectorize.args='conf.level')

##### probability functions of meta_bayes
##### I'm pretty sure d_bayes & mean_bayes work. others not well-tested
## probability density of meta_bayes
d_bayes=function(x,n,d0,mean.prior,sd.prior) {
  d.bayes=meta_bayes(n,d0,mean.prior,sd.prior);
  d.bayes(x)
}
mean_bayes=function(n,d0,mean.prior,sd.prior) 
  integrate(function(x) x*d_bayes(x,n,d0,mean.prior,sd.prior),-Inf,Inf)$value;
  
## cumulative probability of meta_bayes
p_bayes=function(x,n,d0,mean.prior,sd.prior,dlim=c(-2,2),dinc=.005) {
  dx=seq(min(dlim),max(dlim),by=dinc);
  d.bayes=meta_bayes(n,d,mean.prior,sd.prior);
  d.interp=approxfun(dx,d.bayes(dx),yleft=0,yright=0);
  p.bayes=Vectorize(function(x) integrate(d.interp,-10,x)$value);
  p.bayes(x);
}
## quantile function of meta_bayes
q_bayes=function(p,n,d0,mean.prior,sd.prior,dlim=c(-2,2),dinc=.005) {
  dx=seq(min(dlim),max(dlim),by=dinc);
  d.bayes=meta_bayes(n,d0,mean.prior,sd.prior);
  d.interp=approxfun(dx,d.bayes(dx),yleft=0,yright=0);
  p.bayes=Vectorize(function(x) integrate(d.interp,-10,x)$value);
  sapply(p,function(p) uniroot(function(x) p.bayes(x)-p,interval=c(-10,10))$root);
}
## random function of meta_bayes
r_bayes=function(m,n,d0,mean.prior,sd.prior,qlim=c(1e-2,1-1e-2),qinc=1e-2) {
  p=runif(m,min(qlim),max(qlim));
  if (m>1/qinc) {
    dq=seq(min(qlim),max(qlim),by=qinc);
    q.interp=approxfun(dq,q_bayes(dq,n,d0,mean.prior,sd.prior),rule=2);
    q.interp(p);
  } else q_bayes(p,n,d0,mean.prior,sd.prior);
}
####################
## quick hack to compute ci for d2t posterior probability of d.pop
## I think this is what unif is
## adapted from ci_norm
ci_xxx=function(n,d0,simplify=T,conf.level=param(conf.level)) {
  ci=ci_xxx_(n,d0,conf.level);
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}
ci_xxx_=Vectorize(function(n,d0,conf.level) {
  p0=(1-conf.level)/2; p1=1-p0;
  lo=q_d2t(n,p0,d0);
  hi=q_d2t(n,p1,d0);
  c(lo,hi);
},vectorize.args='conf.level')


####################
## generate standard line properties
## n to color
n2col=function(n=param(n.confi)) setNames(RColorBrewer::brewer.pal(max(3,length(n)),'Set1'),n)
## distr to line type. by default, sim is 'dotted', rest are 'solid','dashed', 'dotdash', ...
distr2lty=
  function(distribution,lty.sim='dotted',lty.all=cq(solid,dashed,dotdash,longdash,twodash)) {
  lty=NULL;
  if ('sim' %in% distribution) {
    lty=c(sim=lty.sim);
    distribution=setdiff(distribution,'sim');
  }
  ld=length(distribution);
  if (ld>length(lty.all)) stop('more distributions than line types. is distribution right?');
  c(lty,setNames(lty.all[1:ld],distribution));
}
## distr to line width. by default, sim is 1.5, rest are 2
distr2lwd=function(distribution,lwd.sim=1.5,lwd.all=c(2)) {
  lwd=NULL;
  if ('sim' %in% distribution) {
    lwd=c(sim=lwd.sim);
    distribution=setdiff(distribution,'sim');
  }
  lwd.all=rep(lwd.all,len=length(distribution));
  c(lwd,setNames(lwd.all,distribution));
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

########################################
## below here from stats.R
## some of these functions may get moved back when suitably generalized
## ---- Statistical Functions for Simulated Distribution ----
## NG 19-07-22: VERY INCOMPLETE! just realized I needed this. sigh
##   quick web search suggests this is more complicated than it looks...
## data is vector 
d_sim=function(data) stop('d_sim not yet implemented')
p_sim=function(data,d,...) ecdf(data,...)(d);
q_sim=function(data,p,...) quantile(data,probs=p,...)
r_sim=function(n,data) q_sim(data,runif(n,0,1))

## this one is equivalent to the one below
## ci_sim_=function(data,d,conf.level) {
##   p0=(1-conf.level)/2; p1=1-p0;
##   p_sim=ecdf(data);
##   lo=suppressWarnings(uniroot(function(d0) p_sim(d0)-p0,interval=c(-10,10))$root);
##   hi=suppressWarnings(uniroot(function(d0) p_sim(d0)-p1,interval=c(-10,10))$root);
##   c(lo,hi);
## }
## ci_sim=Vectorize(ci_sim_,vectorize.args=c('d','conf.level'));

## adapted from ci_norm
ci_simNORM=function(data,d,conf.level) {
  p0=(1-conf.level)/2; p1=1-p0;
  lo=quantile(data,probs=p0,type=1);
  hi=quantile(data,probs=p1,type=1);
  ci=rbind(lo,hi);
  colnames(ci)=conf.level;
  ci;
}
## this on is from first principles
ci_sim=function(data,d,conf.level) {
  p0=1-conf.level; p1=1-p0;
  data.lt=data[data<=d];
  data.gt=data[data>=d];
  lo=quantile(data.lt,probs=p0,type=1);
  hi=quantile(data.gt,probs=p1,type=1);
  ci=rbind(lo,hi);
  colnames(ci)=conf.level;
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

## adapted from https://www.cyclismo.org/tutorial/R/confidence.html
ci_norm=function(mean,sd,d,conf.level) {
  p0=(1-conf.level)/2; p1=1-p0;
  lo=qnorm(p0,mean=mean,sd=sd);
  hi=qnorm(p1,mean=mean,sd=sd);
  ci=rbind(lo,hi);
  colnames(ci)=conf.level;
  ci;
}
