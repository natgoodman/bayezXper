#################################################################################
##
## Author:  Nat Goodman
## Created: 19-07-16
##          from confi.R created 19-07-04
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate figures and tables for confi document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
doc_confi=function(sect=NULL,id=NULL,...) {
  ## each section displays one simulation id
  sect.all=cq(unif,
              norm_0.0_0.2,norm_0.2_0.2,norm_0.5_0.2,
              norm_0.0_0.4,norm_0.2_0.4,norm_0.5_0.4 );
  if (missing(sect)&!missing(id)) sect=id;
  if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all,start=FALSE);
  sapply(sect,function(sect) {
    sect_start(sect,sect.all);
    sect_confi(id=sect,...);
  });
  return(sect);
}
## do one section of confi doc. shows results for one simulation id
sect_confi=
  function(id='unif',n=param(n.confi),m=param(m.confi),d0=param(d0.confi),tol=1e-3,
           conf.level=param(conf.level),smooth='spline') {
    ## CAUTION: code assumes single d0
    sim=load_sim_all(n,m,d0,id,tol,prune=T);
    sim.byn=split(sim,sim$n);
    cover=do.call(cbind,lapply(n,function(n) {
      sim=sim.byn[[as.character(n)]];
      ci=ci_d2t(n=n,d=d0,simplify=F,conf.level=conf.level);
      lo=ci[1,]; hi=ci[2,];
      cover=do.call(cbind,lapply(sim$d.pop,function(d.pop) between(d.pop,lo,hi)));
      rowSums(cover)/nrow(sim);
      }));
    colnames(cover)=names(sim.byn);
    rownames(cover)=conf.level;
    title=figtitle('coverage vs. confidence level by n',id=id,d0=d_pretty(d0));
    figname=paste(sep='.','cover_by_n',id,
                  paste(sep=',',paste_nv(d0,d_pretty(d0)),paste_nv(m,m_pretty(m))));
    dofig(plotm_cover,figname,cover=cover,title=title,smooth=smooth);
    ## plot 
    
    ## CAUTION: next figures work on single n
    n=pick(n,1);
    sim=sim.byn[[as.character(n)]];
    ## next figures do special things if d.pop prior is normal.  
    ## CAUTION: code splits id to get d.pop prior, mean, sd. short term hack!
    prior=unlist(strsplit(id,'_'));
    distr.prior=prior[1]; mean.prior=as.numeric(prior[2]); sd.prior=as.numeric(prior[3]);
    ## if d.pop prior is normal, do fixed effect meta-analysis.
    ##   from www.johndcook.com/blog/2012/10/29/product-of-normal-pdfs/ comment #2
    ## compute fixed effect params and multi-assign to current environment
    if (distr.prior=='norm')
      list2env(meta(mean.prior,d0,sd.prior,sd_d2t(n=n,d0=d0)),enivronment());
    title=figtitle(c('histogram of d.pop for d.sdz near',d0),id=id,n=n)
    figname=paste(sep='.','hist',id,
                  paste(sep=',',paste_nv(d0,d_pretty(d0)),paste_nv(m,m_pretty(m)),paste_nv(n)));
    dofig(plothist_dpop,figname,sim=sim,n=n,d0=d0,
          distr.prior=distr.prior,mean.meta=mean.meta,sd.meta=sd.meta,
          title=title);
    title=figtitle(c('QQ plot of d.pop vs. d2t for d.sdz near',d0),id=id,n=n)
    figname=paste(sep='.','qq',id,
                  paste(sep=',',paste_nv(d0,d_pretty(d0)),paste_nv(m,m_pretty(m)),paste_nv(n)));
    dofig(plotqq_dpop,figname,sim=sim,n=n,d0=d0,
          distr.prior=distr.prior,mean.meta=mean.meta,sd.meta=sd.meta,
          title=title);
    return();
  }
