#################################################################################
##
## Author:  Nat Goodman
## Created: 19-09-10
##          from misig/doc_confi.R created 19-07-16
##          from misig/confi.R created 19-07-04
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate figures and tables for confi document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Generate Figures and Tables for effit Blog Post ---
library(nor1mix);
## no sections. only 4 figures
## TODO: turn constants into params
doc_confi=
  function(sect=NULL,
           m=3000,d.obs=0.5,
           mean.true=0.3,mean.false=0,sd.true=0.1,sd.false=0.05,
           n.lo=20,n.hi=200,prop.lo=0.25,prop.hi=0.75,...) {
    mean.mix=c(mean.true,mean.false);
    sd.mix=c(sd.true,sd.false);
    ## draw the figures  
    ## NG 19-09-03: REAL HACK. hardcode xlim,ylim to put all figures on same scale
    xlim=c(-0.2,0.9);
    ylim=c(0,6);
    ## normals - use norMix so same code will work for norm and mixture
    fig_effit('lo_norm',n.lo,1,m,d.obs,mean.true,sd.true,xlim,ylim);
    fig_effit('hi_norm',n.hi,1,m,d.obs,mean.true,sd.true,xlim,ylim);
    ## mixtures
    fig_effit('lo_lo',n.lo,prop.lo,m,d.obs,mean.mix,sd.mix,xlim,ylim);
    fig_effit('hi_lo',n.hi,prop.lo,m,d.obs,mean.mix,sd.mix,xlim,ylim);
    fig_effit('lo_hi',n.lo,prop.hi,m,d.obs,mean.mix,sd.mix,xlim,ylim);
    fig_effit('hi_hi',n.hi,prop.hi,m,d.obs,mean.mix,sd.mix,xlim,ylim);
    ##
    invisible();
  }
## do one figure of effit doc. shows results for one simulation id
fig_effit=function(id,n,prop.true,m,d.obs,mean.mix,sd.mix,xlim,ylim) {
  w=if(prop.true==1) 1 else c(prop.true,1-prop.true);
  mix=norMix(mu=mean.mix,sigma=sd.mix,w=w);
  sim=load_sim(n=n,m=m,d0=d.obs,id=id);
  title=if(prop.true==1) figtitle('Normal prior',n=n,d.obs=d.obs)
        else figtitle('Mixture prior',n=n,prop.true=prop.true,d.obs=d.obs);
   dofig(plot_effit,id,sim=sim,title=title,n=n,d0=d.obs,mix=mix, 
         breaks=25,xlim=xlim,ylim=ylim);
}
## plot histogram, bayesian distributions, medians
## adapted from plothist_dpop
plot_effit=
  function(sim,n,d0,mix,
           title=NULL,cex.title='auto',legend='right',
           xlab='d.pop',ylab='probability density',xlim=NULL,ylim=NULL,
           col.hist='grey90',border.hist='grey80',breaks='Sturges',
           col.distr=RColorBrewer::brewer.pal(3,'Set1'),lwd.distr=2,lty.distr='solid',
           vlty='dashed',vlwd=1,vdigits=2,vcol.dobs='grey50',
           ...){
    distr=cq(bayes,prior);
    ld=length(distr);
    prior=prior_mix(mix);
    median_prior=median_mix(mix);
    sim=sim[sim$n==n,];
    col=setNames(col.distr[1:ld],distr);
    lwd=setNames(rep(lwd.distr,len=ld),distr);
    lty=setNames(rep(lty.distr,len=ld),distr);
    d.pop=sim$d.pop;
    init_bayes(n=n,d0=d0,prior=prior);
    hist.obj=hist(d.pop,breaks=breaks,plot=FALSE);
    if (is.null(xlim)) xlim=range(hist.obj$breaks); # default per R refman
    if (is.null(ylim)) {
      ## compute max y value across hist and all distributions
      ymax=max(hist.obj$density,
               optimize(d_bayes,c(-10,10),maximum=T)$objective,
               optimize(prior,c(-10,10),maximum=T)$objective);
      ylim=c(0,ymax);
    }
    if (is.null(cex.title)|cex.title=='auto') cex.title=cex_title(title);
    ## plot histogram
    plot(hist.obj,freq=F,main=title,cex.main=cex.title,
         col=col.hist,border=border.hist,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,...);
    ## add distributions
    x=seq(xlim[1],xlim[2],len=1000);
    y=cbind(bayes=d_bayes(x),prior=prior(x));
    matlines(x,y,col=col,lty=lty,lwd=lwd);
    ## plot medians and d.obs
    ## TODO: be smarter about label
    vline=c(d.obs=d0,bayes=median_bayes(),prior=median_prior());
    vcol=c(vcol.dobs,col);
    vlab=c(TRUE,FALSE,FALSE);
    vhline(vline=vline,vlab=vlab,vhdigits=vdigits,lty=vlty,col=vcol,lwd=vlwd);
    ## add grid and legend
    grid();
    legend(legend,legend=cq(posterior,prior),title=NULL,col=col,lwd=lwd,lty=lty,cex=0.8,bty='n');
  }
