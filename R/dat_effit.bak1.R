#################################################################################
##
## Author:  Nat Goodman
## Created: 19-07-16
##          from confi.R created 19-07-04
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate data for confi document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
dat_confi=
  function(n=param(n.confi),m=param(m.confi),d0=param(d0.confi),...) {
    cases=expand.grid(n=n,m=m,d0=d0);
    withrows(cases,case, {
      dosim(n=n,m=m,d0=d0,id='unif',d.gen=runif,d.args=list(min=d0-3,max=d0+3),...);
      dosim(n=n,m=m,d0=d0,id='norm_0.0_0.2',d.gen=rnorm,d.args=list(mean=0.0,sd=0.2),...);
      dosim(n=n,m=m,d0=d0,id='norm_0.2_0.2',d.gen=rnorm,d.args=list(mean=0.2,sd=0.2),...);
      dosim(n=n,m=m,d0=d0,id='norm_0.5_0.2',d.gen=rnorm,d.args=list(mean=0.5,sd=0.2),...);
      dosim(n=n,m=m,d0=d0,id='norm_0.0_0.4',d.gen=rnorm,d.args=list(mean=0.0,sd=0.4),...);
      dosim(n=n,m=m,d0=d0,id='norm_0.2_0.4',d.gen=rnorm,d.args=list(mean=0.2,sd=0.4),...);
      dosim(n=n,m=m,d0=d0,id='norm_0.5_0.4',d.gen=rnorm,d.args=list(mean=0.5,sd=0.4),...);
   });
    invisible();
  }
vrnorm=Vectorize(rnorm,"mean");
dosim=
  function(n,m,d0=0.5,id='unif',d.gen=runif,d.args=list(min=d0-3,max=d0+3),m0=1e4) {
    if (param(verbose)) print(paste(sep=' ',
                                    '>>> dosim:',nvq(id,n,m,d0),format(Sys.time(),"%b %d %X")));;
    more=m; i=1;
    while(more>0) {
      m1=min(m0,more);
      group0=replicate(m1,rnorm(n,mean=0));
      d=do.call(d.gen,c(n=m1,d.args));
      group1=vrnorm(n,mean=d);
      mean0=colMeans(group0);
      mean1=colMeans(group1);
      d.raw=mean1-mean0;
      sd0=apply(group0,2,sd);
      sd1=apply(group1,2,sd);
      sd=pooled_sd(sd0,sd1);
      d.sdz=d.raw/sd;
      sim=data.frame(n,d.pop=d,d.sdz,sd,d.raw,mean0,mean1,sd0,sd1,row.names=NULL);
      save_sim(sim,n,m,d0,id,i);
      more=more-m1; i=i+1;
    }
    ## consolidate subfiles into one
    sim=cat_sim(n,m,d0,id);
    invisible(sim);
  }

## file functions. should be in datman but I'm worried about name conflicts...
filename_sim=function(n,m,d0,id,i=NULL) 
  filename(param(simdir),base=paste(sep='_','sim',id),
           tail=paste(sep=',',paste_nv(n),paste_nv(m,m_pretty(m)),paste_nv(d0,d_pretty(d0))),
           suffix=paste(collapse='.',c(if (!is.null(i)) sprintf("%03i",i),'RData')));

save_sim=function(sim,n,m,d0,id,i=NULL) {
  save(sim,file=filename_sim(n,m,d0,id,i));
}
load_sim=get_sim=function(n,m,d0=0.5,id='unif',i=NULL,subfiles=FALSE) {
  if ((is.null(i)&!subfiles)|!is.null(i)) sim=load_(filename_sim(n,m,d0,id,i),'sim')
  ## if (!is.null(i)) sim=load_(filename_sim(n,m,d0,id,i),'sim')
  else {
     pattern=paste0('sim_',id,'\\.',
                    paste(sep=',',paste_nv(n),paste_nv(m,m_pretty(m)),paste_nv(d0,d_pretty(d0))),
                    paste0('\\.','\\d+','\\.','\\RData'));
    files=
      list.files('data/confi/sim',full.names=T,pattern=pattern);
    if (length(files)==0) stop(paste0('no files found: pattern=',pattern));
    sim=do.call(rbind,lapply(files,function(file) load_(file,'sim')))
  }
  invisible(sim);
}
## remove sim subfiles after consolidating
rm_sim=function(n,m,d0=0.5,id='unif',i=NULL,subfiles=T) {
  if ((is.null(i)&!subfiles)|!is.null(i)) file.remove(filename_sim(n,m,d0,id,i))
       pattern=paste0('sim_',id,'\\.',
                    paste(sep=',',paste_nv(n),paste_nv(m,m_pretty(m)),paste_nv(d0,d_pretty(d0))),
                    paste0('\\.','\\d+','\\.','\\RData'));
  files=
      list.files('data/confi/sim',full.names=T,pattern=pattern);
  if (length(files)==0) stop(paste0('no files found: pattern=',pattern));
  file.remove(files);
}
## consolidate sim subfiles, save as one file, rm subfiles
cat_sim=function(n,m,d0=0.5,id='unif') {
  sim=load_sim(n,m,d0,id,subfiles=T);
  save_sim(sim,n,m,d0,id);
  rm_sim(n,m,d0,id);
  invisible(sim);
}
## load all (consolidated) sim files for parameter grid, optionally pruning to region near d0
load_sim_all=get_sim_all=
  function(n=param(n.confi),m=param(m.confi),d0=param(d0.confi),id=param(id.confi),
           tol=param(tol),prune=T) {
    cases=expand.grid(n=n,m=m,d0=d0);
    sim=do.call(rbind,withrows(cases,case,{
      sim=load_sim(n=n,m=m,d0=d0,id=id);
      if (prune) sim=subset(sim,subset=near(d.sdz,d0,tol));
    }));
    invisible(sim);
  }
