c    Include file for the GA optimization program.
c    This file contains only optimization variables
c    Peter Gerstoft
c
c     CHARACTER(len=6) progna
      CHARACTER*6 progna
      INTEGER prtfil
      PARAMETER(progna='GAINV',prtfil=7)
      INTEGER mpar,mdig,mq,mq_post,mphys, ! max dimension
     &     mlay,mdep,mx,mobs,mopt,
     &     mfreq,msrd,
     &     Meof, Meofvar
c---  mq      is population size  
c---  mq_post is population size, OR the size of all models in post
      PARAMETER(mpar=100,mdig=1200,mq=100,mq_post=200000,mphys=40)  
      PARAMETER(mlay=200,mdep=1100,mx=1500,mobs=1000000,mopt=40,msrd=4)
      PARAMETER(mfreq=1001, Meof=20, Meofvar=65)

c     CHARACTER(len=80) title
      CHARACTER*80 title
      INTEGER iopt(mopt),itrans(mopt),isubopt(mopt),nprior
      REAL iapri(mpar)
      CHARACTER optinv(mopt)
      COMMON /options/iopt,itrans,isubopt,optinv
      INTEGER nparm,q,qin,npopin ! present dimensions 
c     pointers from parm 2 phy and lay
      INTEGER par2phy(mpar),par2lay(mpar),par3(mpar) 
      INTEGER ndigit(mpar),nbit(mpar)     
      
      REAL fmin(mpar),fmax(mpar) ! bounds and possible val
      REAL f(mpar),df(mpar),xstar(mpar),xprior(mpar)
c     CHARACTER(len=40) phystxt(mphys),phystxt2(mphys)
      CHARACTER*40 phystxt(mphys),phystxt2(mphys)
      
      COMMON /optpar/title,nparm,q,iapri,nprior,
     &     par2phy,par2lay,par3,phystxt,phystxt2,
     &     ndigit,nbit,qin,npopin,
     &     fmin,fmax,f,df,xstar,xprior
c     COMPLEX   data(mobs),cov(mobs)
c     EQUIVALENCE (data,cov)
      double precision xcov_scale,xcov_sum,xcov_trace(mfreq),
     &     xcov_lamb(mfreq),zcov_noise(mfreq)
      double precision  xcor_sum(mfreq),xrepl_sum(mfreq),
     &     fit_bart_ind(mfreq)
      integer ncov_siz          ! size of covariance matrix
      integer nsensor           ! number of sensor for each point
      COMMON /optfield/ xcov_scale,xcov_sum,xcov_trace,xcov_lamb,
     &     xcor_sum,xrepl_sum,fit_bart_ind,zcov_noise,
     &     ncov_siz,nsensor
c     INTEGER  model  (mpar,mq_post)  ! integer disc of models
c     INTEGER  allmodel(mpar,mq_post) ! integer disc of models
c     real  allfit(mq_post)      ! corresponding fitness
c     real  allfitbart(mfreq,mq_post)
c     real  allfitbart(10,1000)
c     INTEGER parents(mpar,mq)        ! integer disc of models
      INTEGER iseed,iallmodel,iallstart,ipop
      COMMON /optmod/iseed,
     &     iallmodel,iallstart,ipop
      
      REAL rng1,rng2,drng,temp0,temp1
c     the index for 2D ppd 
      INTEGER ippd1,ippd2,       
     &     ilin(Mpar),ncurv, nummodes,npop,niter,mter
      REAL pm,px,pu
      COMMON /samppar/rng1,rng2,drng,temp0,niter,temp1,
     &     ilin,pm,px,pu,npop,mter,
     &     ippd1,ippd2,ncurv, nummodes
c***  for use of shape functions
      INTEGER Neof, Neofvar, 
     &     par2phy_eof(Meofvar),par2lay_eof(Meofvar),
     &     par3_eof(Meofvar),
     &     point_eof(Meofvar),nlimits
      REAL eofcoef(Meofvar,Meof),aeof(Meof),
     &     fmin_lim(Meofvar),fmax_lim(Meofvar)
      common /eof/Neof,Neofvar,par2phy_eof,par2lay_eof,par3_eof,
     &     eofcoef,aeof,
     &     point_eof,fmin_lim,fmax_lim,nlimits
c**** for masking of data
c     complex    weight(mobs)
c     common /weightscom/ weight
c**** for adding noise to the data
      real snr_db
      common /noise_bl/snr_db
c**** for writing of sound speed profile.
      integer iwrite_soundsp
      common /out/iwrite_soundsp
      
      
      real sigmad,sigmar
      common /sigmas/sigmad,sigmar
      
c***********for Metropolis-Hastings sampling
      integer nppd,iforwt
      real numh,epsstop,epscov,kgrow,rankCd 
      common /mhsamp/numh,epsstop,epscov,kgrow,rankCd, 
     &     nppd,iforwt
