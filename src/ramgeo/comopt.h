c    Include file for the GA optimization program.
c    This file contains only optimization variables
c    Peter Gerstoft
c
      CHARACTER*6 progna
      INTEGER prtfil
      PARAMETER(progna='GAINV',prtfil=7)
      INTEGER   mpar,mdig,mq,mq_post,mphys,           ! max dimension
     &          mlay,mdep,mx,mobs,mopt,
     &          mfreq,msrd,
     &          Meof, Meofvar
c---  mq          is population size  
c---  mq_post      is population size, OR the size of all models read in post
      PARAMETER (mpar=20,mdig=1200,mq=100,mq_post=100000,mphys=40)  
      PARAMETER(mlay=40,mdep=1100,mx=1000,mobs=100000,mopt=40,msrd=4)
      PARAMETER(mfreq=1001, Meof=20, Meofvar=40)

      CHARACTER*80 title
      INTEGER iopt(mopt),itrans(mopt),isubopt(mopt),iapri(mpar),nprior
      COMMON /options/iopt,itrans,isubopt
      INTEGER    nparm,q            ! present dimensions 
      INTEGER    par2phy(mpar),par2lay(mpar),par3(mpar) ! pointers from parm 2 phy and lay
      INTEGER    ndigit(mpar),nbit(mpar)     
c      REAL       fmin(mpar),fmax(mpar),fval(0:mdig,mpar) ! bounds and possible val
      REAL       fmin(mpar),fmax(mpar)      ! bounds and possible val
      REAL       f(mpar),df(mpar),xstar(mpar)
      CHARACTER*40 phystxt(mphys),phystxt2(mphys)
      
      COMMON /optpar/ title,nparm,q,iapri,nprior,
     &                par2phy,par2lay,par3,phystxt,phystxt2,
     &                ndigit,nbit,
     &                fmin,fmax,f,df,xstar
      COMPLEX   data(mobs),cov(mobs)
      EQUIVALENCE (data,cov)
      real*8 xcov_scale,xcov_sum,xcov_trace(mfreq),
     &           xcov_lamb(mfreq),zcov_noise(mfreq)
      real*8 xcor_sum(mfreq),xrepl_sum(mfreq),fit_bart_ind(mfreq)
      integer ncov_siz                  ! size of covariance matrix
      integer nsensor                   ! number of sensor for each point
      COMMON /optfield/ data,xcov_scale,xcov_sum,xcov_trace,xcov_lamb,
     &                xcor_sum,xrepl_sum,fit_bart_ind,zcov_noise,
     &                ncov_siz,nsensor
      INTEGER  model  (mpar,mq_post)	! integer discretization of the models
      INTEGER  allmodel  (mpar,mq_post)	! integer discretization of the models
      INTEGER parents(mpar,mq)		! integer discretization of the models
      INTEGER iseed,iallmodel
      COMMON /optmod/model,parents,iseed,allmodel,iallmodel

      REAL rng1,rng2,drng,temp0,niter,temp1,npop,mter
      INTEGER ippd1,ippd2,       ! the index for 2D ppd 
     &        ilin(Mpar),ncurv, nummodes
      REAL pm,px,pu
      COMMON /samppar/rng1,rng2,drng,temp0,niter,temp1,
     &                ilin,pm,px,pu,npop,mter,
     &                ippd1,ippd2,ncurv, nummodes
c***  for use of shape functions
      INTEGER Neof, Neofvar, 
     &        par2phy_eof(Meofvar),par2lay_eof(Meofvar),
     &        par3_eof(Meofvar),
     &        point_eof(Meofvar),nlimits
      REAL    eofcoef(Meofvar,Meof),aeof(Meof),
     &        fmin_lim(Meofvar),fmax_lim(Meofvar)
      common /eof/ Neof,Neofvar,par2phy_eof,par2lay_eof,par3_eof
     &             ,eofcoef,aeof,
     &             point_eof,fmin_lim,fmax_lim,nlimits
c****   for masking of data
      complex    weight(mobs)
      common /weightscom/ weight
c**** for adding noise to the data
      real snr_db
      common /noise_bl/snr_db
c**** for writing of sound speed profile.
      integer iwrite_soundsp
      common /out/ iwrite_soundsp





