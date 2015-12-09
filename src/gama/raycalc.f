      subroutine raycalc(mr,kr,xl,kl,nray,tf,tfz,frf,ssf,ss2f,bsf,bs2f,
     .   plcoh,plinc,range,xh2,jrec,angs,frsl,dbsl,phsl,angb,frbl,
     .   dbbl,phbl,exparr,ntlay,NLAYM,irec)
c
c: this subroutine processes all nray eigenrays (whose characteristics
c: have been loaded into array xl and kl) for the current s/r
c: configuration.
c
      implicit integer*4(i-n)
      include 'common/srlox'
      include 'common/gamaoptions'
      include 'common/pii'
      include 'common/freqcom'
      include 'common/scatarr'
      include 'common/scatcom'
      include 'common/timing'
      include 'common/pathway'
      include 'common/paths'
      include 'common/pathchr'
      include 'common/svp'
      include 'common/bdcom'
      include 'common/bdchar'
      real*4 xh2(20),plinc(nrtot,nfr2),ths,thr
      real*4 angs(0:nasl),frsl(0:nfsl),dbsl(0:nfsl,0:nasl),
     .   phsl(0:nfsl,0:nasl),angb(0:nabl),frbl(0:nfbl),dbbl(0:nfbl,
     .   0:nabl),phbl(0:nfbl,0:nabl)
      complex*8 plcoh(nrtot,nfr2),tf(nzr*nrangx*nfbb),tfz(nffth1),
     .          zzero,plbmp(107)
      integer*4 kl(2,2500),jj0
      complex cc,plcray,de,expfac,trpfac,exparr(nffth1)
      real range(nrtot),xl(7,2500),frf(nffth1),plray(2)
      real*4 ssf(2,nffth1),ss2f(2,nffth1),bsf(2,nffth1),bs2f(2,nffth1),
     .   phloss,zero
      real*4 tshift
      equivalence(plcray,plray)
      data zzero/(0.,0.)/
      data zero/0./
c
c: initialize transfer function array:  
cpln      jj0=((kr-1)*nzr+jzr-1)*nffth1
      if(iiffi .ne. 0) then
         do 940 jcan=1,nffth1
            tfz(jcan)=zzero
940      continue
      endif
c
cpln      write(6,*)'IIBET: ',iibet
      if(iibet .ne. 0) call eigini(kr,kl,xl,nray,jrec)
      jcaus=0
      ibd=0
      nre=0
c: 9-24-92: leave out if much weaker than cutval:
      cutval=.1*dbsr*cutmag
c: loop over all the eigenrays found for this s/r configuration: 
      do 20 nr=1,nray
         if(iidiag .ne. 0) print *,'start of raycalc: kr,nray = ',
     .      kr,nray
         if(ibd .eq. 0 .and. jcaus .gt. 0) then
            jcaus=0 
            goto 20 
         endif
c     
         kcode=kl(1,nr)
         kbdf=kcode/1000000
         if(kbdf .eq. 0) goto 20
         ibd=mod(kcode/1000,10)
c     : BUG: 3-29-91
         nbas=mod(kcode/10000,100)
c     nbas=mod(kcode/10000,10)
         jcaus=mod(kcode/100,10)
         nps=mod(kcode,100)
         kcode=kl(2,nr)
         icpth=mod(kcode,1000)
         mbp=kcode/1000
         ntop=ncpth(icpth,1)/2
         nbotx=ncpth(icpth,3)/2
         a=xl(1,nr) 
         gsl=abs(xl(2,nr))
         t=xl(3,nr) 
         e=xl(4,nr) 
         trm=xl(5,nr)
         trp=xl(7,nr)
c     : the sign of xlist is that of dr.  subtract off a -90-degree phase
c     : shift for rays that have not grazed the last caustic: 
         if(xl(2,nr) .lt. 0.) then
            nps=nps - 1
            if(jcaus .ne. 2) trp=trp + pieh
         endif
         nrayt=nrayt + 1
c
c: if doppler option set, shift frequencies from frx to frq: 
c:(not yet) if(iidop .eq. 1) call dopshif(a,mr,jcaus)
c
c: compute surface scattering loss for ray:  
         call sscalc(ntop,iisl,a,atop,amtop,ctop,ssigma,swind,1,nfr2,
     .        frq,ss,angs,frsl,dbsl,phsl,nfsl,nasl)
c: compute basement scattering loss for ray:
         call sscalc(nbas,iibl,a,abas,ambas,cbas,bsigma,bwind,1,nfr2,
     .        frq,bs,angb,frbl,dbbl,phbl,nfbl,nabl)
         if(jcaus .eq. 2) then
            call sscalc(ntop,iisl,a2,atop,amtop,ctop,ssigma,swind,1,
     .           nfr2,frq,ss2,angs,frsl,dbsl,phsl,nfsl,nasl)
            call sscalc(nbas,iibl,a2,abas,ambas,cbas,bsigma,bwind,1,
     .           nfr2,frq,bs2,angb,frbl,dbbl,phbl,nfbl,nabl)
         endif
         if(iiffi .ne. 0) then
            call sscalc(ntop,iisl,a,atop,amtop,ctop,ssigma,swind,
     .           1,nffth1,frf,ssf,angs,frsl,dbsl,phsl,nfsl,nasl)
            call sscalc(nbas,iibl,a,abas,ambas,cbas,bsigma,bwind,
     .           1,nffth1,frf,bsf,angb,frbl,dbbl,phbl,nfbl,nabl)
            if(jcaus .eq. 2) then
               call sscalc(ntop,iisl,a2,atop,amtop,ctop,ssigma,swind,
     .              1,nffth1,frf,ss2f,angs,frsl,dbsl,phsl,nfsl,nasl)
               call sscalc(nbas,iibl,a2,abas,ambas,cbas,bsigma,bwind,
     .              1,nffth1,frf,bs2f,angb,frbl,dbbl,phbl,nfbl,nabl)
            endif
         endif
c
cc: SECTION OF CODE FOR EIGENRAYS WITH NO BEAM DISPLACEMENT:
cc
         if(ibd .eq. 0) then
cc
c: jcaus=0 indicates a normal eigenray: 
            if(jcaus .eq. 0) then
               gstr=gsl*trm
               if(iiffi .ne. 0) then
                  trpfac=cmplx(cos(trp),sin(trp))
                  de=exp(df*cmplx(fax2*e,twpie*(t - tmall)))
                  expfac=cmplx(1.,0.)
                  do 231 j=1,nffth1
                     exparr(j)=expfac
                     expfac=expfac*de
 231              continue
                  if((iisl .eq. 0) .and. (iibl .eq. 0)) then
cxx   print *,'nr,trpfac,gstr,trp = ',nr,trpfac,gstr,trp*180./pie
                     do 31 j=1,nffth1
                        plcray=gstr*exparr(j)*trpfac
                        tfz(j)=tfz(j) + plcray
 31                  continue
                  else
cxx   print *,'iisl,iibl,nr,trpfac,gstr,trp = ',iisl,iibl,
cxx  .   nr,trpfac,gstr,trp*180./pie
                     do 331 j=1,nffth1
                        phloss=ssf(2,j) + bsf(2,j)
                        plcray=gstr*exparr(j)*trpfac*ssf(1,j)*
     .                       bsf(1,j)*cmplx(cos(phloss),sin(phloss))
                        tfz(j)=tfz(j) + plcray
 331                 continue
                  endif
               endif
               do 32 kfr=1,nfr2
                  fdep=exp(max(-656.,fax2*frq(kfr)*e))
                  tphi=wom(kfr)*t + trp
                  plcray=gstr*fdep*cmplx(cos(tphi),sin(tphi))
                  if((iisl .ne. 0) .or. (iibl .ne. 0)) then
                     plcray=plcray*ss(1,kfr)*bs(1,kfr)
                     if((ss(2,kfr) .ne. 0.) .or. 
     .                  (bs(2,kfr) .ne. 0.)) then
                        phloss=ss(2,kfr) + bs(2,kfr)
                        plcray=plcray*cmplx(cos(phloss),sin(phloss))
                     endif
                  endif
                  plcoh(mr,kfr)=plcoh(mr,kfr) + plcray
                  plmagsq=plray(1)**2 + plray(2)**2
                  plinc(mr,kfr)=plinc(mr,kfr) + plmagsq
                  if(iibmp .eq. 1) plbmp(kfr)=plcray
 32            continue
               attn=sqrt(plmagsq)
c: if attn of ray at eigenray list freq is too small, don't process:  
               if(attn .lt. cutval) goto 20
c
c: jcaus=1 indicates a shadow zone eigenray: 
            elseif(jcaus .eq. 1) then
               nr2=nr+1
               xaifac=xl(1,nr2)
               dr2c=xl(2,nr2)
               sg=sign(1.,dr2c)
               s1fac=xl(3,nr2)
               s2fac=xl(4,nr2)
               rc=xl(5,nr2)
               if(iiffi .ne. 0) then
                  trp2=trp + piequ
                  trpfac=cmplx(cos(trp2),sin(trp2))
                  de=exp(df*cmplx(fax2*e,twpie*(t - tmall)))
                  expfac=cmplx(1.,0.)
                  do 233 j=1,nffth1
                     exparr(j)=expfac
                     expfac=expfac*de
 233              continue
                  if((iisl .eq. 0) .and. (iibl .eq. 0)) then
cvd	cncall
                     do 33 j=2,nffth1
                        w16f=(twpie*frf(j))**(.16667)
                        w23f=w16f**4
                        xairy=w23f*xaifac
                        call airy(-1.*xairy,ai,aip)
                        s1=s1fac*w16f
                        s2=s2fac/w16f
                        plcray=cmplx(s1*ai,s2*aip)*exparr(j)*
     .                  trpfac
                        tfz(j)=tfz(j) + plcray
 33                  continue
                  else
cvd	cncall
                     do 333 j=2,nffth1
                        w16f=(twpie*frf(j))**(.16667)
                        w23f=w16f**4
                        xairy=w23f*xaifac
                        call airy(-1.*xairy,ai,aip)
                        s1=s1fac*w16f
                        s2=s2fac/w16f
                        phloss=ssf(2,j) + bsf(2,j)
                        plcray=cmplx(s1*ai,s2*aip)*exparr(j)*
     .                       trpfac*ssf(1,j)*bsf(1,j)*
     .                       cmplx(cos(phloss),sin(phloss))
                        tfz(j)=tfz(j) + plcray
 333                 continue
                  endif
               endif
               do 34 kfr=1,nfr2
                  xairy=w23(kfr)*xaifac
                  call airy(-1.*xairy,ai,aip)
                  fdep=exp(max(-656.,fax2*frq(kfr)*e))
                  s1=s1fac*fdep*w16(kfr)
                  s2=s2fac*fdep/w16(kfr)
                  cc=cmplx(s1*ai,s2*aip)
                  tphi=wom(kfr)*t + trp + piequ
                  if((iisl .ne. 0) .or. (iibl .ne. 0)) then
                     cc=cc*ss(1,kfr)*bs(1,kfr)
                     if((ss(2,kfr) .ne. 0.) .or. 
     .                  (bs(2,kfr) .ne. 0.)) then
                        tphi=tphi + ss(2,kfr) + bs(2,kfr)
                     endif
                  endif
                  plcray=cc*cmplx(cos(tphi),sin(tphi))
                  plcoh(mr,kfr)=plcoh(mr,kfr) + plcray
                  plmagsq=plray(1)**2 + plray(2)**2
                  plinc(mr,kfr)=plinc(mr,kfr) + plmagsq
                  if(iibmp .eq. 1) plbmp(kfr)=plcray
 34            continue
c: 4-11-90: CHANGED 53-DB EXTENT TO 40-DB EXTENT (Ai(3.33)/Ai(0) = .01):
c           delta=-4.12*(rng - rc)/xairy
               delta=-3.33*(rng - rc)/xairy
c: BUG 10-30-90: nfeig not set here:
cxx         dellam=delta/bdlam(nfeig)
               dellam=delta/frq(nfr2)
               attn=sqrt(plmagsq)
c     print *,'raycalc shadow: ',attn,phr*piedeg
c     print *,xairy,ai,aip
c     print *,s1,s2,tphi*piedeg
c: if attn of ray at eigenray list freq is too small, don't process:  
               if(attn .lt. cutval) goto 20
               gsl=attn/(trm*fdep*ss(1,nfr2)*bs(1,nfr2))
c
c: jcaus=2 indicates a pair of caustic eigenrays (insonified zone): 
            elseif(jcaus .eq. 2) then
               nr2=nr+1
               nps2=mod(kl(1,nr2),100)
               a2=xl(1,nr2)
               gsl2=abs(xl(2,nr2))
               t2=xl(3,nr2)
               e2=xl(4,nr2)
               trm2=xl(5,nr2)
               trp2=xl(7,nr2)
c: check for grazing of caustic: 
               if(xl(2,nr2) .lt. 0.) nps2=nps2 - 1
               tdif=.75*(t2 - t) 
               phdif=.75*(trp2 - trp)
               tavg=(t + t2)/2.
               phavg=(trp + trp2)/2. + piequ
               sg=sign(1.,tdif)
               gstr=gsl*trm
               gstr2=gsl2*trm2
               if(iiffi .ne. 0) then
                  trpfac=cmplx(cos(phavg),sin(phavg))
                  de=exp(df*cmplx(0.,twpie*(tavg - tmall)))
                  expfac=cmplx(1.,0.)
                  do 235 j=2,nffth1
                     exparr(j)=expfac
                     expfac=expfac*de
 235              continue
                  if((iisl .eq. 0) .and. (iibl .eq. 0)) then
cvd	cncall
                     do 35 j=2,nffth1
                        xairy=abs(twpie*frf(j)*tdif + phdif)**(.66667)
                        xaip=xairy**(.25)
                        call airy(-1.*xairy,ai,aip)
                        s1=sqpie*(gstr2 + gstr)*xaip
                        s2=sg*sqpie*(gstr2 - gstr)/xaip
                        fdep=exp(max(-656.,fax2*frf(j)*(e+e2)/2.))
                        plcray=fdep*cmplx(s1*ai,s2*aip)*exparr(j)*
     .                  trpfac
                        tfz(j)=tfz(j) + plcray
 35                  continue
                  else
cvd	cncall
                     do 335 j=2,nffth1
                        xairy=abs(twpie*frf(j)*tdif + phdif)**(.66667)
                        xaip=xairy**(.25)
                        call airy(-1.*xairy,ai,aip)
                        fac2=gstr2*ss2f(1,j)*bs2f(1,j)
                        fac1=gstr*ssf(1,j)*bsf(1,j)
                        s1=sqpie*(fac2 + fac1)*xaip
                        s2=sg*sqpie*(fac2 - fac1)/xaip
                        fdep=exp(max(-656.,fax2*frf(j)*(e+e2)/2.))
                        phloss=(ssf(2,j)+bsf(2,j)+ss2f(2,j)+
     .                  bs2f(2,j))/2.
                        plcray=fdep*cmplx(s1*ai,s2*aip)*exparr(j)*
     .                       trpfac*cmplx(cos(phloss),sin(phloss))
                        tfz(j)=tfz(j) + plcray
 335                 continue
                  endif
               endif
               do 36 kfr=1,nfr2
                  xairy=abs(wom(kfr)*tdif + phdif)**(.66667)
                  xaip=xairy**(.25)
                  call airy(-1.*xairy,ai,aip)
                  fac2=gstr2*ss2(1,kfr)*bs2(1,kfr)
                  fac1=gstr*ss(1,kfr)*bs(1,kfr)
                  s1=sqpie*(fac2 + fac1)*xaip
                  s2=sg*sqpie*(fac2 - fac1)/xaip
                  cc=cmplx(s1*ai,s2*aip)
                  tphi=wom(kfr)*tavg + phavg
                  fdep=exp(max(-656.,fax2*frq(kfr)*(e+e2)/2.))
                  plcray=fdep*cc*cmplx(cos(tphi),sin(tphi))
                  if((iisl .ne. 0) .or. (iibl .ne. 0)) then
                     if((ss(2,kfr) .ne. 0.) .or. 
     .                  (bs(2,kfr) .ne. 0.)) then
                        phloss=(ss(2,kfr)+bs(2,kfr)+ss2(2,kfr)+
     .                       bs2(2,kfr))/2.
                        plcray=plcray*cmplx(cos(phloss),sin(phloss))
                     endif
                  endif
                  plcoh(mr,kfr)=plcoh(mr,kfr) + plcray
                  plmagsq=plray(1)**2 + plray(2)**2
                  plinc(mr,kfr)=plinc(mr,kfr) + plmagsq
                  if(iibmp .eq. 1) plbmp(kfr)=plcray
 36            continue
               attn=sqrt(plmagsq)
c     print *,'raycalc dublet: ',attn,phr*piedeg
c     print *,xairy,ai,aip
c     print *,s1,s2,tphi*piedeg
c: if strength of ray at eig list freq is too small, don't process: 
               if(attn .lt. cutval) goto 20
c
               gsl=attn/(trm*fdep*ss(1,nfr2)*bs(1,nfr2))
               gsl2=attn/(trm2*fdep*ss2(1,nfr2)*bs2(1,nfr2))
               nrayt=nrayt + 1
            else
               print *,'bug in raycalc: jcaus = ',jcaus
            endif
cc
cc: SECTION OF CODE FOR EIGENRAYS WITH BEAM DISPLACEMENT:
cc
         else
cc
            call bdsort(kl,xl,nray,nr,jcaus,mbp,kbdf1,kbdf2,
     .           nreig,nreig2)
            call bdfit(kbdf1,kbdf2)
            kk=kbdf2+1
            if(iiffi .ne. 0) then
               kbdf0=kbdf1
               jmin=max0(1,nint(bdf(kbdf1-1)/df) + 1)
               jmax=min0(nffth1,nint(bdf(kbdf2+1)/df) + 1)
               if((iisl .eq. 0) .and. (iibl .eq. 0)) then
cvd	cncall
                  do 40 j=jmin,jmax
                     fj=(j-1)*df
                     if(fj .gt. bdf(kbdf0)) kbdf0=min0(kk,kbdf0 + 1)
                     t=parma(kbdf0,1)*fj + parmb(kbdf0,1)
                     trp=parma(kbdf0,2)*fj + parmb(kbdf0,2)
                     ph=twpie*fj*(t-tmall) + trp
                     gstr=parma(kbdf0,3)*fj + parmb(kbdf0,3)
                     plcray=gstr*cmplx(cos(ph),sin(ph))
                     if(nparm(kbdf0) .eq. 6) then
                        xairy=parma(kbdf0,4)*fj + parmb(kbdf0,4)
                        s1=parma(kbdf0,5)*fj + parmb(kbdf0,5)
                        s2=parma(kbdf0,6)*fj + parmb(kbdf0,6)
                        call airy(-xairy,ai,aip)
                        plcray=plcray*cmplx(s1*ai,s2*aip)
                     endif
c     write(80+nrayt,660) fj,kbdf0,t,trp*piedeg,gstr,xairy,s1,s2
c660   format(f9.3,1x,i2,6(1x,e11.6))
                     tfz(j)=tfz(j) + plcray
 40               continue
               else
cvd	cncall
                  do 340 j=jmin,jmax
                     fj=(j-1)*df
                     if(fj .gt. bdf(kbdf0)) kbdf0=min0(kk,kbdf0 + 1)
                     t=parma(kbdf0,1)*fj + parmb(kbdf0,1)
                     trp=parma(kbdf0,2)*fj + parmb(kbdf0,2)
                     ph=twpie*fj*(t-tmall) + trp
                     gstr=parma(kbdf0,3)*fj + parmb(kbdf0,3)
                     plcray=gstr*cmplx(cos(ph),sin(ph))
                     if(nparm(kbdf0) .eq. 6) then
                        xairy=parma(kbdf0,4)*fj + parmb(kbdf0,4)
                        s1=parma(kbdf0,5)*fj + parmb(kbdf0,5)
                        s2=parma(kbdf0,6)*fj + parmb(kbdf0,6)
                        call airy(-xairy,ai,aip)
                        plcray=plcray*cmplx(s1*ai,s2*aip)
                     endif
c     write(80+nrayt,660) fj,kbdf0,t,trp*piedeg,gstr,xairy,s1,s2
c660   format(f9.3,1x,i2,6(1x,e11.6))
                     if(jcaus .eq. 2) then
                        avscat=(ssf(1,j)*bsf(1,j) + ss2f(1,j)*
     .                  bs2f(1,j))/2.
                        phloss=(ssf(2,j)+bsf(2,j)+ss2f(2,j)+
     .                  bs2f(2,j))/2.
                        plcray=plcray*avscat*
     .                  cmplx(cos(phloss),sin(phloss))
                     else
                        phloss=ssf(2,j) + bsf(2,j)
                        plcray=plcray*ssf(1,j)*bsf(1,j)*
     .                       cmplx(cos(phloss),sin(phloss))
                     endif
                     tfz(j)=tfz(j) + plcray
 340              continue
               endif
            endif
c
            kbdf0=kbdf1
            flo=bdf(kbdf1-1)
            fup=bdf(kbdf2+1)
            do 50 kfr=1,nfr2
               fj=frq(kfr)
               if((fj .gt. fup) .or. (fj .lt. flo)) goto 50
 52            if(fj .lt. bdf(kbdf0-1)) then
                  kbdf0=kbdf0-1
                  goto 52
               endif
 54            if(fj .gt. bdf(kbdf0)) then
                  kbdf0=kbdf0+1
                  goto 54
               endif
               t=parma(kbdf0,1)*fj + parmb(kbdf0,1)
               trp=parma(kbdf0,2)*fj + parmb(kbdf0,2)
               ph=twpie*fj*t + trp
               gstr=parma(kbdf0,3)*fj + parmb(kbdf0,3)
               plcray=gstr*cmplx(cos(ph),sin(ph))
c     print *,'bd pl calc: fj,t,trp,gstr = ',fj,t,trp*piedeg,gstr
               if(nparm(kbdf0) .eq. 6) then
                  xairy=parma(kbdf0,4)*fj + parmb(kbdf0,4)
                  call airy(-xairy,ai,aip)
                  s1=parma(kbdf0,5)*fj + parmb(kbdf0,5)
                  s2=parma(kbdf0,6)*fj + parmb(kbdf0,6)
c     print *,'bd caustic parms: xairy,s1,s2 = ',xairy,s1,s2
                  plcray=plcray*cmplx(s1*ai,s2*aip)
               endif
               if((iisl .ne. 0) .or. (iibl .ne. 0)) then
                  if(jcaus .eq. 2) then
                     avscat=(ss(1,kfr)*bs(1,kfr)+ss2(1,kfr)*
     .               bs2(1,kfr))/2.
                     plcray=plcray*avscat
                     if((ss(2,kfr) .ne. 0.) .or. 
     .                  (bs(2,kfr) .ne. 0.)) then
                        phloss=(ss(2,kfr) + bs(2,kfr) + ss2(2,kfr)
     .                       + bs2(2,kfr))/2.
                        plcray=plcray*cmplx(cos(phloss),sin(phloss))
                     endif
                  else
                     plcray=plcray*ss(1,kfr)*bs(1,kfr)
                     if((ss(2,kfr) .ne. 0.) .or. 
     .                  (bs(2,kfr) .ne. 0.)) then
                        phloss=ss(2,kfr) + bs(2,kfr)
                        plcray=plcray*cmplx(cos(phloss),sin(phloss))
                     endif
                  endif
               endif
               plcoh(mr,kfr)=plcoh(mr,kfr) + plcray
               plmagsq=plray(1)**2 + plray(2)**2
               plinc(mr,kfr)=plinc(mr,kfr) + plmagsq
               if(iibmp .eq. 1) plbmp(kfr)=plcray
 50         continue
            attn=sqrt(plmagsq)
c: if attn of ray at eigenray list freq is too small, don't process:  
            if(attn .lt. cutval) goto 20
c: if eigenray list, bar graph, or picture, set the eigenray char
c: to those at f=freq.  bdsort set nreig to the ray at that freq:
            if((iibet .ne. 0) .or. (iipic .ne. 0)) then
               if(nreig .eq. 0) goto 67
               a=xl(1,nreig)
               gsl=abs(xl(2,nreig))
               t=xl(3,nreig)
               e=xl(4,nreig)
               trm=xl(5,nreig)
               trp=xl(7,nreig)
               kcode=kl(1,nreig)
               jcaus0=jcaus
               jcaus=mod(kcode/100,10)
               if(jcaus .eq. 0 .and. jcaus0 .eq. 0) jcaus=jcaus0
               nps=mod(kcode,100)
               if(xl(2,nreig) .lt. 0.) then
                  nps=nps - 1
               if(jcaus .ne. 2) trp=trp + pieh
            endif
            if(jcaus .eq. 1) then
               nr2=nreig+1
               xaifac=xl(1,nr2)
               dr2c=xl(2,nr2)
               rc=xl(5,nr2)
               delta=-3.33*(rng - rc)/(bdom23(nfeig)*xaifac)
               dellam=delta/bdlam(nfeig)
               gsl=attn/(trm*expo(frq(nfr2),e)*ss(1,nfr2)*bs(1,nfr2))
            elseif(jcaus .eq. 2) then
               nr2=nreig+1
               nps2=mod(kl(1,nr2),100)
               a2=xl(1,nr2)
               gsl2=abs(xl(2,nr2))
               t2=xl(3,nr2)
               e2=xl(4,nr2)
               trm2=xl(5,nr2)
               trp2=xl(7,nr2)
               if(xl(2,nr2) .lt. 0.) nps2=nps2 - 1
               gsl=attn/(trm*expo(frq(nfr2),e)*ss(1,nfr2)*bs(1,nfr2))
               gsl2=attn/(trm2*expo(frq(nfr2),e2)*ss2(1,nfr2)*
     .              bs2(1,nfr2))
               nrayt=nrayt + 1
            endif
         endif
cc
      endif
cc
cc: END OF SEPARATE TREATMENT OF BEAM DISPLACED RAYS.
c
c
 199  continue
      if(iibet .ne. 0) then
         if(plcray .ne. zzero) then
            phr=atan2(plray(2),plray(1))*piedeg
         else
            phr=0.
            print *,'Warning: Ray with zero magnitude found: '
            print *,'kr,nr,nray = ',kr,nr,nray
         endif
         trp=mod(trp,twpie)*piedeg
         if(abs(trp) .gt. 180.) trp=trp - sign(360.,trp)
         nre=nre+1
         dbattn=9999.99
         if(attn .ne. 0.) dbattn=20.*log10(attn)
      endif
c     : plot eigenray picture if desired: 
      if(iipic .ge. 1) then
         if((jcaus .eq. 1) .and. (rc .gt. 1.25*rng)) then
c     print *,'out of bounds caustic ray not plotted',
c     .              ' for range = ',range(mr)
         else
            bdomega=bdom(nfeig)
            call rayplot(a,range(mr),jcaus,ntlay,NLAYM)
         endif
      endif
c     
      if(iieig .eq. 1) then
         if(iidiag .gt. 1)  print *,'calling eigwrit ...'
         call eigwrit(nre,jcaus,a,t,e,trm,trp,gsl,
     .        dbattn,nps,rc,dr2c,delta,dellam,ss(1,nfr2),bs(1,nfr2),
     .        phr,ibd,ntlay,NLAYM)
         if(iidiag .gt. 1)  print *,'done eigwrit ...'
      endif
      if(iibar .eq. 1) call bargr(nre,jcaus,t,e,dbattn,a,phr)
      if(jcaus .ge. 2) then
         nre=nre+1
         icpth=mod(kl(2,nr2),1000)
         ntop=ncpth(icpth,1)/2
         nbotx=ncpth(icpth,3)/2
         mbp=kl(2,nr2)/1000
         trp2=mod(trp2,twpie)*piedeg
         if(trp2 .gt. 180.) trp2=trp2 - 360.
         if(iieig .eq. 1) call eigwrit(nre,jcaus,a2,t2,e2,trm2,trp2,
     .        gsl2,dbattn,nps2,rc,dr2c,delta,dellam,
     .        ss2(1,nfr2),bs2(1,nfr2),phr,ibd,ntlay,NLAYM)
         if(iipic .ge. 1) call rayplot(a2,range(mr),jcaus,
     .        ntlay,NLAYM)
         if(iibar .eq. 1) call bargr(nre,jcaus,t2,e2,dbattn,a2,phr)
      endif
c
 67   continue
c     : write out eigenray angles and complex fields to beam pattern file:
      if(iibmp .eq. 1) then
         ths=acos(a*cs)*piedeg
         if(fdir(icpth) .eq. 'd') ths=-1.*ths
         thr=acos(a*cr)*piedeg
         if(ldir(icpth) .eq. 'u') thr=-1.*thr
         write(18) ths,thr,(plbmp(jcan),jcan=1,nfr)
      endif
 20   continue
c
c: write out impulse response and/or fft files: 
      if(iiffi .ne. 0) then
c: take complex conj of transfer function to agree with alliant fft:
cpln         do 944 jcan=1,nffth1
cpln            tfz(jcan)=conjg(tfz(jcan))
cpln 944     continue
c     : make sure imaginary parts of first and last FFT bin are zero:
cpln         tfz(1)=cmplx(real(tfz(1)),0.)
cpln         tfz(1)=cmplx(real(tfz(1)),zero)
cpln         tfz(nffth1)=cmplx(real(tfz(nffth1)),0.)
cpln         tfz(nffth1)=cmplx(real(tfz(nffth1)),zero)
         do 110 j=1,nffth1
            tfmag=abs(tfz(j))
            if(tfmag .ne. 0.) then
               tfmax=max(tfmax,tfmag)
               tfmin=min(tfmin,tfmag)
            endif
 110     continue
         if(iiir .ne. 0) then
c: taper upper end of transfer function with a raised cosine: 
            do 100 j=ntap,nffth1
               tfz(j)=tfz(j)*(1. + .999*cos(dtap*(j-ntap)))/2.
 100        continue
         endif
         if(iifft .ne. 0) then
c: write out nfft/2 + 1 complex numbers to the fft file:
cpln            irec=irec + 1
cpln            write(16,rec=irec) (xh2(jcan),jcan=1,20),
cpln     .         (tfz(jcan),jcan=1,nffth1)
cpln            write(16,*) (xh2(jcan),jcan=1,20)
cpln
cpln            write(6,*)'zindx, rindx: ',jzr,kr
cpln            write(6,*)'nzr,nrangx,nfbb: ',nzr,nrangx,nfbb
cpln            pause
cpln            do jcan=1,nffth1
cpln            do jcan=nf1,nf2
            do jcan=1,nfbb
c               
               jj0=((kr-1)*nzr+jzr-1)*nfbb+jcan
cpln               write(6,*)'kr,jzr,jcan,jj0: ',kr,jzr,jcan,jj0
               tshift=(jcan+nf1-2)*twpie*df*xh2(1)
               tf(jj0)=tfz(jcan+nf1-1)*
     .         cmplx(cos(tshift),sin(tshift))
               tfz(jcan+nf1-1)=
     .         tfz(jcan+nf1-1)*cmplx(cos(tshift),sin(tshift))
cpln             write(19,*)real(tfz(jcan+nf1-1)),aimag(tfz(jcan+nf1-1))
c
            end do
         endif
         if(iiir .eq. 1) then
c: take inverse fft to get impulse response:
cpln            tfz(1)=cmplx(real(tfz(1)),0.)
            tfz(1)=cmplx(real(tfz(1)),zero)
            call csfft1d(nfft,tfz,1)
c: write out nfft real # (nffth complex #) to impulse response file:
            write(50) (tfz(jcan),jcan=1,nffth)
         endif
      endif
c
      if(iibet .ne. 0) call eigend(plinc(mr,nfr2),plcoh(mr,nfr2))
c
      return
      end 
