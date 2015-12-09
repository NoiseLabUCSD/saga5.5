      subroutine smfind(al,ah,af,r,drf,gsl,t,e,trm,trp,sbscat,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
c
c: this subroutine is called from smcaus to find eigenrays near 
c: caustics.
      implicit integer*4(i-n)
      include 'common/srlox'
      include 'common/pathway'
      include 'common/scatarr'
      include 'common/scatcom'
      include 'common/caustix'
      include 'common/freqcom'
      include 'common/svp'
      real al(3),ah(3)
      real*4 angs(0:nasl),frsl(0:nfsl),dbsl(0:nfsl,0:nasl),
     .   phsl(0:nfsl,0:nasl),angb(0:nabl),frbl(0:nfbl),dbbl(0:nfbl,
     .   0:nabl),phbl(0:nfbl,0:nabl)
c
      rtol=rtol/10.
      kcub=1
      if(al(2) .lt. ah(2)) then
         a1=al(1)
         a2=ah(1)
         r1=al(2)
         r2=ah(2)
         dr1=al(3)
         dr2=ah(3)
      else
         a2=al(1)
         a1=ah(1)
         r2=al(2)
         r1=ah(2)
         dr2=al(3)
         dr1=ah(3)
      endif
      kloop=0
10    kloop=kloop + 1
      if(kloop .ge. 6) kcub=0 
      if(kloop .gt. 50) then
         if(abs(diff) .lt. 100.*rtol) goto 77
         print *,'kloop > 50 in smfind.'
         print *,'af,1/af,rf,range = ',af,1./af,rf,r,diff,rtol
         goto 77
      endif
c: bisect to get next guess for af: 
      if(kcub .eq. 1) then
         call polfit(a1,a2,r1,r2,dr1,dr2)
         call cubroot(r,a1,a2,af,kcub)
      else
         af=(a1+a2)/2.
      endif
      call rdrcalc(af,rf,drf,ibd,rbd,drbd)
      diff=rf - r
c: update a1 or a2 if not close enough: 
      if(abs(diff) .gt. rtol) then
         if(diff .gt. 0.) then
            a2=af
            r2=rf
            dr2=drf 
         else
            a1=af
            r1=rf
            dr1=drf 
         endif
         goto 10
      endif
c: calculate eigenray characteristics:  
77    call timatt(af,e,t,ibd) 
c: EKW 5-3-93: rtc_calc handles attenuation correctly:
ccc   call trcalc(af,trm,trp,nps)
      call rtc_multi(af,trm,trp,nps)
      ps=1. - (af*cs)**2
      pr=1. - (af*cr)**2
      gsl=sqrt(af*cs*cr/(sqrt(ps*pr)*abs(rf*drf)))
      call sscalc(ntop,iisl,af,atop,amtop,ctop,ssigma,swind,nfmin,nfmin,
     .   frq,ss,angs,frsl,dbsl,phsl,nfsl,nasl)
      call sscalc(nbas,iibl,af,abas,ambas,cbas,bsigma,bwind,nfmin,nfmin,
     .   frq,bs,angb,frbl,dbbl,phbl,nfbl,nabl)
      sbscat=ss(1,nfmin)*bs(1,nfmin)
c
      rtol=10.*rtol
      return
      end 
