      subroutine dubloma(amono,anx,amonox,a1,a2,kf,nlom,xlist,
     .   klist,nlist,range,wk,ap,kkk,xxx,yyy,nda,tmin,dbmax,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
c
c: this subroutine is called when an double extremum is detected
c: in found, as called from rcheck.  the two extrema are found by
c: calling loma twice, and the monotonic intervals are check by calling
c: rcheck. amono, amonox, a1, and a2 are updated so that the next 
c: shadow call is correct.
c
      implicit integer*4(i-n)
      include 'common/dublcom'
      include 'common/saga_gama_fail'
      real amono(3),anx(3),amonox(3),a1(3),a2(3),a3(3),a4(3),a5(3),a6(3)
c
c     print *,'dubloma called: '
c: find first caustic:
cpln      write(6,*)'LOMA #2'
      call loma(amono,abad,a3,a4,kf1,nlom)
c: check range interval up to first caustic:
cpln      write(6,*)'dubloma: rcheck #1'
      call rcheck(amono,a3,nlom,xlist,klist,nlist,range,wk,ap,kkk,xxx,
     .   yyy,nda,tmin,dbmax,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
      if(iifail.eq.1) return
c: if a previous caustic had been found, do caustic corrections:
      if(nlom .gt. 1) then
         call shadow(a1,a2,kf,amonox,a3,iifc,
     .   xlist,klist,nlist,range,wk,ap,nda,tmin,dbmax,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl,0)
      end if
c: find second caustic:
cpln      write(6,*)'LOMA #3'
      call loma(abad,anx,a5,a6,kf2,nlom)
c: check range interval from first caustic to second caustic:
cpln      write(6,*)'dubloma: rcheck #2'
      call rcheck(a4,a5,nlom,xlist,klist,nlist,range,wk,ap,kkk,xxx,
     .   yyy,nda,tmin,dbmax,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
      if(iifail.eq.1) return
c: do caustic corrections for first caustic:
      call shadow(a3,a4,kf1,a2,a5,iifc,
     .   xlist,klist,nlist,range,wk,ap,nda,tmin,dbmax,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl,0)
c: do last rcheck back in rays:
c     call rcheck(a6,anx,nlom,xlist,klist,nlist,range,wk,ap,kkk,xxx,
c    .   yyy,nda,tmin,dbmax,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
      do 940 jcan=1,3
         amono(jcan)=a6(jcan)
         amonox(jcan)=a4(jcan)
         a1(jcan)=a5(jcan)
         a2(jcan)=a6(jcan)
940   continue
cmay  amono(1:3)=a6(1:3)
cmay  amonox(1:3)=a4(1:3)
cmay  a1(1:3)=a5(1:3)
cmay  a2(1:3)=a6(1:3)
c
      return
      end 
