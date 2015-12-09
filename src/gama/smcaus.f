      subroutine smcaus(a1,al,ah,a2,rsm,s1,s2,xairy,nogo,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
c
c: this subroutine calculates s1 and s2 at the range rsm slightly
c: on the insonified zone side of the caustic.
c
      implicit integer*4(i-n)
      include 'common/pii'
      include 'common/freqcom'
      include 'common/caustix'
      include 'common/bdcom'
      real a1(3),al(3),ah(3),a2(3)
c
c     print *,'a1,al,ah,a2 = ',a1,al,ah,a2
c added pln 241100
      nogo=0
      if((rsm-a1(2))*(rsm-al(2)) .ge. 0.) then
         nogo=1
         return
      elseif((rsm-a2(2))*(rsm-ah(2)) .ge. 0.) then
         nogo=1
         return
      endif
      call smfind(a1,al,asm1,rsm,drsm1,gsl1,t1,e1,trm1,trp1,sbscat1,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
      call smfind(ah,a2,asm2,rsm,drsm2,gsl2,t2,e2,trm2,trp2,sbscat2,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
      tdif=.75*(t2 - t1) 
      phdif=.75*(trp2 - trp1)
      sg=sign(1.,tdif)
      gstr1=gsl1*trm1
      gstr2=gsl2*trm2
      xairy=abs(bdom(kbdf)*tdif + phdif)**(.66667)
      xaip=xairy**(.25)
      fdep1=expo(bdf(kbdf),e1)
      fdep2=expo(bdf(kbdf),e2)
c: NOTE: phases of bottom loss terms cannot be included here:
      fac1=fdep1*gstr1*sbscat1
      fac2=fdep2*gstr2*sbscat2
      s1=sqpie*(fac2 + fac1)*xaip
      s2=sg*sqpie*(fac2 - fac1)/xaip
c     print *,'smcaus: xairy,s1,s2 = ',xairy,s1,s2
c
      return
      end
