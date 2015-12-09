      subroutine winefit(k,c1,c2,b1,b2,h,dc2max,jok)
c
c: this subroutine fits the sound speeds and gradients to a Weinber
c: profile.
c
      implicit integer*4(i-n)
      include 'common/svp'
c
      jok=0
      sumd=(c1-c2)*(c1+c2)/h
      cc=c1*c2
      ccsq=cc**2
      bet=sumd/ccsq
      g00=-2.*b1/c1**3
      alp=-2.*b2/c2**3
      rad=sumd**2 - 4.*b1*b2*cc
      if(rad .lt. 0.) then  
c: check if rad<0 due to round-off error:
         check=abs(rad)/sumd**2
         if(check .lt. 1.e-6) then
            rad=0.
         else
            print *,'rad < 0 in winefit: rad,check = ',rad,check
            return
         endif
      endif
      sqrad=sqrt(rad)/ccsq
      facp=(bet+sqrad)/alp
      facm=(bet-sqrad)/alp
c: fit the curve in the downward direction.
      g0(k,2)=g00
      g1p=(bet*facp**2 - g00)/h
      g2p=(facp-1.)/h
      g1m=(bet*facm**2 - g00)/h
      g2m=(facm-1.)/h
      call cdcwine(c1,g0(k,2),g1p,g2p,0.,cp,dcp,dcp2,2)
      call cdcwine(c1,g0(k,2),g1p,g2p,h,cp,dcp,dcph2,2)
      call cdcwine(c1,g0(k,2),g1m,g2m,0.,cm,dcm,dcm2,2)
      call cdcwine(c1,g0(k,2),g1m,g2m,h,cm,dcm,dcmh2,2)
      dc2maxp=max(abs(dcp2),abs(dcph2))
      dc2maxm=max(abs(dcm2),abs(dcmh2))
      if(dc2maxp .lt. dc2maxm) then
         g1(k,2)=g1p
         g2(k,2)=g2p
         sg=1.
         dc2max=dc2maxp
      else
         g1(k,2)=g1m
         g2(k,2)=g2m
         sg=-1.
         dc2max=dc2maxm
      endif
c: fit the curve in the upward direction now:
      bet=-1.*bet
      g0(k,1)=-1.*alp
      alp=-1.*g00
      fac=(bet+sg*sqrad)/alp
      g1(k,1)=(bet*fac**2 - g0(k,1))/h
      g2(k,1)=(fac - 1.)/h
c: do not allow g2=0 or else expressions for travel time will blow up:
      if(abs(g2(k,2)) .lt. 1.e-10 .or. abs(g2(k,1)) .lt. 1.e-10) return
      jok=1
c
      return
      end 
