      subroutine blug(geo,h,kt,bp,klay,nadd)
c
      implicit none
      include 'Parms_com'
      integer*4 kt,j,klay,nadd,kdo
      real*8 geo(2,5,NLMAX),h(NLMAX),bp(2),
     .   fac,rad,delc,beta,g0,bden,zmax,cmax,cdifmax
c
      if(kt .eq. 2) then
         g0=geo(2,1,klay)
         beta=bp(1)
         fac=geo(1,1,klay)*(1. + beta)
         rad=fac**2 + 2.*g0*fac*h(klay)
         geo(2,1,klay)=dsign(1.d0,fac)*dsqrt(rad) - beta*geo(1,1,klay)
      elseif(kt .eq. 3) then
c: kt=3 means that cp2 is given instead of gradient g for blug:
         beta=bp(1)
         delc=geo(2,1,klay) - geo(1,1,klay)
         g0=delc*(1. + delc/(2.*geo(1,1,klay)*(1. + beta)))/h(klay)
      elseif(kt .eq. 4) then
c: kt=4 means that cp2 is given 2nd, g instead of beta at end:
         g0=bp(1)
         delc=geo(2,1,klay) - geo(1,1,klay)
         bden=g0*h(klay) - delc
         if(bden .eq. 0.) then
            print *,'blug layer ',j,'illegal: cp1,cp2,g gives ',
     .         'a linear profile.'
            stop
         endif
         beta=delc**2/(2.*geo(1,1,klay)*bden) - 1.
      endif
      if(g0 .eq. 0. .or. rad .lt. 0.) then
         print *,'illegal blug profile: j,g0,rad = ',j,g0,rad,
     .      geo(1,1,klay),geo(2,1,klay),beta
         stop
      endif
c
      kdo=klay
      nadd=0
99    call c_dif_max(geo(1,1,kdo),geo(2,1,kdo),h(kdo),beta,
     .   cdifmax,cmax,zmax)
      if(cdifmax .gt. bp(2)) then
c: Create a new layer by splitting current layer:
         call add_layer(kdo,kdo+1,klay,zmax,cmax,geo,h)
         nadd=nadd + 1
      else
         kdo=kdo + 1
      endif
      if(kdo .le. klay) goto 99
c
      return
      end
ccc
      subroutine c_dif_max(c1,c2,h,beta,cdifmax,cmax,zmax)
c
      implicit none
      include 'Parms_com'
      integer*4 j,nvalx
      real*8 c1,c2,h,beta,g0,delc,fac,sg,g0fac2,fac3,cdifmax,zpt,cai,
     .   cblug,rad,cmax,zmax,cdif,facsq,alph
c
      delc=c2 - c1
      g0=delc*(1. + delc/(2.*c1*(1. + beta)))/h
      fac=c1*(1. + beta)
      sg=dsign(1.d0,fac)
      facsq=fac**2
      g0fac2=2.*g0*fac
      fac3=beta*c1
      alph=(1./c2**2 - 1./c1**2)/h
c: Find maximum difference between BLUG profile and 1/c**2 profile:
      nvalx=19
      cdifmax=-1.
      do j=1,nvalx
         zpt=j*h/(nvalx+1)
         cai=c1/dsqrt(1.d0 + alph*c1*c1*zpt)
         rad=facsq + g0fac2*zpt
         cblug=sg*dsqrt(rad) - fac3
         cdif=abs(cai - cblug)
         if(cdif .gt. cdifmax) then
            zmax=zpt
            cmax=cblug
            cdifmax=cdif
         endif
      enddo
c
      return
      end
ccc
      subroutine add_layer(k,kp1,ktot,zpt,cpt,geo,h)
c
      implicit none
      integer*4 k,kp1,ktot,j,j1,j2,ji
      real*8 zpt,cpt,geo(2,5,kp1),h(kp1),fac,hold
c
c: Copy any layers below to one layer down:
      do j=ktot,kp1,-1
         j1=j+1
         h(j1)=h(j)
         do ji=1,2
            do j2=1,5
               geo(ji,j2,j1)=geo(ji,j2,j)
            enddo
         enddo
      enddo
c
      hold=h(k)
c: Make thickness of new layer the leftover amount:
      h(kp1)=h(k) - zpt
c: Make thickness of old layer the depth of max sound speed diff:
      h(k)=zpt
c: Make bottom of new layer same as bottom of old layer:
      do j2=1,5
         geo(2,j2,kp1)=geo(2,j2,k)
      enddo
c
c: Enter blug value of sound speed into profile:
      geo(2,1,k)=cpt
c: Linearly interpolate (except cp) to get bottom of old layer:
ccc   fac=zpt/hold
c: Interpolate other parameters in same non-linear way as BLUG for cp:
      fac=(cpt - geo(1,1,k))/(geo(2,1,kp1) - geo(1,1,k))
      do j2=2,5
         geo(2,j2,k)=geo(1,j2,k) + fac*(geo(2,j2,k)-geo(1,j2,k))
      enddo
c
c: Make top of new layer same as bottom of old layer:
      do j2=1,5
         geo(1,j2,kp1)=geo(2,j2,k)
      enddo
c
c: Increment # of layers:
      ktot=ktot + 1
c
      return
      end
