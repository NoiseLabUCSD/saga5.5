      subroutine rzlay(jprof,a,ainv,zn,cn,cnp,g0,g1,g2,rmax,
     .   dr,dz,npt,nptmax,kturn,iilast,radd,zadd)
c
      implicit integer*4(i-n)
      include 'common/kwine'
      include 'common/depth'
      include 'common/bottom' 
      real dr(nptmax),dz(nptmax)
c
c: EKW 9-21-92: Included iilast,radd,zadd in order to avoid huge
c: ray files.  Does not require one point to be plotted for each layer
c: in ocean.
c
c: compute constants which do not change while tracing the ray in
c: the current layer.
      if(jprof .eq. 3) then
         psi=sqrt(1.-(a*cn)**2)
         phi1=psi/cn
         f1=g1 + g2**2*phi1**2
         f3=g1 - g0*g2/2.
         rabf1=sqrt(abs(f1))
      else
         psi=sqrt(1.-(a*cn)**2)
      endif
c: find the vertical distance, zmax, the ray travels in the layer.
      if((ainv .ge. cnp) .or. ((jprof .eq. 1) .and. (g0 .eq. 0.)))
     .      then
         kturn=0
         zmax=zn
      else
         kturn=1
         if(jprof .eq. 1) then
            zmax=(ainv - cn)/g0
         elseif(jprof .eq. 2) then
            fac=cn*(1.+g1)
            zmax=((ainv + g1*cn)**2 - fac**2)/(2.*g0*fac)
         elseif(jprof .eq. 3) then
            zmax=zwine(sg,zn,g0,g1,g2,phi1,f1,g1-g0*g2)
         endif
      endif
c
c: find the horizontal distance rdif covered in the layer and use the 
c: appropriate number of plotting points in the layer, npt.
      rdif=rnglay(jprof,a,cn,g0,g1,g2,zmax)
      nptmin=0
      if(kturn .eq. 1 .or. iilast .eq. 1) nptmin=1
      npt=min0(int(.8*nptmax),max0(int(.9*nptmax*(radd+rdif)/rmax),
     .   nptmin))
      if(npt .eq. 0) then
         radd=radd + rdif
         zadd=zadd + zmax
      else
         fac=zmax/float(npt)
         do 20 k=1,npt
            dzlay=float(k)*fac
            dz(k)=zadd + dzlay
            dr(k)=radd + rnglay(jprof,a,cn,g0,g1,g2,dzlay)
20       continue
         radd=0.
         zadd=0.
      endif
c
      return
      end 
