      subroutine mode_branch(iiref,k_st,iiccwx,iimstx,iilk0,
     .   iidup,ndup,ndup_max)
c
c: Finds modes along a given branch (|R1R2|=1 contour) starting at k_st
c: in the direction given by iiccwx (1=to left in k plane, -1=to right).
c: For iidup=1, check found modes for duplicated with modes kn(1:nm_ok).
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 nm_ok,iiref,iiccwx,iimstx,iilk0,iidup,iduct,ndup,
     .   ndup_max,nm1
      complex*16 k_st,k0,rr0(3,4),r1r2(3,4)
c
      kn(nmode)=k_st
      iiccw=iiccwx
      iimst=iimstx
      iilk=iilk0
      nhigh_max=3
      if(iilk .gt. 0) nhigh_max=1
c
      iifail=0
      iidone=0
      nhigh=0
      nrise=0
      ndup=0
      nm_ok=nm_put
c
c: Start at k of reference depth and find |R1R2|=1 contour:
cc    kcut=dcmplx((1.d0-1.d-8)*dreal(xkref),dimag(xkref))
      if(iiref .eq. 1) then
c: Starting at reference depth, move slightly to left:
cc       k0=dcmplx((1.d0-1.d-8)*dreal(k_st),dimag(k_st))
         if(nsvmin .eq. nlay) then
            k0=cdsqrt(k_st*k_st - .5d0*etasq(nlay))
            call sheet_init(k0,1,iish,iish_ref)
            if(dreal(k0) .gt. dreal(k_st)) iish_ref(1)=-1
cc          print *,'HSP k0 = ',dreal(k0)/kw0,dimag(k0)*8685.9
         else
c: Move to left by 1/1000 of nominal mode spacing:
            k0=k_st - min(.001d0*pie/Htot,.1d0*dreal(k_st))
            call sheet_init(k0,1,iish,iish_ref)
         endif
      else
c: Starting at point on contour, start on k_st:
         k0=k_st
      endif
      call r1r2_calc(k0,rr0,2,0)
cc    call fix_path(k0,rr0,0)
cc    if(iidone .eq. 1) goto 88
      if(iimt .eq. 1) call mode_traj(k0,rr0,0)
c
      nm1=nmode+1
10    continue
         call eig_findm(k0,rr0,r1r2,nm1)
         if(iifail .eq. 1) then
cMRF            print *,'Mode finding failure in duct ',kduct,
cMRF     .         ' at depth ',zduct(kduct)
cMRF            print *,'Continuing with modes found ...'
            iifail=0
            return
         elseif(iifail .eq. 2) then
c: Mode finding failed due to island mode at XKREF; don't warn since this
c: is normal under most circumstances:
            iifail=0
            return
         endif
c: Return if branch line mode found on main contour and BLMs will
c: be found later (for all-fluid when ref depth can be placed in 
c: halfspace):
         if(iiblm .gt. 0) return
         if(iidone .eq. 1) then
c: Subtract off any modes that were too weak at end of mode search:
            if(nhigh .gt. 1) nm_put=nm_put + 1 - nhigh
            return
         endif
         if(iidup .eq. 1) then
            call duct_dupl(nm_ok,ndup,r1r2,phi,dphi,iduct)
            if(ndup .ge. ndup_max) return
         else
            nm_put=nm_put + 1
         endif
      goto 10
c
      end
