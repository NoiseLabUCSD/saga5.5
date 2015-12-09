      subroutine fix_path(k0,rr0,iileft)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 iileft,j
      complex*16 k0,rr0(3,4),klast,lnlast,dlnlast,branch,delk,kcut0
      real*8 ph_des0,dmag_left
c
      k0=kcut
      kcut0=kcut
      call sheet_init(k0,1,iish,iish_ref)
c
c: phcut is phase along |R1R2|=1 contour at point where contour crossed
c: branch cut:
      phcut=dimag(lncut)
c
      if(iimt .ne. 0) call mode_traj(k0,rr0,0)
      call r1r2_calc(k0,rr0,2,0,jjfail)
      if(jjfail .gt. 0) return
      lnlast=rr0(1,4)
      if(dabs(dreal(lnlast)) .gt. 0.05d0) then
c: Move along phase contour to |R1R2|=1 contour:
         klast=k0
         dlnlast=rr0(2,4)
         dmag_left=-dreal(lnlast)
         ph_des0=dimag(lnlast)
         call traj_phase(ph_des0,dmag_left,klast,lnlast,
     .      dlnlast,k0,rr0,.05d0,0.19635d0)
         if(jjfail.gt.0) return
         if(iifail .eq. 1) then
            print *,'Failure in traj_phase from fix_path: ',nmode,k0,rr0
            call stop_run
            return
         endif
         if(iidiag .ge. 2) print *,'traj_phase from fix_path: ',
     .      k0,(rr0(j,4),j=1,2)
         call sheet_look(0,branch)
         if(branch .ne. (0.d0,0.d0)) then
            delk=k0 - branch
            if(dreal(delk) .lt. 0.d0 .or. dimag(delk) .lt. 0.d0) then
               k0=kcut0
               call sheet_init(k0,1,iish,iish_ref)
      print *,'Calling cut_cross2 after bad traj_phase: ',k0/kw0
               call r1r2_calc(k0,rr0,2,0,jjfail)
               if(jjfail .gt. 0) return
               call cut_cross2(k0,rr0)
            endif
         endif
      endif
      if(iileft .eq. 1 .and. dimag(rr0(2,4)) .ge. 0.d0) then
c: Found contour, but heading to right in complex k-plane.  Must find
c: contour crossing farther up:
         call cut_cross(k0,rr0)
cc       print *,'cut_cross: ',k0,(rr0(j,4),j=1,2)
         if(iidone .eq. 1) then
            nmode=nmode - 1
            print *,'no crossing found '
            return
         elseif(iifail .eq. 1) then
            nmode=nmode - 1
            print *,'Failure in cut_cross from fix_path: ',
     .         nmode,k0,rr0
            call stop_run
            return
         endif
      endif
      if(iimt .ne. 0) call mode_traj(k0,rr0,0)
c: Set iimst so that we go to right (backing up) if phase is < pi/2:
      iimst=-1
c
      return
      end
ccc
      subroutine cut_cross(kmid,rrmid)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 kcub,ntry,nwt,jpt,jpt2,iilo
      complex*16 kmid,rrmid(3,4),khi,rrhi(3,4)
      real*8 kilo,kihi,maglo,maghi,dmaglo,dmaghi,c1,c2,c3,c4,delk,kx,
     .   k_real,lnmid,kimid,magmid,dmagmid,wtlo(9),wthi(9),kimid2
      data nwt/9/
      data wtlo/.5,.75,.25,.875,.125,.9375,.0625,.01,.99/
      data wthi/.5,.25,.75,.125,.875,.0625,.9375,.99,.01/
c
c: Found contour, but heading to right in complex k-plane.  Must find
c: contour crossing farther up:
      k_real=dreal(kmid)
      call ki_copy(kmid,rrmid,kilo,maglo,dmaglo)
c: Make sure maglo is not slightly positive (cubic root won't work):
      if(maglo .gt. 0.d0) maglo=-1.d-10
c: Find maximum imaginary value for k of interest:
      if(rmin .ge. 999.) then
c: Use .1 km as arbitrary rmin here:
         khi=dcmplx(k_real,dble(db_cut)/(8685.9*(.1)))
      elseif(kim_max .lt. 1.d100) then
         khi=dcmplx(k_real,kim_max)
      else
         khi=dcmplx(k_real,dkim)
      endif
      if(dimag(khi) .le. kilo) then
c: Return if higher value of ki not larger than starting value:
         iidone=1
         return
      endif
c
      call r1r2_calc(khi,rrhi,2,0,jjfail)
      if(jjfail .gt. 0) return
      call ki_copy(khi,rrhi,kihi,maghi,dmaghi)
c
      if(maghi .lt. 0.d0) then
c: Value is still negative, so no crossing
         if(dmaghi .lt. 0.d0) then
c: But derivative indicates there could be a maximum with positive value:
c: Sample at several points between, searching for positive values
c: or positive derivatives:
            do jpt=1,nwt
               kimid=wtlo(jpt)*kilo + wthi(jpt)*kihi
               kmid=dcmplx(k_real,kimid)
               call r1r2_calc(kmid,rrmid,2,0,jjfail)
               if(jjfail .gt. 0) return
               call ki_copy(kmid,rrmid,kimid,magmid,dmagmid)
               if(magmid .gt. 0.d0) then
                  call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
cc    print *,'cut_cross found pos value '
                  goto 20
               elseif(dmagmid .gt. 0.d0) then
c: Found positive derivative between kimid and kihi, but still 
c: need positive value:
                  do jpt2=1,nwt
                     kimid2=wtlo(jpt2)*kimid + wthi(jpt2)*kihi
                     kmid=dcmplx(k_real,kimid2)
                     call r1r2_calc(kmid,rrmid,2,0,jjfail)
                     if(jjfail .gt. 0) return
                     call ki_copy(kmid,rrmid,kimid,magmid,dmagmid)
                     if(magmid .gt. 0.d0) then
                        call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
cc    print *,'cut_cross found pos value after pos deriv'
                        goto 20
                     endif
                  enddo
cc    print *,'cut_cross found pos deriv, but not pos value'
                  iidone=1
                  return
               endif
            enddo
            iidone=1
            return
         else
            iidone=1
            return
         endif
      endif
c
20    continue
      ntry=0
cccc  iilo=0
      iilo=1
10    continue
      ntry=ntry + 1
      if(iilo .eq. 0) then
c: On first try, try to get away from branch point:
         kmid=dcmplx(k_real,0.999*kilo+.001*kihi)
      elseif(ntry .gt. 30) then
         print *,'Failure in cut_cross'
         iifail=1
         return
      elseif(ntry .lt. 18 .and. mod(ntry,3) .ne. 0) then
         call cub_fit_new(kilo,kihi,maglo,maghi,dmaglo,dmaghi,
     .      c1,c2,c3,c4,delk)
         call cub_root_new(c1,c2,c3,c4,kx,kcub)
         if(kcub .eq. 0) print *,'kcub = 0 from cut_cross'
         kmid=dcmplx(k_real,kilo + kx*delk)
      else
         kmid=dcmplx(k_real,0.5d0*(kilo+kihi))
      endif
      call r1r2_calc(kmid,rrmid,2,0,jjfail)
      if(jjfail .gt. 0) return
      lnmid=dreal(rrmid(1,4))
      if(dabs(lnmid) .lt. .05d0) then
cc       print *,'cut_cross succcessful: ',ntry
         return
      elseif(lnmid .lt. 0.d0) then
         call ki_copy(kmid,rrmid,kilo,maglo,dmaglo)
         iilo=1
      else
         call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
      endif
      goto 10
c
      end
ccc
      subroutine ki_copy(klo,rrlo,kilo,maglo,dmaglo)
c
c: Copies real values from complex values for use in cut_cross.
c
      implicit none
      complex*16 klo,rrlo(3,4)
      real*8 kilo,maglo,dmaglo
c
      kilo=dimag(klo)
      maglo=dreal(rrlo(1,4))
      dmaglo=-dimag(rrlo(2,4))
c
      return
      end
ccc
      subroutine move_up_cut2(k0,rr0,kw,emagmax,iifail,nctot,iidiag,
     .           jjfail)
c
      implicit none
      integer*4 iifail,nctot,ntry,ntry2,iidiag,jjfail
      complex*16 k0,rr0(3,4),kneg,kpos,rrneg(3,4),rrpos(3,4),kneg2
      real*8 kw,emagmax,delki,delk
c
c: Move up in complex plane until |R1*R2| has pos gradient:
cxx   delki=.05d0*kw
      delki=.01d0*kw
      kneg=k0
5     call r1r2_calc(kneg,rrneg,2,0,jjfail)
      if(jjfail .gt. 0) return
      if(dreal(rrneg(1,4)) .gt. 0.d0) then
         delk=.1*delki
         ntry2=0
54       ntry=0
55       ntry=ntry+1
         kneg2=dcmplx(dreal(kneg),dimag(kneg)+delk)
         call r1r2_calc(kneg2,rrneg,2,0,jjfail)
         if(jjfail .gt. 0) return
         if(ntry .gt. 8) then
            kneg=dcmplx(dreal(kneg),dimag(kneg)+.1*delki)
            ntry2=ntry2 + 1
            if(ntry2 .gt. 20) then
               print *,'Failure to pick up path in move_up_cut2: ',
     .            k0/kw,kneg/kw
               iifail=1
               return
            endif
            goto 54
         elseif(dreal(rrneg(1,4)) .gt. 0.d0) then
            delk=0.5d0*delk
            goto 55
         endif
         kneg=kneg2
         if(iidiag .ge. 2) print *,'Pos re(ln) done: ',k0/kw,
     .      kneg/kw,rrneg(1,4),rrneg(2,4)
      endif
c
      kpos=k0
      ntry=0
10    kpos=dcmplx(dreal(kpos),dimag(kpos) + delki)
      ntry=ntry + 1
      if(ntry .gt. 50) then
         print *,'Failure to find kpos in move_up_cut2: ',k0/kw,kpos/kw
         iifail=1
         return
      endif
      call r1r2_calc(kpos,rrpos,2,0,jjfail)
      if(jjfail .gt. 0) return
      if(dreal(rrpos(1,4)) .lt. 0.d0) then
         delki=1.25*delki
         goto 10
      endif
cxx   print *,'move_up2 start: ',kneg/kw,rrneg(1,4),rrneg(2,4),
cxx  .   kpos/kw,rrpos(1,4),rrpos(2,4)
      ntry=0
c
c: We now have kneg with neg Re[ln(r1r2)], kpos with pos Re[ln(r1r2)]:
20    k0=dcmplx(dreal(k0),.5d0*(dimag(kneg)+dimag(kpos)))
      ntry=ntry + 1
      if(ntry .gt. 20) then
         print *,'Unsuccessful move_up_cut2: ',
     .         kneg/kw,rrneg(2,4),kpos/kw,rrpos(2,4)
         iifail=1
         return
      endif
      call r1r2_calc(k0,rr0,2,0,jjfail)
      if(jjfail .gt. 0) return
c: Note that change in Re[f(k)] along imaginary k axis is -IM[df/dk]:
      if(abs(dreal(rr0(1,4))) .lt. emagmax) then
         if(dimag(rrpos(2,4)) .lt. 0.d0 .and.
     .      dimag(rrneg(2,4)) .lt. 0.d0) then
cxx   print *,'cut2 done: ',rrpos(1,4),rrneg(1,4),
cxx  .   kpos/kw,kneg/kw,k0/kw,rr0(1,4),emagmax
            return
         else
            if(iidiag .ge. 2) then
               print *,'emagmax OK, but derivatives wrong sign: ',
     .            kneg/kw,rrneg(2,4),kpos/kw,rrpos(2,4)
            endif
         endif
      endif
      if(dreal(rr0(1,4)) .gt. 0.d0) then
         kpos=k0
         rrpos(1,4)=rr0(1,4)
         rrpos(2,4)=rr0(2,4)
      else
         kneg=k0
         rrneg(1,4)=rr0(1,4)
         rrneg(2,4)=rr0(2,4)
      endif
cxx   print *,'move_up2 loop: ',kneg/kw,rrneg(1,4),rrpos(2,4),
cxx  .   kpos/kw,rrpos(1,4),rrpos(2,4)
      goto 20
c
      end
ccc
      subroutine cut_cross2(kmid,rrmid)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 kcub,ntry,nwt,jpt,jpt2
      complex*16 kmid,rrmid(3,4),khi,rrhi(3,4)
      real*8 kilo,kihi,maglo,maghi,dmaglo,dmaghi,c1,c2,c3,c4,delk,kx,
     .   k_real,lnmid,kimid,magmid,dmagmid,wtlo(9),wthi(9),kimid2
      data nwt/9/
      data wtlo/.5,.75,.25,.875,.125,.9375,.0625,.01,.99/
      data wthi/.5,.25,.75,.125,.875,.0625,.9375,.99,.01/
c
c: Found contour, but heading to right in complex k-plane.  Must find
c: contour crossing farther up:
      k_real=dreal(kmid)
      call ki_copy(kmid,rrmid,kilo,maglo,dmaglo)
c: Make sure maglo is not slightly positive (cubic root won't work):
ccc   if(maglo .gt. 0.d0) maglo=-1.d-10
c: Find maximum imaginary value for k of interest:
      if(rmin .ge. 999.) then
c: Use .1 km as arbitrary rmin here:
         khi=dcmplx(k_real,dble(db_cut)/(8685.9*(.1)))
      elseif(kim_max .lt. 1.d100) then
         khi=dcmplx(k_real,kim_max)
      else
         khi=dcmplx(k_real,dkim)
      endif
      if(dimag(khi) .le. kilo) then
c: Return if higher value of ki not larger than starting value:
         iidone=1
         return
      endif
c
      call r1r2_calc(khi,rrhi,2,0,jjfail)
      if(jjfail .gt. 0) return
      call ki_copy(khi,rrhi,kihi,maghi,dmaghi)
c
      if(maghi .lt. 0.d0) then
c: Value is still negative, so no crossing
         if(dmaghi .lt. 0.d0) then
c: But derivative indicates there could be a maximum with positive value:
c: Sample at several points between, searching for positive values
c: or positive derivatives:
            do jpt=1,nwt
               kimid=wtlo(jpt)*kilo + wthi(jpt)*kihi
               kmid=dcmplx(k_real,kimid)
               call r1r2_calc(kmid,rrmid,2,0,jjfail)
               if(jjfail .gt. 0) return
               call ki_copy(kmid,rrmid,kimid,magmid,dmagmid)
               if(magmid .gt. 0.d0) then
                  call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
cc    print *,'cut_cross found pos value '
                  goto 20
               elseif(dmagmid .gt. 0.d0) then
c: Found positive derivative between kimid and kihi, but still 
c: need positive value:
                  do jpt2=1,nwt
                     kimid2=wtlo(jpt2)*kimid + wthi(jpt2)*kihi
                     kmid=dcmplx(k_real,kimid2)
                     call r1r2_calc(kmid,rrmid,2,0,jjfail)
                     if(jjfail .gt. 0) return
                     call ki_copy(kmid,rrmid,kimid,magmid,dmagmid)
                     if(magmid .gt. 0.d0) then
                        call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
cc    print *,'cut_cross found pos value after pos deriv'
                        goto 20
                     endif
                  enddo
cc    print *,'cut_cross found pos deriv, but not pos value'
                  iidone=1
                  return
               endif
            enddo
            iidone=1
            return
         else
            iidone=1
            return
         endif
      endif
c
20    continue
      ntry=0
10    continue
      ntry=ntry + 1
      if(ntry .eq. 1) then
c: On first try, try to get away from branch point:
         kmid=dcmplx(k_real,0.98*kilo+.02*kihi)
      elseif(ntry .gt. 25) then
         print *,'Failure in cut_cross'
         iifail=1
         return
      elseif(mod(ntry,3) .ne. 0) then
         call cub_fit_new(kilo,kihi,maglo,maghi,dmaglo,dmaghi,
     .      c1,c2,c3,c4,delk)
         call cub_root_new(c1,c2,c3,c4,kx,kcub)
         if(kcub .eq. 0) print *,'kcub = 0 from cut_cross'
         kmid=dcmplx(k_real,kilo + kx*delk)
      else
         kmid=dcmplx(k_real,0.5d0*(kilo+kihi))
      endif
      call r1r2_calc(kmid,rrmid,2,0,jjfail)
      if(jjfail .gt. 0) return
      lnmid=dreal(rrmid(1,4))
      if(dabs(lnmid) .lt. .05d0) then
cc       print *,'cut_cross succcessful: ',ntry
         return
      elseif(lnmid .lt. 0.d0) then
         call ki_copy(kmid,rrmid,kilo,maglo,dmaglo)
      else
         call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
      endif
      goto 10
c
      end
