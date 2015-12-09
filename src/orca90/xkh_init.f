      subroutine xkh_init(ndv)
c
      use parms_com
      use i_o_com
      use gen_com
c
      integer ndv,j,ii,jj,jp1,jm1,ii2
      complex*16 ik,xbx
c
c: Keep track of previous ratio in case we need to back up in eig_find:
      xkhratp=xkhrat
      call iish_xfer(iish,iish_ref,iish0,iishr0)
c: Variable used to check for sheet changes (use kw0, which does not change
c: when the duct changes):
      xkhrat=xkh/kw0
      xkhsq=xkh*xkh
c
      ik=dcmplx(-dimag(xkh),dreal(xkh))
      ikcon(1)=ik
      ikcon(2)=(0.d0,1.d0)
      ikcon(3)=(0.d0,0.d0)
c
      if(isp(nlay) .eq. 0) then
c: For Airy halfspaces, work outwards from reference depth:
c: Reference layer:
         j=nsvmin
         ii=isvmin
         ii2=3-isvmin
         call sheet(xksq(ii,j),xkhsq,xk(ii,j),xkh,w,gami(1,ii,j),
     .      ndv,iich_ref,iish_ref(1),xkhratp,xkhrat,xkrat_ref(1))
         if(isp(j) .eq. 0) then
            call gami_calc(xkh,xkhsq,xk(ii2,j),xksq(ii2,j),w,
     .         gami(1,ii2,j),ndv)
         else
            do jj=1,ndv
               gami(jj,ii2,j)=gami(jj,ii,j)
            enddo
         endif
c
c: Above reference layer:
         do j=nsvmin-1,2,-1
            jp1=j+1
            call gami_lay(mm(j),isp(j),xkh,xkhsq,gami(1,2,j),
     .         gami(1,1,jp1),gami(1,1,j),xk(2,j),xksq(2,j),
     .         xk(1,j),xksq(1,j),w,ndv)
            if(iisol(j) .eq. 1) then
               call gami_lay(mm(j),iss(j),xkh,xkhsq,beti(1,2,j),
     .            beti(1,1,jp1),beti(1,1,j),xb(2,j),xbsq(2,j),
     .            xb(1,j),xbsq(1,j),w,ndv)
            endif
         enddo
c: Below reference layer:
         do j=nsvmin+1,nlay-1
            jm1=j-1
            call gami_lay(mm(jm1),isp(j),xkh,xkhsq,gami(1,1,j),
     .         gami(1,2,jm1),gami(1,2,j),xk(1,j),xksq(1,j),
     .         xk(2,j),xksq(2,j),w,ndv)
            if(iisol(j) .eq. 1) then
               call gami_lay(mm(jm1),iss(j),xkh,xkhsq,beti(1,1,j),
     .            beti(1,2,jm1),beti(1,2,j),xb(1,j),xbsq(1,j),
     .            xb(2,j),xbsq(2,j),w,ndv)
            endif
         enddo
c
c: Lower halfspace:
         call gami_hsp(nlay,1,ndv,mm(nlay-1),gami(1,2,nlay-1))
c: Upper halfspace:
         call gami_hsp(1,2,ndv,mm(1),gami(1,1,2))
      else
c: For homogeneous halfspaces, work inwards from halfspaces:
c: Lower halfspace:
         call gami_hsp(nlay,1,ndv,mm(nlay-1),gami(1,2,nlay-1))
c: Upper halfspace:
         call gami_hsp(1,2,ndv,mm(1),gami(1,1,2))
c: Above reference layer:
         do j=2,nsvmin-1
            jm1=j-1
            call gami_lay(mm(jm1),isp(j),xkh,xkhsq,gami(1,1,j),
     .         gami(1,2,jm1),gami(1,2,j),xk(1,j),xksq(1,j),xk(2,j),
     .         xksq(2,j),w,ndv)
            if(iisol(j) .eq. 1) then
               call gami_lay(mm(jm1),iss(j),xkh,xkhsq,beti(1,1,j),
     .            beti(1,2,jm1),beti(1,2,j),xb(1,j),xbsq(1,j),
     .            xb(2,j),xbsq(2,j),w,ndv)
            endif
         enddo
c: Below reference layer:
         do j=nlay-1,nsvmin+1,-1
            jp1=j+1
            call gami_lay(mm(j),isp(j),xkh,xkhsq,gami(1,2,j),
     .         gami(1,1,jp1),gami(1,1,j),xk(2,j),xksq(2,j),
     .         xk(1,j),xksq(1,j),w,ndv)
            if(iisol(j) .eq. 1) then
               call gami_lay(mm(j),iss(j),xkh,xkhsq,beti(1,2,j),
     .            beti(1,1,jp1),beti(1,1,j),xb(2,j),xbsq(2,j),
     .            xb(1,j),xbsq(1,j),w,ndv)
            endif
         enddo
c: Reference layer:
         if(nsvmin .ne. nlay .and. nsvmin .ne. 1) then
            j=nsvmin
            ii=isvmin
            ii2=3-isvmin
            call sheet_ref(xksq(ii,j),xkhsq,xk(ii,j),xkh,w,
     .         gami(1,ii,j),ndv,iich_ref,iish_ref(1),xkhratp,xkhrat,
     .         xkrat_ref(1))
            if(isp(j) .eq. 0) then
               call gami_calc(xkh,xkhsq,xk(ii2,j),xksq(ii2,j),w,
     .            gami(1,ii2,j),ndv)
            else
               do jj=1,ndv
                  gami(jj,ii2,j)=gami(jj,ii,j)
               enddo
            endif
         endif
      endif
c
c: Vertical wavenumber at reference depth:
      gamiref=gami(1,isvmin,nsvmin)
c
      do j=1,jsol(2,1)-1
         if(mm(j) .eq. 1) then
            jp1=j + 1
            call pquv_calc(rhorat(j),w,ik,xbsqinv(2,j),xbsqinv(1,jp1),
     .         xkh,Pcon(1,j),Qcon(1,j),Ucon(1,j),Vcon(1,j),ndv)
         endif
      enddo
      do j=jsol(1,1),nlay-1
         if(mm(j) .eq. 1) then
            jp1=j + 1
            call pquv_calc(rhorat(j),w,ik,xbsqinv(1,jp1),xbsqinv(2,j),
     .         xkh,Pcon(1,j),Qcon(1,j),Ucon(1,j),Vcon(1,j),ndv)
         endif
      enddo
c
      do ii=1,2
c: Alay,Blay(1:2,ii),ikcon are elements and derivatives 
c: of the T^(-1) matrix at the top of the first solid layer 
c: in bottom (ii=1) and top (ii=2) layering.
         if(allf(ii) .eq. 0) then
            xbx=xbsqinv(ii,jsol(ii,1))
            Alay(1,ii)=-1.d0 + 2.d0*xkhsq*xbx
            Blay(1,ii)=-2.d0*ik*xbx
            if(ndv .ge. 2) then
               Alay(2,ii)=4.d0*xkh*xbx
               Blay(2,ii)=(0.,-2.d0)*xbx
               if(ndv .ge. 3) then
                  Alay(3,ii)=-4.d0*xkhsq*xbx/w
                  Blay(3,ii)=-2.d0*Blay(1,ii)/w
               endif
            endif
         endif
      enddo
c
      return
      end
ccc
      subroutine gami_hsp(j,ii,ndv,mmx,gamip)
c
      use parms_com
      use i_o_com
      use gen_com
c
      integer j,ii,ndv,mmx,jj
      complex*16 gamip(3)
c
      if(isp(j) .eq. 1) then
c: Homogeneous halfspace:
         call sheet(xksq(ii,j),xkhsq,xk(ii,j),xkh,w,gami(1,ii,j),
     .      ndv,iich,iish(ii,1),xkhratp,xkhrat,xkrat(ii,1))
      else
c: Airy halfspace:
         if(j .ne. nsvmin) then
            if(mmx .eq. 0) then
               do jj=1,ndv
                  gami(jj,ii,j)=gamip(jj)
               enddo
            else
               call gami_calc(xkh,xkhsq,xk(ii,j),xksq(ii,j),w,
     .            gami(1,ii,j),ndv)
            endif
         endif
      endif
      if(iisol(j) .eq. 1) then
         if(iss(j) .eq. 1) then
            call sheet(xbsq(ii,j),xkhsq,xb(ii,j),xkh,w,beti(1,ii,j),
     .         ndv,iich,iish(ii,2),xkhratp,xkhrat,xkrat(ii,2))
         else
            call gami_calc(xkh,xkhsq,xb(ii,j),xbsq(ii,j),w,
     .         beti(1,ii,j),ndv)
         endif
      endif
c
      return
      end
ccc
      subroutine gami_lay(mmx,ispx,xkh,xkhsq,gami1,gamip,gami2,
     .   xk1,xk1sq,xk2,xk2sq,w,ndv)
c
      implicit none
      integer*4 mmx,ispx,ndv,jj
      complex*16 xkh,xkhsq,gami1(3),gamip(3),gami2(3),xk1,xk1sq,
     .   xk2,xk2sq
      real*8 w
c
      if(mmx .eq. 0) then
         do jj=1,ndv
            gami1(jj)=gamip(jj)
         enddo
      else
         call gami_calc(xkh,xkhsq,xk1,xk1sq,w,gami1,ndv)
      endif
c
      if(ispx .eq. 0) then
         call gami_calc(xkh,xkhsq,xk2,xk2sq,w,gami2,ndv)
      else
         do jj=1,ndv
            gami2(jj)=gami1(jj)
         enddo
      endif
c
      return
      end
ccc
      subroutine sheet_init(k0,ii_ones,iish_in,iish_ref_in)
c
      use parms_com
      use gen_com
      integer*4 ii_ones,iish_in(4),iish_ref_in(2)
      complex*16 k0
c
      xkhratp=k0/kw0
      xkhrat=k0/kw0
      if(ii_ones .eq. 1) then
         iish(1,1)=1
         iish(1,2)=1
         iish(2,1)=1
         iish(2,2)=1
         iish_ref(1)=1
         iish_ref(2)=1
      else
         call iish_xfer(iish_in,iish_ref_in,iish,iish_ref)
      endif
c
      return
      end
ccc
      subroutine gami_calc(xkh,xkhsq,xk,xksq,omega,gami,ndv)
c
c: Computes vertical wavenumber and its derivative w.r.t k using 
c: the Pekeris cut.
      implicit none
      integer ndv
      complex*16 xkh,xkhsq,xk,xksq,gami(3),gamma,gamsq
      real*8 omega
c
      gamsq=xksq - xkhsq
      gamma=cdsqrt(gamsq)
c: Pekeris cut:
      if(dreal(xkh) .gt. dreal(xk)) then
         if(dimag(gamma) .lt. 0.d0) then
            gamma=-gamma
         endif
      endif
c
      gami(1)=dcmplx(-dimag(gamma),dreal(gamma))
      if(ndv .ge. 2) then
         if(gami(1) .ne. (0.d0,0.d0)) then
            gami(2)=xkh/gami(1)
            if(ndv .ge. 3) gami(3)=-xksq/(omega*gami(1))
         else
            gami(2)=1.d100*xkh
            if(ndv .ge. 3) gami(3)=1.d100*xkh
         endif
      endif
c
      return
      end
ccc
      subroutine sheet(xksq,xkhsq,xk,xkh,omega,gami,ndv,iich,iish,
     .   xkhratp,xkhrat,xkrat)
c
c: Computes vertical wavenumber and its derivative w.r.t k using 
c: Pekeris cut.  Sheets are changed if iich=1 and present and previous
c: k-values straddle the branch cut. Positive sheet always taken when
c: iich=-1
      implicit none
      integer ndv,iich,iish
      complex*16 xksq,xkhsq,xk,xkh,gamma,gami(3),gamsq,
     .   xkhratp,xkhrat,xkrat
      real*8 omega
c
      gamsq=xksq - xkhsq
      gamma=cdsqrt(gamsq)
c: Check if present and previous points straddle branch cut:
      if(iich .eq. 1) then
c: Check sheet changing in cos(theta)=k/kw0 space instead of k-space
c: so that frequency changes won't affect the process:
         call cross_cut_pek(xkhratp,xkhrat,xkrat,iish)
      endif
c: Get on positive sheet of Pekeris branch cut (do EJP only when 
c: k to right of branch point):
      if(dreal(xkhrat) .gt. dreal(xkrat)) then
         if(dimag(gamma) .lt. 0.) then
            gamma=-gamma
         endif
      endif
c
c: Apply current sheet to gamma:
      gamma=iish*gamma
      gami(1)=dcmplx(-dimag(gamma),dreal(gamma))
      if(ndv .ge. 2) then
         if(gami(1) .ne. (0.,0.)) then
            gami(2)=xkh/gami(1)
         else
            gami(2)=1.d100*xkh
         endif
         if(ndv .ge. 3) then
            if(gami(1) .ne. (0.d0,0.d0)) then
               gami(3)=-xksq/(omega*gami(1))
            else
               gami(3)=1.d100*xkh
            endif
         endif
      endif
c
      return
      end
ccc
      subroutine sheet_ref(xksq,xkhsq,xk,xkh,omega,gami,ndv,iich_ref,
     .   iish_ref,xkhratp,xkhrat,xkrat)
c
c: Computes vertical wavenumber and its derivative w.r.t k using 
c: Pekeris cut.  Sheets are changed if iich_ref=1 and present and previous
c: k-values straddle the branch cut.
      implicit none
      integer ndv,iich_ref,iish_ref
      complex*16 xksq,xkhsq,xk,xkh,gamma,gami(3),gamsq,
     .   xkhratp,xkhrat,xkrat
      real*8 omega
c
      gamsq=xksq - xkhsq
      gamma=cdsqrt(gamsq)
c: Check if present and previous points straddle branch cut:
      if(iich_ref .eq. 1) then
         call cross_cut_pek(xkhratp,xkhrat,xkrat,iish_ref)
      endif
c: Get on positive sheet of Pekeris branch cut (do EJP only when 
c: k to right of branch point):
      if(dreal(xkhrat) .gt. dreal(xkrat)) then
         if(dimag(gamma) .lt. 0.d0) then
            gamma=-gamma
         endif
      endif
c
c: Apply current sheet to gamma:
      gamma=iish_ref*gamma
      gami(1)=dcmplx(-dimag(gamma),dreal(gamma))
      if(ndv .ge. 2) then
         if(gami(1) .ne. (0.,0.)) then
            gami(2)=xkh/gami(1)
         else
            gami(2)=1.d100*xkh
         endif
         if(ndv .ge. 3) then
            if(gami(1) .ne. (0.d0,0.d0)) then
               gami(3)=-xksq/(omega*gami(1))
            else
               gami(3)=1.d100*xkh
            endif
         endif
      endif
c
      return
      end
ccc
      subroutine cross_cut_pek(kp,k,kx,iishq)
c
c: Checks to see if the Pekeris branch cut has been crossed.
c: kp,k are previous and current values of xkh, kx is branch point.
      use parms_com
      use i_o_com
      use gen_com
      integer*4 iishq
      complex*16 kp,k,kx,delk
      real*8 ypt
c
      if((dreal(kp) .le. dreal(kx) .and. dreal(k) .gt. dreal(kx)) .or.
     .   (dreal(kp) .gt. dreal(kx) .and. dreal(k) .le. dreal(kx))) then
         delk=k - kp
         ypt=dimag(kp) + (dreal(kx)-dreal(kp))*dimag(delk)/dreal(delk)
         if(ypt .gt. dimag(kx)) then
            iishq=-iishq
            iicut=1
            if(iiccw*dreal(delk) .lt. 0.d0) then
               kcut=k*kw0
            else
               kcut=kp*kw0
            endif
            if(iidiag .ge. 2) then 
               print *,'pek sheet changed: ',-iishq,' to ',iishq,
     .            dreal(kp),dreal(k),dreal(kx),
     .            dimag(kp)*8685.9*kw0,dimag(k)*8685.9*kw0,
     .            dimag(kx)*8685.9*kw0
            endif
         endif
      endif
c
      return
      end
ccc
      subroutine pquv_calc(rhorat,w,ik,xbsqinv1,xbsqinv2,xkh,
     .   P,Q,U,V,ndv)
c
      implicit none
      integer*4 ndv
      complex*16 ik,xbsqinv1,xbsqinv2,xkh,P(3),Q(3),U(3),V(3),F,
     .   i1mp,k_w
      real*8 rhorat,w
c
      F=2.d0*(xbsqinv1 - rhorat*xbsqinv2)
      U(2)=dcmplx(-dimag(F),dreal(F))
      U(1)=xkh*U(2)
      Q(2)=-2.d0*xkh*F
      Q(1)=0.5d0*xkh*Q(2) + 1.d0
      P(2)=-Q(2)
      P(1)=0.5d0*xkh*P(2) + rhorat
      i1mp=dcmplx(0.d0,1.d0)*(1.d0 - P(1))
      V(2)=-ik*P(2) + i1mp
      V(1)=xkh*i1mp
      if(ndv .ge. 3) then
         k_w=-xkh/w
         P(3)=k_w*P(2)
cxx      Q(3)=k_w*Q(2)
         Q(3)=-P(3)
         U(3)=2.d0*k_w*U(2)
         V(3)=-ik*P(3)
      endif
c
      return
      end
ccc
      subroutine iish_xfer(iish_in,iish_ref_in,iish_out,iish_ref_out)
c
      implicit none
      integer*4 iish_in(2,2),iish_out(2,2),iish_ref_in(2),
     .   iish_ref_out(2)
c
      iish_out(1,1)=iish_in(1,1)
      iish_out(2,1)=iish_in(2,1)
      iish_out(1,2)=iish_in(1,2)
      iish_out(2,2)=iish_in(2,2)
      iish_ref_out(1)=iish_ref_in(1)
      iish_ref_out(2)=iish_ref_in(2)
c
      return
      end
ccc
      subroutine iish_code(iish,iish_ref,iicode,ii)
c
c: Codes iish(2,2) into a single integer iicode for storage (ii=1).
c: or decodes iicode into iish(2,2) (ii=-1).
c
      implicit none
      integer*4 iish(4),iish_ref(2),iicode,ii,ipow(4),ipow2(2),j,iic
      data ipow/-1,-2,-4,-8/,ipow2/-16,-32/
c
      if(ii .eq. 1) then
         iicode=min(iish_ref(2),0)*32 + min(iish_ref(1),0)*16 + 
     .      min(iish(4),0)*8 + min(iish(3),0)*4 + 
     .      min(iish(2),0)*2 + min(iish(1),0)
      else
         iish_ref(1)=1
         iish_ref(2)=1
         iish(1)=1
         iish(2)=1
         iish(3)=1
         iish(4)=1
c: Short cut:
         if(iicode .eq. 0) return
         iic=iicode
c: Decode the sheets:
         do j=2,1,-1
            if(iic .le. ipow2(j)) then
               iish_ref(j)=-1
               iic=iic - ipow2(j)
            endif
         enddo
         do j=4,1,-1
            if(iic .le. ipow(j)) then
               iish(j)=-1
               iic=iic - ipow(j)
            endif
         enddo
      endif
c
      return
      end
ccc
      subroutine sheet_look(kkch,branch)
c
      use parms_com
      use i_o_com
      use gen_com
      integer kkch,ii0
      complex*16 branch
c
      branch=(0.d0,0.d0)
      ii0=0
      if(isp(1) .eq. 1) call sheet_ch(iish(1,1),xkbp(1,1),
     .   branch,kkch,ii0)
      if(iss(1) .eq. 1) call sheet_ch(iish(1,2),xkbp(1,2),
     .   branch,kkch,ii0)
      if(isp(nlay) .eq. 1) call sheet_ch(iish(2,1),xkbp(2,1),
     .   branch,kkch,ii0)
      if(iss(nlay) .eq. 1) call sheet_ch(iish(2,2),xkbp(2,2),
     .   branch,kkch,ii0)
c
      return
      end
ccc
      subroutine sheet_ch(iish,xkbp,branch,kkch,ii0)
c
      implicit none
      integer iish,kkch,ii,ii0,jj
      complex*16 xkbp,branch
c
      if(iish .lt. 0) then
c: If more than one negative sheet, choose cut with branch point farthest
c: to the right:
         if(dreal(xkbp) .gt. dreal(branch)) then
            branch=xkbp
            if(kkch .eq. 1) then
               iish=1
c: If more than on -1 sheet, change left ones back to -1:
               if(ii0 .gt. 0) iish=-1
               ii0=ii
            endif
         endif
      endif
c
      return
      end
ccc
      subroutine xkh_backup
c
c: Called from k-plane search routines when we want to back up to previously
c: called value of k.
c
      use parms_com
      use i_o_com
      use gen_com
c
      iish(1,1)=iish0(1,1)
      iish(2,1)=iish0(2,1)
      iish(1,2)=iish0(1,2)
      iish(2,2)=iish0(2,2)
      iish_ref(1)=iishr0(1)
      iish_ref(2)=iishr0(2)
      xkhrat=xkhratp
c
      return
      end
