      subroutine clean_array
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
      include 'lab_com'
      include 'debug_com'

      integer j,jj,izero
      real*4 zero
      real*8 dzero
      complex*16 dczero
      data dczero/(0.d0,0.d0)/
      data dzero/0.d0/
      data zero/0.0/
      data izero/0/

cpln      nzs=0
cpln      nrec=0
cpln      nzmf=0
cpln      nsrc=0
      nzsr=0
      nlay=izero
      nlayb=izero
      nlayt=izero
      nduct=izero
      nzmxtot=izero
      npt=0
      nmin_dln=izero
      kduct=izero
      kduct0=izero
      isvmin=izero
      nsvmin=izero
      jsurf=izero
      jobot=izero
      ii_xi_mv(1)=izero
      ii_xi_mv(2)=izero
c This has caused me a lot of problems
c Has to be reset for multible calls of bb_brute
      cphlo=dzero
c
c Variables used in traj_sdp
c
      k_sdp=dczero
      ln_sdp=dczero
      dln_sdp=dczero
      k_spt=dczero
      ln_spt=dczero
      dln_spt=dczero
c
      do j=1,NLMAX
         zduct(j)=dzero
         isp(j)=izero
         iss(j)=izero
         iiww(j)=izero
         zdep(j)=dzero
cpln         nzmx(j)=izero
cpln         jzmx(j)=izero
         do jj=1,2
            aisoln(jj,j)=izero
         end do
      end do

      do j=1,NSRMAX
         jsrmx(j)=izero
         jsr2j(j)=izero
         mx_m(j)=izero
cpln         zmx(j)=dzero
cpln         zrec(j)=dzero
cpln         zsrc(j)=dzero
         zsr(j)=dzero
cpln         phix(j)=dczero
cpln         dphix(j)=dczero
      end do
      
      do j=1,5
         do jj=1,2
            jflu(jj,j)=izero
         end do
      end do
      
      do j=1,NDMAX
         mzduct(j)=izero
         zpeak(j)=dzero
         do jj=1,5
            jduct(jj,j)=izero
         end do
      end do
      zpeak(NDMAX+1)=dzero

c      do j=1,32
c         jmin_dln(j)=izero
c         dln_deep(j)=zero
c      end do
c
c      do j=1,150
c         k_cont(j)=dczero
c         ln_cont(j)=dczero
c         dln_cont(j)=dczero
c      end do
c
      return
      end

