      subroutine mode_findm
c
c: Finds eigenvalues and computes mode functions.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
      integer*4 iidone
      complex*16 k0,rr0(3,4),rr(3,4)
c
      nmode=0
      call freq_init
      kn(0)=xkref
c: Move around unit circle in complex R1R2-plane in CCW direction:
      iiccw=1
      nm_miss=0
      nhigh=0
      iilk=0
      iifail=0
c
c: Flags to change sheets when crossing Pekeris branch cuts:
      iich=1
c
      phfac0=dmax1(2.d0,dble(phfac))
c: Set ph_step, the phase step initially used along the |R1R2|=1 contour:
      ph_step=twpie/nint(phfac0)
c: Set mag_step, the magnitude step initially used along the phase contour:
      mag_step=0.1d0
c
      k0=(1.d0-1.d-8)*xkref
      call sheet_init(k0,1,iish,iish_ref)
      call fix_path(k0,rr0,iidone)
c
      iidone=0
      if(iimt .ne. 0) call mode_traj_bp(k0,rr0)
c
      iimst=0
      do while(iidone .eq. 0)
         nmode=nmode+1
         ncalc(nmode)=0
         call eig_findm(k0,rr0,kn(nmode),rr,iidone)
         if(jjfail.eq.1) return
         if(iifail .eq. 1) return
      enddo
c
      return
      end
