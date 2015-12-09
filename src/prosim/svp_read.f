c:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************
c: *******************************
c: *     REVISION (1996):        *
c: *         E M G  GROUP        *
c: *     S A C L A N T C E N     *
c: *******************************
      subroutine svp_read
c
      implicit none
      include 'Parms_com'
      include 'i_o_svp_com'
      include 'i_o_opt_com'
      include 'i_o_1b_com'
      include 'sector_env_com'
      include 'depth_com'
      common /out_com3/ nzsr, nlay, NSCTOR,wclmn
      integer*4 nzsr,nlay, NSCTOR,wclmn
      common /alphafrq/ attfq
      logical attfq
c: Local variables:
c      integer*2 jc
c      integer*4 jd,nshft,templuinp
      integer*4 templuinp
c      integer*4 nline,j,j1,j2,iiblug(-1:4),iierr,ndel,nr
      integer*4 nline,j,iierr
c      real*8 zdel(NLMAX),hbb(NLMAX),geobb(2,5,NLMAX),bpbb(2,NLMAX)
      character*64 eline
c fbv
      integer*4 ISECT
cfmc
       integer  MINMOD, MODCUT
c       real*8 z11_fmc, z12_fmc, c11_fmc, c12_fmc

       real*8 dummy1, dummy2, seclen
c       real*8 h0_fmc, h1_fmc, dummy1, dummy2, NWORD, seclen
cc       REAL*8 H0BEG, H1BEG
c       real*8 R1_FMC, B1_FMC
c       real*8 R2_FMC, B2_FMC, C2_FMC, B2S_FMC, C2S_FMC
       real*8 DSVPDH,dsvpdh_tol

cc       COMMON /NN/  H0BEG, H1BEG, NFREQ_FMC, MSP_FMC

      data eline/'INVALID INPUT IN SVP FILE: '/
c      data iiblug/0,0,0,2,2,2/
      data dsvpdh_tol/0.01/
c
c      write(6,*)'FILE UNIT NO: ',luinp
c      open(luinp, err=500,status= 'OLD')
      templuinp = luinp
      iierr=0
      nline=0
      svp_ver= 1.0
      rho_svp= 1.0
cpln      attfq=.true.
      attfq=.false.
      alpha_svp= 0.0
c
      svp_title(1:64)=' '
c     
      read(luinp,*)nsect
c
      READ(luinp,*,ERR=1525) MINMOD, cphmax, MODCUT

C REVIEW 1525 CONTINUE
 1525 CONTINUE

      do isect=1,nsect
c: Sound speed profile in ocean:
         READ(luinp,*,err=480)r_h0(isect), dummy1, dummy2, seclen
c
         write(6,*)'Water depth: ',r_h0(isect)
c     
c The distance from the source is stored in a common variable:
         secleng(ISECT) = seclen*1000.0
c
         write(6,*)'Sector length: ',secleng(isect)
c
cc      H0BEG= h0_fmc
         nline = 1
         read(luinp,*,end=510,err=510)R_Z0(1,ISECT),R_C0(1,ISECT)
         write(6,*)'Sound Speed in water at sector: ',isect
         write(6,*)R_Z0(1,ISECT),R_C0(1,ISECT)
         call check_val_r8(r_z0(1,isect),0.d0,0.d0,eline,27,
     .        nline,'zsvp(1) MUST BE 0',17,iierr)
         call check_val_r8(r_c0(1,isect),1.d-10,1.d+10,eline,27,
     .        nline,'csvp(1)',7,iierr)
         do j=2,nlmax
            read(luinp,*,end=510,err=510)R_Z0(j,ISECT),R_C0(j,ISECT)
c     
            write(6,*)R_Z0(j,ISECT),R_C0(j,ISECT)
c
c zsvp(1),csvp(1)
            nline= nline+1
cc fbv         if(j.eq.1) then
cc fbv            call check_val_r8(r_z0(j,isect),0.d0,0.d0,eline,27,
cc fbv     .        nline,'zsvp(1) MUST BE 0',17,iierr)
cc fbv         else
            call check_val_r8(r_z0(j,isect),r_z0(j-1,isect),1.d10,
     .           eline,27,nline,'zsvp(1) MUST BE 0',17,iierr)
cc fbv         end if
            call check_val_r8(r_c0(j,isect),1.d-10,1.d+10,eline,27,
     .           nline,'csvp(1)',7,iierr)
c
cc fbv         if(j.gt.1) then
            DSVPDH=(R_C0(j,ISECT)-R_C0(j-1,ISECT))/
     .           (R_Z0(j,ISECT)-R_Z0(j-1,ISECT))
            if((dabs(DSVPDH).lt.dsvpdh_tol).and.
     .           (dabs(DSVPDH).gt.0.d0)) then
               write(6,*)'sound speed gradient too small ',
     .              'between layer ',j-1,' and ',j
               write(6,*)'Gradient: ',DSVPDH,' m/s/m'
cc fbv               stop
            end if
cc fbv         end if
            if( r_z0(j,isect) .eq. r_h0(isect) )   go to 1200
         end do
         print *, ' sub svp_read, too many points in svp '
         stop
 1200    continue
c
         nsvp= j
         write(*,*)'Number of sound speed points', nsvp
         R_ND0(ISECT)=nsvp
         secz(isect)=r_h0(isect)
c
c: Bottom layering:
cpln ONLY IVWK
cpln         read(luinp,*)  R_H1(isect), R_R1(isect)
         read(luinp,*)  R_H1(isect), R_R1(isect), R_BETA(1,isect)
         write(6,*)' Sediment values'
         write(6,*) R_H1(isect), R_R1(isect), R_BETA(1,isect)
         if( R_H1(isect) .gt. 0.0 ) then
            do j= 1, NLMAX
cpln ONLY IVWK
cpln               read(luinp,*,end=510,err=510) R_Z1(j,isect),
cpln     .              R_C1(j,isect), R_BETA(j,isect)
               read(luinp,*,end=510,err=510) R_Z1(j,isect),
     .              R_C1(j,isect)
               write(6,*)R_Z1(j,isect),R_C1(j,isect)
               if( r_z1(j,isect) .eq. r_h1(isect) )   go to 1400
            enddo
            print *,' sub opt_read, e^ successo un casino '
            stop
 1400       continue
            nlayb=j
         else
            nlayb=0
         end if
cpln ONLY ISO speed bottom layering
cpln ONLY ASCOT01
cpln         if(nlayb.eq.2)
cpln     .     R_C1(2,isect)=R_C1(1,isect)
         R_ND1(ISECT)=nlayb

c: Lower halfspace:
cpln ONLY IVWK
cpln         read(luinp,*)R_R2(isect),R_BETA(nlayb+1,isect),
cpln     .                R_C2(isect)
         read(luinp,*)  R_R2(isect), R_BETA(2,isect), R_C2(isect)
         write(6,*)' bottom values ' 
         write(6,*)R_R2(isect), R_BETA(2,isect), R_C2(isect)
c     
cpln ONLY IVWK
cpln         read(luinp,*)  R_BETA(nlayb+2,isect),  R_C2S(isect)
cpln         write(6,*)  R_BETA(nlayb+2,isect),  R_C2S(isect)
         read(luinp,*)  R_BETA(3,isect),  R_C2S(isect)
         write(6,*)  R_BETA(3,isect),  R_C2S(isect)
c
      end do
 480  continue
      ISECT=ISECT-1
      write(6,*)'ISECT: ',ISECT
c     write(6,*)'AFTER 480 continue'
c
      NSCTOR = ISECT
c
      if(iierr .ne. 0) stop
c
      do ISECT = NSCTOR, 2, -1
         secleng(ISECT)=secleng(ISECT)-secleng(ISECT-1)
      enddo
      secz(NSCTOR+1)=secz(NSCTOR)
c
      f_max=fmax
      if (multcw) then
         f_max=fcw(1)
         do j=1,nfcw
            f_max=dmax1(dble(fcw(j)),f_max)
         enddo
      end if
c     
      return
c 500  print *,'Error opening input file'
c      stop
 510  print *,'Endo or error reading input file at line ',nline
      stop
      end
