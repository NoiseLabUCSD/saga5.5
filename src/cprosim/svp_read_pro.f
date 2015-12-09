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
      subroutine svp_read_pro
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'debug_com'
      include 'sector_env_com'
      include 'depth_com'
      integer*4 ISECT,NSCTOR,NSECT
c
c: Local variables:
c      integer*2 jc
c      integer*4 jd,nshft,templuinp
      integer*4 templuinp
c      integer*4 nline,j,j1,j2,iiblug(-1:4),iierr,ndel,nr
      integer*4 nline,j,iierr,nfcr
c      real*8 zdel(NLMAX),hbb(NLMAX),geobb(2,5,NLMAX),bpbb(2,NLMAX)
c fbv
      integer  MINMOD, MAXMOD, MODCUT
      real*8 dummy1, dummy2, seclen
      character*64 eline

cc       COMMON /NN/  H0BEG, H1BEG, NFREQ_FMC, MSP_FMC

      data eline/'INVALID INPUT IN SVP FILE: '/
c      data iiblug/0,0,0,2,2,2/
c
c
      write(lusvp,95)
95    format(/'### SVP FILE INFORMATION ###')
c
c      open(luinp, err=500,status= 'OLD')
c      open(luinp, file='cprosim_in',err=500,status= 'OLD')
      templuinp = luinp
      iierr=0
      nline=0
      svp_ver= 2.0
      rho_svp= 1.0
      alpha_svp= 0.0
      iidebug=0
c
c

      nfcr = nfcw
      READ(luinp,*)NSECT
      READ(luinp,*,ERR=1525) MINMOD, MAXMOD, MODCUT
c      READ(luinp,*,ERR=1525) MINMOD, cphmax, MODCUT
      write(6,*)'IICW: ',iicw
      write(6,*)'NSECT: ',NSECT
      write(6,*)'MINMOD, cphmax, MODCUT: ',MINMOD, cphmax, MODCUT
c
C REVIEW 1525 CONTINUE
 1525 CONTINUE

 490  CONTINUE

c: Sound speed profile in ocean:
      DO ISECT=1,NSECT
         READ(luinp,*)   r_h0(isect), dummy1, dummy2, seclen
c
c
c The distance from the source is stored in a common variable:
         secleng(ISECT) = seclen*1000.0
c
         if(iiwrite .eq.0) then
            write(6,*)'Water depth: ',r_h0(isect)
            write(6,*)'Sector length: ',secleng(isect)
         end if
c
cc      H0BEG= h0_fmc
         nline = 1
         read(luinp,*,end=510,err=510)R_Z0(1,ISECT),R_C0(1,ISECT)
         call check_val_r8(r_z0(1,isect),0.d0,0.d0,eline,27,
     .        nline,'zsvp(1) MUST BE 0',17,iierr)
         call check_val_r8(r_c0(1,isect),1.d-10,1.d+10,eline,27,
     .        nline,'csvp(1)',7,iierr)
         if(iidebug .eq. 1) then
            write(6,*)'Sound Speed in water at sector: ',isect
            write(6,*)R_Z0(1,ISECT),R_C0(1,ISECT)
         end if
         do j=2,nlmax
            read(luinp,*,end=510,err=510)R_Z0(j,ISECT),R_C0(j,ISECT)
c
            if(iidebug .eq. 1)
     .           write(6,*)R_Z0(j,ISECT),R_C0(j,ISECT)
c     
c     nline= nline+1
            call check_val_r8(r_z0(j,isect),r_z0(j-1,isect),1.d10,
     .           eline,27,nline,'zsvp(1) MUST BE 0',17,iierr)
            call check_val_r8(r_c0(j,isect),1.d-10,1.d+10,eline,27,
     .           nline,'csvp(1)',7,iierr)
c     
            if( r_z0(j,isect) .eq. r_h0(isect) )   go to 1200
         end do
         print *, ' sub svp_read, too many points in svp '
         stop
 1200    continue
c     
         nsvp= j
c     pln      write(*,*)'Number of sound speed points', nsvp
         R_ND0(ISECT)=nsvp
         secz(isect)=r_h0(isect)
c     
c     : Bottom layering:
         read(luinp,*)  R_H1(isect), R_R1(isect),R_BETA(1,isect)
         if(iidebug .eq. 1) then
            write(6,*)' Sediment values'
            write(6,*) R_H1(isect), R_R1(isect),R_BETA(1,isect)
         end if
         if( R_H1(isect) .gt. 0.0 ) then
            do j= 1, NLMAX
cpln Read attenuation in bottom:
cpln               read(luinp,*,end=510,err=510) R_Z1(j,isect),
cpln     .              R_C1(j,isect),R_BETA(j,isect)
cpln     .              write(6,*) R_Z1(j,isect), R_C1(j,isect), 
cpln     .                         R_BETA(j,isect)
               read(luinp,*,end=510,err=510) R_Z1(j,isect),
     .              R_C1(j,isect)
               if(iidebug .eq. 1)
cpln: Write attenuation in bottom
cpln     .              write(6,*) R_Z1(j,isect), R_C1(j,isect),
cpln     .                         R_BETA(j,isect)
     .              write(6,*) R_Z1(j,isect), R_C1(j,isect)
               if( r_z1(j,isect) .eq. r_h1(isect) )   go to 1400
            enddo
            print *,' sub opt_read, e^ successo un ...... '
            stop
 1400       continue
            nlayb=j
         else
            nlayb=0
         end if
cpln ONLY ASCOT01
cpln         if(nlayb.eq.2)
cpln     .     R_C1(2,isect)=R_C1(1,isect)
         R_ND1(ISECT)=nlayb

c: Lower halfspace:
cpln: Attenuation in sediment
cpln         read(luinp,*)  R_R2(isect), R_BETA(nlayb+1,isect), R_C2(isect)
cpln         read(luinp,*)  R_BETA(nlayb+2,isect),  R_C2S(isect)
         read(luinp,*)  R_R2(isect), R_BETA(2,isect), R_C2(isect)
         read(luinp,*)  R_BETA(3,isect),  R_C2S(isect)
c
         if(iidebug .eq. 1) then
            write(6,*)'Sub-bottom values ' 
cpln: Attenuation in sediment
cpln            write(6,*)R_R2(isect), R_BETA(nlayb+1,isect), R_C2(isect)
cpln            write(6,*)  R_BETA(nlayb+2,isect),  R_C2S(isect)
            write(6,*)R_R2(isect), R_BETA(2,isect), R_C2(isect)
            write(6,*)  R_BETA(3,isect),  R_C2S(isect)
         end if
      end do
c
c Use this line for DEC/HP/INTEL (unix) workstations
 480  continue
      ISECT=ISECT-1
      write(6,*)'ISECT: ',ISECT
c
c Use this line for SUN (sunos) workstations
c 480  continue
c
      NSCTOR = ISECT
c
      if(iierr .ne. 0) stop
c
      do ISECT = NSCTOR, 2, -1
       secleng(ISECT) = secleng(ISECT) - secleng(ISECT-1)
      enddo
      secz(NSCTOR+1)=secz(NSCTOR)
c
      f_max=fmax
      if (iicw .eq. 1) then
         f_max=fcw(1)
         do j=1,nfcw
            f_max=dmax1(dble(fcw(j)),f_max)
         enddo
         call hpsort_r4(nfcr,fcw)
         if(iiwrite .eq. 0)
     &      write(6,*)'Frequency band: ',fcw(1),fcw(nfcr)
      end if
c
      return
 500  print *,'Error opening input file'
      stop
 510  print *,'Endo or error reading input file at line ',nline
      stop
      end
