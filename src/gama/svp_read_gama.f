      subroutine svp_read_gama
c
      implicit integer*4(i-n)
      include 'common/parms_conv'
c
c: Local variables:
      integer*4 luinp,nlmax
c      integer*4 nline,j,j1,j2,iiblug(-1:4),iierr,ndel,nr
      integer*4 nline,j,iierr,nfcr
c      real*8 zdel(NLMAX),hbb(NLMAX),geobb(2,5,NLMAX),bpbb(2,NLMAX)
c fbv
      integer  MINMOD, MAXMOD, MODCUT
      real*8 dummy1, dummy2, dummy3, seclen

      data nlmax/10000/
c
      write(57,95)
95    format(/'### SVP FILE INFORMATION ###')
c
      iidebug=1
      luinp=1
c
      READ(luinp,*)NSECT,NTYPE
      READ(luinp,*,ERR=1525) MINMOD, MAXMOD, MODCUT
c      READ(luinp,*,ERR=1525) MINMOD, cphmax, MODCUT
      write(6,*)'IICW: ',iicw
      write(6,*)'NSECT, NTYPE: ',NSECT,NTYPE
      write(6,*)'MINMOD, cphmax, MODCUT: ',MINMOD, cphmax, MODCUT
c
C REVIEW 1525 CONTINUE
 1525 CONTINUE

 490  CONTINUE

c: Sound speed profile in ocean:
      DO ISECT=1,NSECT
         READ(luinp,*) r_h0(isect), dummy1, dummy2, seclen
c
c
c The distance from the source is stored in a common variable:
c
         if(iiwrite .eq.0) then
            write(6,*)'Water depth: ',r_h0(isect)
         end if
c
cc      H0BEG= h0_fmc
         nline = 1
         read(luinp,*,end=510,err=510)R_Z0(1,ISECT),R_C0(1,ISECT)
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
            if( r_z0(j,isect) .eq. r_h0(isect) )   go to 1200
         end do
         print *, ' sub svp_read, too many points in svp '
         stop
 1200    continue
c     
         nsvp= j
c     pln      write(*,*)'Number of sound speed points', nsvp
         R_ND0(ISECT)=nsvp
c     
c     : Bottom layering:

         if(NTYPE.eq.1) then
           read(luinp,*) R_H1(isect), R_R1(1,isect), R_BETA(1,isect)
         else
           read(luinp,*) R_H1(isect)
         end if

         if(iidebug .eq. 1) then
            write(6,*)' Sediment values'
            if(NTYPE.eq.1) then
               write(6,*) R_H1(isect), R_R1(1,isect), R_BETA(1,isect)
            else
               write(6,*) R_H1(isect)
            end if
         end if

         if( R_H1(isect) .gt. 0.0 ) then
            do j= 1, nlmax
               if(NTYPE.eq.1) then
                  read(luinp,*,end=510,err=510) 
     >                         R_Z1(j,isect),R_C1(j,isect)
               elseif(NTYPE.eq.2) then
                  read(luinp,*,end=510,err=510) 
     >                         R_Z1(j,isect),R_C1(j,isect),
     >                         R_BLUG1(j,isect),R_R1(j,isect),
     >                         R_BETA(j,isect),R_BLUG2(j,isect)
               elseif(NTYPE.eq.10) then
                  read(luinp,*,end=510,err=510) 
     >                         R_Z1(j,isect),R_C1(j,isect),
     >                         R_R1(j,isect),R_BETA(j,isect)
               end if

               if(iidebug .eq. 1) then
                  if(NTYPE.eq.1) then
                     write(6,*) R_Z1(j,isect), R_C1(j,isect)
                  elseif(NTYPE.eq.2) then
                     write(6,*) R_Z1(j,isect), R_C1(j,isect),
     >                    R_BLUG1(j,isect),R_BLUG2(j,isect)
                  elseif(NTYPE.eq.10) then
                     write(6,*) R_Z1(j,isect), R_C1(j,isect),
     >                    R_R1(j,isect),R_BETA(j,isect)
                  end if
                  if( r_z1(j,isect) .eq. r_h1(isect) )   go to 1400
               end if
            enddo
            print *,' sub opt_read, e^ successo un ...... '
            stop
 1400       continue
            nlayb=j
         else
            nlayb=0
         end if
         R_ND1(ISECT)=nlayb


c: Lower halfspace:
         read(luinp,*)  R_R2(isect), R_BETA(nlayb,isect), R_C2(isect)
         read(luinp,*)  R_BETA(nlayb+1,isect),  R_C2S(isect)
c
         if(iidebug .eq. 1) then
            write(6,*)'Sub-bottom values ' 
            write(6,*)R_R2(isect), R_BETA(nlayb,isect), R_C2(isect)
            write(6,*)  R_BETA(nlayb+1,isect),  R_C2S(isect)
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
      if(iierr .ne. 0) stop 'Stop in svp_read_gama'
c
      return
 510  print *,'Endo or error reading input file at line ',nline
      stop
      end
