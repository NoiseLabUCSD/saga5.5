       subroutine write_info
c     reads and interpretates the input file to the genetic algorithm 
c     optimization program
c     PETER GERSTOFT, 1992
c
      implicit integer*4(i-n)
      INCLUDE '../comopt.h'
      INCLUDE '../comforw.h'
      include './common/Parms_com'
      include './common/parms_conv'
c From saga_getdat1.f
      include './common/gamaoptions'
      include './common/srlox'
      include './common/pii'
      include './common/saga_freqcom'
      include './common/saga_caustix'
      include './common/svp'
      include './common/plotcom'
c Diff between saga_getdat1 and saga_getdat2.f
      include './common/depth'
      include './common/bdcom'
      include './common/traxcom'
c Diff between saga_getdat3 and saga_getdat1.f+saga_getdat2.f
      include './common/rthcom' 
      include './common/pcodes' 
      include './common/pcode2' 
      include './common/paths'
      include './common/bottom' 
      include './common/scatcom'
      include './common/charcom' 
c Diff between saga_getsvp and saga_getdat1.f+saga_getdat2.f+saga_getdat3.f
      include './common/attcon' 
      include './common/cbotcom'
c Diff between gama and saga_getdat1.f+saga_getdat2.f+saga_getdat3.f
c +saga_getsvp.f
      include './common/timing' 
      include './common/scatarr'
      include './common/causbad'
      include './common/laydata'
      include './common/tilt'
c
      integer i,j,jj,ndum,ierr,irflag,lluout
      real rdstep,rderr
      character*80 dumch
      lluout=55
c
      write(lluout,*)
      write(lluout,*)'No of sectors: ',NSCTOR
      write(lluout,*)
      write(lluout,*)'Channel information: Water column'
      write(lluout,*)cp1(-1),cs1(-1),rho1(-1),akp1(-1),aks1(-1)
      write(lluout,*)'No of layers: ',nsvp
      write(lluout,*)zsvp(0),csvp(0),rho2(0),akp2(0)
      do i=1,nsvp
         write(lluout,*)zsvp(i),csvp(i)
      end do

      write(lluout,*)
      write(lluout,*)'Channel information: Sediment'
      write(lluout,*)'Bottom type: ',ntype
      write(lluout,*)'No of layers: ',ntot
      do i=1,ntot
         write(lluout,*)z(i),cp1(i),cs1(i),rho1(i),akp1(i),aks1(i)
         write(lluout,*)z(i),cp2(i),cs2(i),rho2(i),akp2(i),aks2(i)
         write(lluout,*)
      end do
      write(lluout,*)
      write(lluout,*)'Channel information: Sub-bottom'
      write(lluout,*)z(ntot+1),cp1(ntot+1),cs1(ntot+1),
     .           rho1(ntot+1),akp1(ntot+1),aks1(ntot+1)
c
      write(lluout,*)
      write(lluout,*)'Source-receiver geometry'
      write(lluout,*)nzs,srloc(1),dzs
      write(lluout,*)nrangx,rtmp(1),drtmp
      write(lluout,*)zr1,nx,ny,nzr,dx,dy,dz(1)
      write(lluout,*)tilth,tiltv,dtilth,dtiltv

      write(lluout,*)
      write(lluout,*)'Frequency info'
      write(lluout,*)freq,cutdb,kseg,rterp
      write(lluout,*)nfr
      do i=1,nfr
        write(lluout,*)frqsaga(i)
      end do
       write(lluout,*)nfft,flo,fhi,fs

       close(lluout)
      return
      end


