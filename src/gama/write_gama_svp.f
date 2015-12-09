       subroutine write_gama_svp
c     reads and interpretates the input file to the genetic algorithm 
c     optimization program
c     PETER GERSTOFT, 1992
c
      implicit integer*4(i-n)
c      INCLUDE '../comopt.h'
c      INCLUDE '../comforw.h'
      include './common/Parms_com'
      include './common/parms_conv'
c From saga_getdat1.f
      include './common/gamaoptions'
      include './common/srlox'
      include './common/pii'
      include './common/saga_freqcom'
      include './common/caustix'
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
      lluout=98
c
 1    format('$$$ SVP FILE FOR GAMARAY PROGRAM ON ALLIANT $$$',/,
     & '$$$ NOTE: program reads lines following starred lines.',  
     &     ' do not delete',/,
     &     '          any starred lines. do not omit data,',
     &     ' even if not needed. $$$')
 2    format('-----------------------------------------------------',
     &     '------------------')
 3    format('*(1) TITLE (up to 64 characters):')
 4    format('GAMA INPUT FROM SAGA')
 5    format('GEOACOUSTIC PARAMETERS OF HALFSPACE ABOVE OCEAN:',/,
     &        '*(2) cp; cs; rho; kp; ks;')
 6    format('SOUND VELOCITY PROFILE IN OCEAN:',/,
     &       '*(3) nsvp = # of depth, sound velocity points')
 7    format('*(4) depth, sound velocity points')
 8    format('GEOACOUSTIC PARAMETERS OF OCEAN BOTTOM LAYERS:',/,
     &     '*(5) number of bottom layers')
 9    format('LAYER PROFILES: type: 1=linear, 2=blug, 3=blug,'
     &     '4=blug; examples:',/,
     &     '1  thick  cp1 cp2  cs1 cs2  rho1 rho2  kp1 kp2  ks1 ks2',/,
     &     '2  thick  cp1   g  cs1 cs2  rho1 rho2  kp1 kp2  ks1 ks2'
     &     'beta rblug',/,
     &     '3  thick  cp1 cp2  cs1 cs2  rho1 rho2  kp1 kp2  ks1 ks2'
     &     ' beta rblug',/,
     &     '4  thick  cp1 cp2  cs1 cs2  rho1 rho2  kp1 kp2  ks1 ks2'
     &     ' g rblug',/,
     &     '[beta=-999 ==> set beta such that g continuous from'
     &     ' previous layer,',/,
     &     'rblug=0   ==> use usual plane wave r/t coeff at bottom'
     &     'of layer]',/,'*(6)')
 10   format('GEOACOUSTIC PARAMETERS OF SUBSTRATE:',/,
     &     '*(7) cp; cs; rho; kp; ks; in substrate')
 100  format(5F10.4)
 101  format(2F10.4)
 102  format(I4,F8.4,10F10.4)
 103  format(4F10.4)

      write(lluout,1)
      write(lluout,2)
      write(lluout,3)
      write(lluout,4)
      write(lluout,2)
      write(lluout,5)
      write(lluout,100)cp1(-1),cs1(-1),rho1(-1),akp1(-1),aks1(-1)
      write(lluout,2)
      write(lluout,6)
      write(lluout,'(I4)')nsvp+1
      write(lluout,7)
      write(lluout,103)zsvp(0),csvp(0),rho2(0),akp2(0)
      do i=1,nsvp
         write(lluout,101)zsvp(i),csvp(i)
      end do
      write(lluout,2)
      write(lluout,8)
      write(lluout,'(I4)')ntot
      write(lluout,9)
      do i=1,ntot
         write(lluout,102)kprof(i),z(i),cp1(i),cp2(i),cs1(i),cs2(i),
     &        rho1(i),rho2(i),akp1(i),akp2(i),aks1(i),aks2(i)
      end do
      write(lluout,2)
      write(lluout,10)
      write(lluout,100)cp1(ntot+1),cs1(ntot+1),
     &           rho1(ntot+1),akp1(ntot+1),aks1(ntot+1)
c
      close(lluout)
      
      return
      end


