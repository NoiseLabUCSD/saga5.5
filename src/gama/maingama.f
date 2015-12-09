      program maingama
c
      implicit integer*4(i-n)
      include './common/Parms_com'
c From saga_getdat1.f
      include 'common/gamaoptions'
      include 'common/srlox'
      include 'common/pii'
      include 'common/freqcom'
      include 'common/caustix'
      include 'common/svp'
      include 'common/plotcom'
c Diff between saga_getdat1 and saga_getdat2.f
      include 'common/depth'
      include 'common/bdcom'
      include 'common/traxcom'
c Diff between saga_getdat3 and saga_getdat1.f+saga_getdat2.f
      include 'common/rthcom' 
      include 'common/pcodes' 
      include 'common/pcode2' 
      include 'common/paths'
      include 'common/bottom' 
      include 'common/scatcom'
      include 'common/charcom' 
c Diff between saga_getsvp and saga_getdat1.f+saga_getdat2.f+saga_getdat3.f
      include 'common/attcon' 
      include 'common/cbotcom'
c Diff between gama and saga_getdat1.f+saga_getdat2.f+saga_getdat3.f
c +saga_getsvp.f
      include './common/timing' 
      include './common/scatarr'
      include './common/causbad'
      include './common/laydata'

      external time
      integer*4 time
      real*4 c60
      data c60/60./
cmay  call usage('gama                ')
      ctab=char(9)
      isec1=time()
      call idate(kdat(1),kdat(2),kdat(3))
      call itime(ktim)
c
      cptim=dtime(cpsec)
c: get input data (opt and svp files).
      call saga_getdat1(nline1)
      call saga_getsvp(nline2)
      call saga_getdat2(nline1)
c:2/1/95
      iibad=0
      call mem_lim(nffth1,NFFTH1M,eline,27,'NFFTH1M',7,iibad)
      call mem_lim(nrtot*nfr2,NRNFMAX,eline,27,'NRNFMAX',7,iibad)
      call mem_lim(nfr2,NFR2MAX,eline,27,'NFR2MAX',7,iibad)
      call mem_lim(nrangx,NRANGXM,eline,27,'NRANGXM',7,iibad)
      call mem_lim(nrtot,NRTOTM,eline,27,'NRTOTM',6,iibad)
      call mem_lim(ntot,NLAYM,eline,27,'NLAYM',5,iibad)

      call saga_getdat3(nline1)
      nlayer=nlay
      call mem_lim(nasl,NASLM,eline,27,'NASLM',5,iibad)
      call mem_lim(nfsl,NFSLM,eline,27,'NFSLM',5,iibad)
      call mem_lim(nasl*nfsl,NAFSLM,eline,27,'NAFSLM',6,iibad)
      call mem_lim(nrth,NRTHMAX,eline,27,'NRTHMAX',7,iibad)
      if(iibad .eq. 1) stop

      do ii=1,5
        call gama
      end do

      cptim=dtime(cpsec)
      write(6,*)
      write(6,*)'SVP FILE: ',lsvpn,svpn(1:lsvpn)
      write(6,*)'OPT FILE: ',loptn,optn(1:loptn)
      write(6,*)'OUT FILE: ',loutn,outn(1:loutn)
      write(6,*)
      write(8,320) svpn(1:lsvpn),optn(1:loptn),outn(1:loutn),
     .   int(cptim/c60),int(mod(cptim,c60)),nrayt
320   format(/'FINISHED GAMARAY RUN FOR '/'SVP FILE = ',a/       
     .   'OPT FILE = ',a/'OUTPUT FILE ROOT = ',a
     .   /'CP MINUTES TAKEN = ',i4,':',i2.2,';  TOTAL ',
     .   'NUMBER OF RAYS PROCESSED = ',i6)

      isec2=time() - isec1
      cptim=etime(cpsec)
      print *,'isec1,isec2 = ',isec1,isec2,time()
      write(8,330) int(cptim/c60),int(mod(cptim,c60)),
     .   int(isec2/60),int(mod(isec2,60)),nraytot
330   format(//'TOTAL CP MINUTES TAKEN FOR RUN = ',i4,':',i2.2,
     .   '; WALL CLOCK TIME = ',i4,':',i2.2,/
     .   'TOTAL NUMBER OF RAYS PROCESSED = ',i6)
      close(8)
      stop
      end
