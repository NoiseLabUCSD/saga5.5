       subroutine write_gama_opt
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
      include './common/arraytype'
c
      integer i,j,jj,ndum,ierr,irflag,lluout
      real rdstep,rderr
      character*80 dumch
      lluout=99
c
 1    format('$$$ OPTION FILE FOR GAMARAY PROGRAM ON ALLIANT $$$',/,
     & '$$$ NOTE: program reads lines following starred lines.',  
     &     ' do not delete',/,
     &     '          any starred lines. do not omit data,',
     &     ' even if not needed. $$$')
 2    format('-----------------------------------------------------',
     &     '------------------')
 3    format('USER OPTIONS: frequency (for eig list, bd); ',
     &     'ray cutoff margin in db;',/,
     &     'ocean svp segments (1=linear, 2=curved); ocean ',
     &     'sv tolerance;',/,
     &     'range interpolation allowed',/,
     &     '*(1) freq  cutdb  iiseg  svtol  rterp)')
 4    format('OPTIONS (1=yes, 0=no):  ray picture; ',
     &     'bar graph; eigenray list;',/,
     &     'include beam displacement; include caustic corrections;',
     &     'output',/,
     &     'data file; output diagnostic messages; metric(1) ',
     &     'or English(0)',/,
     &     '*(2) pic; bar; eig; caus; bd; dat; diag; metric')
 5    format('OPTIONS: iipl(pl vs r); iitf(trans fun vs f); '
     &     'iibmp(beam patfile;',/,
     &     'iiblt(bot loss table); iipl,iitf: 0=none; 1=coh; ',
     &     '2=incoh; 3=both; 4=list',/,
     &     '*(3) iipl; iitf; iibmp;  iiblt  nang  ang1  ang2')
 6    format('FREQUENCIES: nfr>0 ==> give f1,f2,f3...; nfr<0 ==> ',
     &     'give f1,delf',/,
     &     '*(4) nfr  f1  f2(delf)  ...')
 7    format('FFT FILE (iifft), IMPULSE RESPONSE FILE (iiir) options',/,
     &     '[fs=sample f; fmin,fmax=min,max freq of interest in ',
     &     'source spec]',/,
     &     '*(5) iifft  iiir   nfft      fs    fmin   fmax')
 8    format('SOURCE DEPTHS: nzs>0 ==> give zs1,zs2,zs3...; nzs<0 ==> ',
     &     'give zs1,delzs',/,
     &     '*(6) nzs  zs1  zs2(delzs) ...  (m or ft)')
 9    format('SRC-REC RANGES: nr>0 ==> give r1,r2,r3... in m; nr<0 ==> ',
     &     'give r1,delr;',/,
     &     'nr=0 ==> read source track info on next 2 lines',/,
     &     '*(7) nr  r1  r2(dr) ...        (km or nmi) vtilt htilt')
 10   format('*(8) SOURCE TRACK: NLEG(# linear legs); IIPIC ',
     &     '(option to plot track)')
 11   format('LEG SPECIFICATION: V in m/s or kt; NPT>1; T,DT in min;',/,
     &     'CPA,X,Y in km or nmi; PHI in deg',/,
     &     '(for TYPE=1, V=0 ==> T is dist in km; for TYPE=2, ',
     &     'V=0 ==> T2 ignored)',/,
     &     '1 IICONT   V   T1   T2  +DT/-NPT    CPA,PHI    ',
     &     'DUM,DUM (POLAR FORM)',/,
     &     '*(9) 2 IICONT   V   T1  (T2) +DT/-NPT     X1,Y1      ',
     &     'X2,Y2    (X-Y FORM)')
 12   format('*(10) REC ARRAY SPEC: IIARR(1=uniform,2=nonuniform) ',
     &     '(dist in m or ft)')
 13   format('IIARR=1: uniformly spaced receiver array',/,
     &     '*(11) zr1  nx  ny  nz   dx    dy    dz  [dx,dy,dz are ',
     &     'element spacings]')
 14   format('IIARR=2: receiver array depths. nzr>0 ==> give zr1,',
     &     'zr2,... ;',/,
     &     'nzr<0 ==> give zr1,dz2,dz3... (dz are element spacings, ',
     &     'pos down)',/,
     &     '*(12) nzr   zr1  zr2(dz2) zr3(dz3) ...')
 15   format('IIARR=2: specify (x,y) receiver positions for each of ',
     &     'the nzr depths:',/,
     &     '*(13) nxy  x(1),y(1)  x(2),y(2) ... x(nxy),y(nxy) ',
     &     '[nzr lines like this]')
 16   format('OCEAN PATH SPECIFICATION: # standard sets; # ',
     &     'individual paths',/,
     &     '*(14) nstand   nindiv')
 17   format('standard sets: range of bottom interactions; ray ',
     &     'type: "r"=refracting,',/,
     &     '"w"=waterborne, "p"=penetrating, "a"=all ',
     &     '[ex: 0  2  "a"]',/,
     &     '*(15) kb1  kb2  "type"        [nstand lines like this]')
 18   format('individual paths: initial dir ("u" or "d"), #t, #b, ',
     &     'ray type',/,
     &     '*(16) "ini dir"  #t  #b  "type"  [nindiv lines like ',
     &     'this]')
 19   format('BOTTOM PATHS ALLOWED: give maximum # of down-up ',
     &     'traversals T in',/,
     &     'the NLAY bottom layers on EACH bottom interaction ',
     &     '(T usually <= 2)',/,
     &     '*(17) T(1)  T(2)  T(3) ... T(NLAY) [NLAY integer ',
     &     'numbers]')
 20   format('SCATTERING MODEL FOR OCEAN SURFACE, BOTTOM ',
     &     'SUBSTRATE: ii, parm',/,
     &     'none:               ii=0   parm=0',/,
     &     'eckart:             ii=1   parm=wind speed in knots',/,
     &     'modified eckart:    ii=2   parm=wind speed in knots',/,
     &     'beckmann-spizz:     ii=3   parm=wind speed in knots',/,
     &     'uniform:            ii=4   parm=loss in dB per bounce',/,
     &     'loss table:         ii=5   parm="file name" in /inp/',
     &     'gama directory',/,
     &     '*(18) iisl  slparm    FOR OCEAN SURFACE SCATTERING LOSS')
 21   format('*(19) iibl  blparm    FOR BOTTOM SUBSTRATE SCATTERING ',
     &     'LOSS')
 22   format('R-THETA PLOT: iirth (0=no, 1=src angle, 2=rec angle,',/,
     &     '4=MATLAB file of a, R, t [no parameters needed]);',/,
     &     'th1,th2=start, end grazing angles (deg); r1,r2=range '
     &     'interval',/,
     &     '*(20) iirth    th1  th2      r1       r2')

 100  format(2F10.4,I4,2F10.4)
 101  format(8I4)
 102  format(5I4,2F10.4)
 103  format(I4,300F8.2)
 104  format(3I5,3F10.2)
 105  format(I4,4F10.6)
 106  format(2I4,8F8.2)
 107  format(F8.2,3I4,3F8.4)

      write(lluout,1)
      write(lluout,2)
      write(lluout,3)
      write(lluout,100)freq,cutdb,kseg,svtol,rterp
      write(lluout,2)      
      write(lluout,4)
      write(lluout,101)iipic,iibar,iieig,iicaus,iibd,iidat,iidiag,iimet
      write(lluout,2)
      write(lluout,5)
      write(lluout,102)iipl,iitf,iibmp,iiblt,nang,ang1,ang2
      write(lluout,2)
      write(lluout,6)
      write(lluout,103)nfr,(frqsaga(i),i=1,nfr)
      write(lluout,2)
      write(lluout,7)
      write(lluout,104)iifft,iiir,nfft,fs,flo,fhi
      write(lluout,2)
      write(lluout,8)
      write(lluout,103)nzs,(srloc(i),i=1,nzs)
      write(lluout,2)
      write(lluout,9)
      write(lluout,105)-nrangx,rtmp(1),drtmp,dtiltv,dtilth
      write(lluout,10)
      write(lluout,'(2I4)')nleg,iitrx
      write(lluout,11)
      do i=1,nleg
         write(lluout,106)iitype(nl),iicont(nl),vs(nl),
     &         t1(nl),t2(nl),dt(nl),cpa(nl),phid(nl),x2(nl),y2(nl)
      end do
      write(lluout,2)
      write(lluout,12)
      write(lluout,'(I4)')iiarr
      write(lluout,13)
      write(lluout,107)zr1,nx1,ny,nzr,dx,dy,dz(1)      
      write(lluout,14)
      write(lluout,103)nzr,(dz(i),i=1,abs(nzr))
      write(lluout,15)
      irec=0
      do i=1,nzr
         write(lluout,103)nxy(i),(xyz(j,1),
     &         xyz(j,2),j=irec+1,irec+nxy(i))
         irec=irec + nxy(i)
      end do
      write(lluout,2)
      write(lluout,16)
      write(lluout,'(2I4)')nst,ncp0
      write(lluout,17)
      do i=1,nst
         write(lluout,'(2I4,2X,a1)')kb1(i),kb2(i),cz(i)
      end do
      write(lluout,18)
      do i=1,ncp0
         write(lluout,'(2Xa,2I4,2X,a)')fd(i),ktop(i),kbot(i),cz2(i)
      end do
      write(lluout,2)
      write(lluout,19)
      write(lluout,'(100I4)')(maxtrav(i)/2,i=1,ntot)
      write(lluout,2)
      write(lluout,20)
      write(lluout,'(I4,F8.2)')iisl,swind
      write(lluout,21)
      write(lluout,'(I4,F8.2)')iibl,bwind
      write(lluout,2)
      write(lluout,22)
      write(lluout,'(I4,4F10.2)')iirth,thmn,thmx,rmnn,rmxx
c
      close(lluout)
      
      return
      end


