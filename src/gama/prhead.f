      subroutine prhead(nf)
c
      implicit integer*4(i-n)
      include 'common/caustix'
      include 'common/timing'
      include 'common/scatcom'
      include 'common/gamaoptions'
      include 'common/depth'
      include 'common/charcom'
      include 'common/svp'
      include 'common/paths'
      include 'common/srlox'
      include 'common/freqcom'
      include 'common/bottom'
      character*3 labsc(0:5)
      data labsc/'NO ','ECK','M-E','B-S','UNI','TAB'/
c
      write(nf,'(a,a)') 'SVP FILE = ',svpn(1:lsvpn)
      write(nf,'(a,a)') 'OPT FILE = ',optn(1:loptn)
      write(nf,'(a,a)') 'OUTPUT FILE ROOT = ',outn(1:loutn)
      write(nf,100) ktitle(1:lktit),kdat(2),kdat(1),kdat(3),
     .   (ktim(jcan),jcan=1,3)
100   format('SVP TITLE = ',a/'DATE = ',2(i2.2,'/'),i4,
     .   ';  TIME = ',i2.2,2(':',i2.2))
      write(nf,110) cutdb,min0(2,kseg),svtol,int(rterp+.499),
     .   iicaus,iibd,iimet
110   format('EIG CUTOFF MARGIN = ',f7.2,' dB; IISEG = ',i1,         
     .   '; SVTOL = ',f6.3/'RTERP = ',i5,'; IICAUS = ',i1,
     .   '; IIBD = ',i1,'; METRIC = ',i1)
      if(iisl .ne. 5) then
         write(nf,115) labsc(iisl),swind
115      format('SURFACE SCAT = ',a3,'; WIND = ',f5.1,' KT')
      else
         write(nf,116) labsc(iisl),slfile(1:lslf)
116      format('SURFACE SCAT = ',a3,'; SURFACE LOSS FILE = ',a)
      endif
      if(iibl .ne. 5) then
         write(nf,117) labsc(iibl),bwind
117      format('SUBSTRATE SCAT = ',a3,'; WIND = ',f5.1,' KT')
      else
         write(nf,118) labsc(iibl),blfile(1:lblf)
118      format('SUBSTRATE SCAT = ',a3,'; SUBSTRATE LOSS FILE = ',a)
      endif
      write(nf,120) ncp,(maxtrav(jcan)/2,jcan=1,ntot)
120   format('# OCEAN PATHS = ',i3, 
     .   '; MAX # TRAV IN BOTTOM LAYERS = ',42(i2,1x))
      write(nf,130) nzs,nzr,nrec,nrangx,zlev(4),zlev(4+nlev)-zlev(4)
130   format('NZS = ',i2,'; NZR = ',i2,'; # REC = ',i2,
     .   '; # SRC SEQ = ',i4/'WATER DEPTH = ',f9.2,
     .   ' M; SED THICK = ',f8.2,' M')
      write(nf,135) csvp(0),csvp(nsvp),cp1(1),cp2(ntot),cp1(ntot+1)
135   format('SOUND SPEEDS (M/S): OCEAN SURF = ',f7.2,'; OCEAN BOT = ',
     .  f7.2/'                    SED SURF = ',f7.2,'; SED BOT = ',
     .   f7.2,'; SUBSTRATE = ',f7.2)
      write(nf,140) iifft,iiir,nfft,fs
140   format('IIFFT = ',i1,'; IIIR = ',i1,'; NFFT = ',i6,
     .   '; FS = ',f9.2,' HZ')
      write(nf,150) df,flo,fhi,min(twin,999.99)
150   format('DF = ',f8.3,' HZ; FMIN,FMAX = ',f7.0,1X,f7.0,     
     .   ' HZ; TIME WINDOW DUR = ',f8.4,' S')
c
      return
      end
