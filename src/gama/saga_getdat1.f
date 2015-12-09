      subroutine saga_getdat1(nline)
c
c: this subroutine reads the option file up to the number of
c: frequencies desired,nfr.
c
      implicit integer*4(i-n)
      include 'common/gamaoptions'
      include 'common/srlox'
      include 'common/pii'
      include 'common/freqcom'
      include 'common/saga_caustix'
      include 'common/caustix'
      include 'common/svp' 
      include 'common/plotcom'
c

      pie=acos(-1.) 
      pierad=pie/180.
      piedeg=1./pierad
      piequ=pie/4.
      pieh=pie/2.
      twpie=2.*pie
      pietth=2.*pie/3.
      piefth=2.*pietth
      sqpie=sqrt(pie)
      sqpinv=1./sqpie
      sq2pie=sqrt(twpie)
      fax2=-1.*log(10.)/20000.
      fax1=fax2/twpie
      nline=0
      fmin=fminsaga
      fmax=fmaxsaga
cpln      call star(1,nline)
c: read frequency for eig list and bd calcs and eigenray cutoff margin: 
cpln      read(1,*,end=500,err=500) freq,cutdb,kseg,svtol,rterp
      freq=100
      cutdb=60
cpln      cutdb=80
      kseg=1
      svtol=0
cpln Range interpolation (EW suggests 200-400 m)
      rterp=0.
      if(kseg .eq. 2) then
         kseg=3
      elseif(kseg .ne. 1) then
         print *,'illegal ocean svp segments; must be 1 (linear) or ',
     .      '2 (curved):',kseg
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      if(svtol .lt. 0.) then
         print *,'illegal sound vel error tolerance: ',svtol
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      if(cutdb .lt. 0.) then
         print *,'illegal cutdb (should be pos # of db): ',cutdb
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      if(cutdb .lt. 10.) then 
         print *,'warning: eigenrays weaker than cutdb of the ',
     .      'strongest will be ignored.  cutdb = ',cutdb
         print *,'Check line number ',nline,' in .opt file...'
      endif
      if(rterp .lt. 0.) then
         print *,'RTERP is the range interpolation allowed when ',   
     .      'finding eigenrays.  Must be > 0.'
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      cutmag=10**(-cutdb/20.)
c
      if(freq .le. 0.) then
         print *,'illegal frequency: ',freq
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      womega=twpie*freq
c
cpln      call star(1,nline)
cpln      read(1,*,end=500,err=500) iipic,iibar,iieig,iicaus,iibd
      iipic=0
      iibar=0
      iieig=0
      iicaus=0
      iibd=0
c: EKW (3-26-92) iipic neg means use that many points per ray to plot:
      if(iipic .lt. 0) then
         nptmax=min0(200,iabs(iipic))
         print *,'Negative iipic: will plot ',nptmax,' points per ray',
     .      ' (Max is 200).'
         iipic=1
      elseif(iipic .gt. 0) then
         nptmax=120
      endif
      if((iipic .lt. 0) .or. (iipic .gt. 1) .or. (iibar .lt. 0) .or.
     .   (iibar .gt. 1) .or. (iieig .lt. 0) .or. (iieig .gt. 1) .or.
     .   (iibd .gt. 1) .or. (iibd .lt. 0) .or. (iicaus .lt. 0) .or.
     .   (iicaus .gt. 2)) then
         print *,'options in opt file must be either 1 (yes) ',
     .      'or 0 (no): ',iipic,iibar,iieig,iibd,iicaus
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
cpln      backspace(1)
cpln      read(1,*,end=480,err=480) iipicx,iibar,iieig,iicaus,iibd,
cpln     .   iidat,iidiag,iimet
      iipicx=0
      iibar=0
      iieig=0
      iicaus=0
      iibd=0
      iidat=0
      iidiag=0
      iimet=1
      if(iimet .ne. 0 .and. iimet .ne. 1) then
         print *,'metric option must be 1 (metric) or 0 (english)',iimet
         stop
      endif
      if(iidiag .ne. 0) print *,'iidiag not zero: will print ',
     .   'diagnostic messages ...'
      if((iidat .lt. 0) .or. (iidat .gt. 1)) then
         print *,'illegal eigenray data file option: ',iidat
         stop
      endif
      goto 482
480   iidiag=0
      iidat=0
      iimet=1
cpln      backspace(1)
482   continue
      if(iimet .eq. 0) rterp=.30480*rterp
cpln      call star(1,nline)
cpln      read(1,*,err=500,end=500) iipl,iitf,iibmp,iiblt,nang,ang1,ang2
      iipl=4
      iitf=0
      iibmp=0
      iiblt=0
      nang=0
      ang1=0
      ang2=90
      if((iipl .lt. 0) .or. (iipl .gt. 4) .or. (iitf .lt. 0) .or.
     .      (iitf .gt. 4) .or. (iiblt .lt. 0) .or. (iiblt .gt. 1)
     .      .or. (iibmp .lt. 0) .or. (iibmp .gt. 1)) then
         print *,'illegal iipl,iitf,iibmp, or iiblt option: ',iipl,
     .      iitf,iibmp,iiblt
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      if(iiblt .eq. 1) then
         if((nang .le. 0) .or. (ang1 .lt. 0.) .or. (ang2 .lt. ang1))
     .      then
            print *,'illegal nang,ang1,ang2 for iiblt option: ',nang,
     .         ang1,ang2
            print *,'Check line number ',nline,' in .opt file...'
            stop
         endif
         if(ang1 .lt. .5) then
            print *,'angle 1 for bottom loss changed to .5 degrees ',
     .         'to avoid very long ranges '
            ang1=.5
         endif
      endif
      if(iibmp .eq. 1 .and. iicaus .eq. 1) then
         print *,'caustic corrections option set to 0 since iibmp=1'
         iicaus=0
      endif
      iipltf=iipl + iitf + iiblt + iibmp
cpln      call star(1,nline)
cpln      read(1,*,end=500,err=500) nfr
cpln      nfr=1
cpln      backspace(1) 
      if(iipltf .eq. 0) then
         print *,'nfr in line ',nline,' set to zero since pl,tf,blt,',
     .      'bmp options = 0.'
         nfr=0 
      endif
      if(abs(nfr) .gt. 105) then
         print *,'# freq for pl or tf plots cannot exceed 105. ',nfr
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      if(iipltf .ne. 0) then
         if(nfr .gt. 0) then
cpln            read(1,*,err=500,end=500) nfr,(frq(j),j=1,nfr)
cpln            frq(1)=60
         elseif(nfr .lt. 0) then
cpln            read(1,*,err=500,end=500) nfr,frq(1),dfrq 
cpln            frq(1)=60
cpln            nfr=-1*nfr
cpln            do 110 kfr=2,nfr
cpln               frq(kfr)=frq(kfr-1) + dfrq
cpln 110         continue
         else
            iipl=0
            iitf=0
            iiblt=0
            print *,'iipl,iitf,iiblt set to 0 since nfr=0'
         endif
      endif
c
c: read in input for fft files [alliant only]: 
cpln      call star(1,nline)
cpln      read(1,*,end=500,err=500) iifft,iiir,nfft,fs,flo,fhi
      iifft=2
      iiir=0
cpln      nfft=8192
cpln      fs=6000
cpln      flo=220
cpln      fhi=800
c added pln 241100
      flotmp=flo
      fhitmp=fhi
      if((iifft .gt. 6) .or. (iifft .lt. 0) .or. (iiir .gt. 1) .or.
     .   (iiir .lt. 0)) then
         print *,'illegal iifft,iiir,nfft in opt file: ',iifft,iiir
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      iiffi=iifft + iiir
      if(iiffi .ne. 0) then
         if(fs .le. 0. .or. flo .le. 0. .or. fhi .lt. flo) then
            print *,'illegal fs, flo, or fhi in opt file: ',fs,flo,fhi
            print *,'Check line number ',nline,' in .opt file...'
            stop 'Peter Nielsen'
         endif
         if(nfft .le. 0) then
            print *,'illegal nfft in opt file: ',nfft
            print *,'Check line number ',nline,' in .opt file...'
            stop
         endif
cpln         kp2=int(log10(float(nfft))/log10(2.) + .49)
cpln         nfft=2**kp2
cpln         twintot=.92*float(nfft)/fs
         twintot=float(nfft)/fs
c: set twin, the time window at the receiver, according to the options:
         if((iiir .eq. 1) .or. (mod(iifft,3) .eq. 1)) then
c: for any kind of correlation processing, make sure last half of
c: time series will be zero:
            twin=twintot/2.
c: NEW. IIFFT=3 10-1-90:
         elseif((iiir .eq. 2) .or. (mod(iifft,3) .eq. 0)) then
c: for iifft=3,6 or impulse responses, allow arrivals in full window:
            twin=twintot
         else
c: for only propagation loss fft's, allow all ray arrivals:
            twin=1.e13
         endif
      else
         twin=1.e13
         nfft=4
         fs=4.
      endif
      nffth=nfft/2
      nffth1=nffth + 1
      df=fs/nfft
c     ntap=19*nffth/20
      ntap=9*nffth1/10
      dtap=pie/(nffth1 - ntap) 
      nfft2=nfft*2
      nfft1=nfft + 1
      nfft12=2*nfft1
c
      return
c
cpln 500   call readerr(1,nline)
      return
      end 
