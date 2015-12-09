      subroutine saga_getdat2(nline)
c
c: this subroutine gets the second batch of input data.
c
      implicit integer*4(i-n)
      include 'common/gamaoptions'
      include 'common/svp'
      include 'common/caustix'
      include 'common/pii'
      include 'common/freqcom'
      include 'common/srlox'
      include 'common/depth'
      include 'common/bdcom'
      include 'common/traxcom'
      include 'common/tilt'
      include 'common/arraytype'
c: nbdfmax is the max # of freq at which to calc bd eigs, fsfac is the
c: factor that controls the spacing of the frequecies:
      common /info/ienv_info,lu_env,iiwrite
      integer ienv_info,lu_env,iiwrite
      data nbdfmax/19/,fsfac/2.5/
c
c: iibet is zero when iibar=iieig=0
      iibet=iibar+iieig
c

      fmin=1.e38
      do 10 kfr=1,nfr
         if(frq(kfr) .lt. fmin) then
            fmin=frq(kfr)
            nfmin=kfr
         endif
10    continue
      if(fmin .le. 0.) then
         print *,'illegal freq for prop loss plots: ',fmin
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      nfr2=nfr
      iiadd=0
      if((iibet .ne. 0) .or. (iipic .ne. 0) .or. (nfr .eq. 0)) iiadd=1
      if(iiffi .ne. 0) then
         if(flo .lt. fmin) then
            if((iiadd .eq. 1) .and. (flo .lt. freq)) then
               nfr2=nfr2 + 1
               frq(nfr2)=flo
               nfmin=nfr2
               fmin=flo
            endif
         endif
      endif
      if(iiadd .eq. 1) then
         nfr2=nfr2 + 1
         frq(nfr2)=freq
         if(freq .lt. fmin) then
            nfmin=nfr2
            fmin=freq
         endif
      endif
c: put freq in frq array for calculating prop loss: 
c: find low,high f of interest for caustic calcs and range error tol: 
      if(iiffi .ne. 0) then
c: limit fhi to the Nyquist frequency:
         fhi=min(fhi,fs/2.)
         fmax=fhi
         fmaxx=max(fhi,fs/2.) + 1.
      else
         fmax=0
      endif
      do 11 kfr=1,nfr2
         fmax=max(fmax,frq(kfr))
         wom(kfr)=twpie*frq(kfr)
         w16(kfr)=wom(kfr)**(.16667)
         w23(kfr)=wom(kfr)**(.66667)
         frsq(kfr)=frq(kfr)**2
11    continue
      fmaxx=max(fmaxx,fmax + 1.)
c: wmin used in caustic correction subroutines: 
      wmin=twpie*fmin
      wm16=wmin**(.16667)
      wm23=wmin**(.66667)
c
      bdf(0)=0.
      nfeig=0
      if((fmax .eq. fmin) .or. (iibd .eq. 0)) then
         nbdf=1
         bdf(1)=fmin
         nfeig=1
         goto 360
      elseif((iiffi .eq. 0) .and. (nfr2 .eq. 2)) then
         nbdf=2
      else
c: compute the freq at which the total waveguide depth is 50 lambda:
         f50=max(min(fmax,50.*csvp(nsvp)/zlev(4+nlev)),2.*fmin)
280      nbdf=nint(log10(f50/fmin)/log10(fsfac))+1
         if(nbdf .gt. nbdfmax) then
            fsfac=fsfac + .25
            print *,'changing fsfac in getdat1 ... '
            goto 280
         endif
         do 290 j=2,nbdf-1
            bdf(j)=fmin*fsfac**(j-1)
290      continue
      endif
      bdf(1)=fmin
      bdf(nbdf)=fmax
      if((iibet .ne. 0) .or. (iipic .ne. 0)) then
         do 350 j=1,nbdf
            if((freq .gt. bdf(j-1)) .and. (freq .lt. bdf(j))) then
               do 352 jj=nbdf+1,j+1,-1
                  bdf(jj)=bdf(jj-1)
352            continue
               nfeig=j
               bdf(nfeig)=freq
               nbdf=nbdf + 1
               goto 360
            elseif(freq .eq. bdf(j)) then
               nfeig=j
            endif
350      continue
      endif
360   continue
      bdf(nbdf+1)=2.*max(fmax,fs/2.)
      do 355 j=1,nbdf
         bdom(j)=twpie*bdf(j)
         bdom16(j)=bdom(j)**(.166667)
         bdom23(j)=bdom(j)**(.666667)
355   continue
      if((iibd .eq. 1) .and. (nbdf .gt. 1)) then
         write(8,370)'FREQUENCIES AT WHICH BD EIGENRAYS CALCULATED: ',
     .      (bdf(jcan),jcan=1,nbdf)
370      format(a,20(f7.1,1x))
      endif
c
cpln      call star(1,nline)
cpln      read(1,*,end=500,err=500) nzs
cpln      nzs=1
      if((abs(nzs) .gt. 60) .or. (nzs .eq. 0)) then
         print *,'illegal nzs in opt file: ',nzs
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
cpln      backspace(1)
      if(nzs .gt. 0) then
cpln         read(1,*,end=500,err=500) nzs,(srloc(j),j=1,nzs)
cpln         nzs=1
cpln         srloc(1)=54.331
         if(iimet .eq. 0) then
            do 940 jcan=1,nzs
               srloc(jcan)=.30480*srloc(jcan)
940         continue
         endif
         do 230 j=1,nzs
            if(srloc(j) .lt. 0.) then
               print *,'negative src depth ',srloc(j),
     .           ' measured from ocean bottom: ',zsvp(nsvp) + srloc(j)
               srloc(j)=zsvp(nsvp) + srloc(j)
            endif
230      continue
      else
cpln         read(1,*,end=500,err=500) nzs,srloc(1),dzs 
         nzs=-1
         srloc(1)=-1
         dzs=-1
         nzs=abs(nzs)
         if(iimet .eq. 0) srloc(1)=.30480*srloc(1)
         if(iimet .eq. 0) dzs=.30480*dzs
         if(srloc(1) .lt. 0.) then
            print *,'negative src depth ',srloc(1),
     .        ' measured from ocean bottom: ',zsvp(nsvp) + srloc(1)
            srloc(1)=zsvp(nsvp) + srloc(1)
         endif
         do 30 jzs=2,nzs
            srloc(jzs)=srloc(jzs-1) + dzs
30       continue
      endif
c
cpln      call star(1,nline)
cpln      read(1,*,end=500,err=500) nrangx
cpln      nrangx=-128
c      nrangx=-33
c nrangx is No of source-receiver ranges
      if(nrangx .eq. 0) then
         iirr=0
cpln         call star(1,nline)
cpln         read(1,*,end=500,err=500) nleg,iitrx
         nleg=1
         iitrx=1
         if((nleg .lt. 1) .or. (nleg .gt. 50)) then
            print *,'illegal nleg (1<=nleg<=50): ',nleg
            print *,'Check line number ',nline,' in .opt file...'
            stop
         endif
cpln         call star(1,nline)
c: nrleg is the # of source positions on the current leg;
c: nrangx is the total number of source track positions: 
         nrangx=0
         do 40 nl=1,nleg
cpln            read(1,*,end=500,err=500) iitype(nl),iicont(nl),vs(nl),
cpln     .         t1(nl),t2(nl),dt(nl),cpa(nl),phid(nl),x2(nl),y2(nl)
            if(iimet .eq. 0) then
               vs(nl)=.51480*vs(nl)                  
               cpa(nl)=1.8520*cpa(nl)
               if(iitype(nl) .eq. 2) phid(nl)=1.8520*phid(nl)
               x2(nl)=1.8520*x2(nl)
               y2(nl)=1.8520*y2(nl)
            endif
c: For vs not zero, set t2 according to distance traveled and velocity:
            if(vs(nl) .ne. 0 .and. iitype(nl) .eq. 2) then
               t2(nl)=t1(nl) + 1000.*sqrt((x2(nl)-x1(nl))**2 +
     .            (y2(nl)-y1(nl))**2)/(vs(nl)*60.)
            endif
            if(dt(nl) .lt. 0.) then
               nrleg(nl)=iabs(nint(dt(nl)))
               if(nrleg(nl) - abs(dt(nl)) .ne. 0.) then
                  print *,'NPT on source track must be integer: ',nl,
     .               dt(nl)
                  stop
               endif
            else
               nrleg(nl)=nint((t2(nl) - t1(nl))/dt(nl)) + 1
            endif
            if(iicont(nl) .ne. 0 .and. iicont(nl) .ne. 1) then
               print *,'illegal iicont(0 or 1) in track: ',nl,iicont(nl)
               stop
            endif
            if((iicont(nl) .eq. 0 .and. nrleg(nl) .lt. 1) .or.
     .         (iicont(nl) .eq. 1 .and. nrleg(nl) .lt. 2)) then
               print *,'NPT for iicont=0 leg must be > 0, ',
     .   'NPT for iicont=1 leg must be > 1: ',nl,dt(nl),nrleg(nl)
               stop
            endif
            nrangx=nrangx + nrleg(nl)
            if(iicont(nl) .eq. 1) nrangx=nrangx - 1
40       continue
         iicont(1)=0
      else
cpln         backspace(1)
         if(nrangx .gt. 0) then
            iirr=1
            if(nrangx .gt. 101) then
               print *,'# ranges limited to 101 when you list them.'
               print *,'Check line number ',nline,' in .opt file...'
               stop
            endif
cpln            read(1,*,end=500,err=500) nrangx,(rtmp(j),j=1,nrangx)
cpln            nrangx=128
cpln            rtmp(1)=-1
         else
            iirr=-1
cpln            read(1,*,end=500,err=500) nrangx,rtmp(1),drtmp
cpln            nrangx=-128
cplc            nrangx=-33
cpln            rtmp(1)=.306457
cpln            drtmp=0.002
            nrangx=iabs(nrangx)
         endif
cpln         call star(1,nline)
cpln         call star(1,nline)
      endif
c
c: read in receiver array positions: 
cpln      call star(1,nline)
cpln      read(1,*,end=500,err=500) iiarr
      iiarr=1
cpln      call star(1,nline)
c: read in uniform array locations: 
      if(iiarr .eq. 1) then
cpln         read(1,*,end=500,err=500) zr1,nx,ny,nzr,dx,dy,dz(1)
cpln         zr1=56.299
         if(tilth) then
c            write(6,*)'dtilth: ',dtilth
c            write(6,*)'nrangx, nzr: ',nrangx,nzr
c            write(6,*)'drtmp0: ',drtmp0
c            nzr=nrangx
            dz(1)=dtilth/float(nzr-1)
            drtmp=dsqrt((drtmp0)**2
     >            -(0.001*abs(dtilth)/dble(nzr-1))**2)
c            drtmp=drtmp0
            if(iiwrite.gt.0) then
               write(6,*)
               write(6,*)'****************************'
               write(6,*)'  Horizontal tilt included  '
               write(6,*)'  DR, NRANGX: ',drtmp,nrangx
               write(6,*)'  DZ, NZRH  : ',dz(1),nzr
               write(6,*)'****************************'
               write(6,*)
            end if
         end if
c
         if(tiltv) then
c            nrangx=nzr
c            nrangx=-1
            nrangx=1
            drtmp=dtiltv/(nzr-1)*0.001
            dz(1)=dsqrt((dz0(1))**2
     >           -(abs(dtiltv)/dble(nzr-1))**2)
            if(iiwrite.gt.0) then
               write(6,*)
               write(6,*)'****************************'
               write(6,*)'  Vertical tilt included    '
               write(6,*)'  DR, NRANGX: ',drtmp,nrangx
               write(6,*)'  DZ, NZRH  : ',dz(1),nzr
               write(6,*)'****************************'
               write(6,*)
            end if
         end if

         nx1=1
         ny=1
cpln         nzr=1
         dx=0.
         dy=0
cpln         dz(1)=0
         if(iimet .eq. 0) then
            zr1=.30480*zr1
            dx=.30480*dx
            dy=.30480*dy
            dz(1)=.30480*dz(1)
         endif
         if(zr1 .lt. 0.) then
            print *,'negative rec depth ',zr1,' measured from ocean ',
     .         'bottom: ',zsvp(nsvp) + zr1
            zr1=zsvp(nsvp) + zr1
         endif
         nxny=nx1*ny 
         if((nx1 .lt. 1) .or. (ny .lt. 1) .or. (nzr .lt. 1)
     .      .or. (nzr .gt. 128) .or. (nxny .gt. 32)) then
            print *,'illegal nx,ny,nz for iiarr=1 in opt file: ',
     .         nx1,ny,nzr
            print *,'Check line number ',nline,' in .opt file...'
            stop
         endif
         dxtot=float(nx1-1)*dx 
         dytot=float(ny-1)*dy 
         nrec=nxny*nzr
         jxy=0
         do 20 jzr=1,nzr
c: zr1 is depth of top of array; recs are numbered beginning from top:
            zr0=zr1 + float(jzr-1)*dz(1)
            nxy(jzr)=nxny
            kxy(jzr)=jxy
c: fill xyz(jxy,i) with (x,y) coordinates from center of array: 
            do 22 jy=1,ny
               do 24 jx=1,nx1
                  jxy=jxy + 1 
                  xyz(jxy,1)=float(jx-1)*dx - dxtot/2.
                  xyz(jxy,2)=float(jy-1)*dy - dytot/2.
                  xyz(jxy,3)=zr0
24             continue
22          continue
20       continue
cpln         call star(1,nline)
cpln         call star(1,nline)
c: read in nonuniform array locations:  
      elseif(iiarr .eq. 2) then
cpln         call star(1,nline)
c: CHANGED FORMAT OF IIARR SPECIFICATION 4-10-90: 
c        read(1,*,end=500,err=500) zr1,nzr,(dz(j),j=1,nzr)
cpln         read(1,*,end=500,err=500) nzr
         if((abs(nzr) .gt. 60) .or. (nzr .eq. 0)) then
            print *,'illegal nzr for iiarr=2 in opt file: ',nzr
            print *,'Check line number ',nline,' in .opt file...'
            stop
         endif
cpln         backspace(1)
         if(nzr .gt. 0) then
cpln            read(1,*,end=500,err=500) nzr,(dz(j),j=1,nzr)
            if(iimet .eq. 0) then
               do 942 jcan=1,nzr
                  dz(jcan)=.30480*dz(jcan)
942            continue
            endif
            do 232 j=1,nzr
               if(dz(j) .lt. 0.) then
                  print *,'negative rec depth ',dz(j),
     .  ' measured from ocean bottom: ',zsvp(nsvp) + dz(j)
                  dz(j)=zsvp(nsvp) + dz(j)
               endif
232         continue
         else
cpln            read(1,*,end=500,err=500) nzr,(dz(j),j=1,abs(nzr))
            nzr=abs(nzr)
            if(iimet .eq. 0) then
               do 944 jcan=1,nzr
                  dz(jcan)=.30480*dz(jcan)
944            continue
            endif
            if(dz(1) .lt. 0.) then
               print *,'negative rec depth ',dz(1),
     .            ' measured from ocean bottom: ',zsvp(nsvp) + dz(1)
               dz(1)=zsvp(nsvp) + dz(1)
            endif
            do 32 jzr=2,nzr
               dz(jzr)=dz(jzr-1) + dz(jzr)
32          continue
         endif
c: nrec will be the total number of receivers in the array: 
         nrec=0
cpln         call star(1,nline)
c: read nzr lines of: # of receivers and their (x,y) locations: 
         do 14 jzr=1,nzr
            kxy(jzr)=nrec
cpln            read(1,*,end=500,err=500) nxy(jzr),(xyz(j,1),
cpln     .         xyz(j,2),j=nrec+1,nrec+nxy(jzr))
            do 946 jcan=nrec+1,nrec+nxy(jzr)
               xyz(jcan,3)=dz(jzr)
946         continue
            if(iimet .eq. 0) then
               do 948 jcan=nrec+1,nrec+nxy(jzr)
                  xyz(jcan,1)=.30480*xyz(jcan,1)
                  xyz(jcan,2)=.30480*xyz(jcan,2)
948            continue
            endif
            nrec=nrec + nxy(jzr)
            if((nrec .gt. 128) .or. (nxy(jzr) .lt. 1)) then
               print *,'illegal nxy ',nxy(jzr),' for jzr = ',jzr
               print *,'# rec cannot exceed 128: ',nrec
               print *,'Check line number ',nline,' in .opt file...'
               stop 
            endif
14       continue
      else
         print *,'illegal rec array type iiarr in opt file: ',iiarr
         print *,'Check line number ',nline-1,' in .opt file...'
         stop
      endif
c
      if(iiblt .eq. 1) then
         nzs=1
         nzr=1
         nrec=1
         nrangx=nang
         print *,'NOTE: src/rec ranges being set for bottom loss calcs'
         rterp=0.
      endif
c: nrtot is the total number of ranges in range() for each zs.
      nrtot=nrec*nrangx
      nrnf=max0(nrtot,nfr2)

      return
      end 
