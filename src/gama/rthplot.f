      subroutine rthplot(iirth,ntlay,NLAYM,facax,rtmat)
c
c: this subroutine plots range vs. theta plots.
c
      implicit integer*4(i-n)
      include 'common/rthcom' 
      include 'common/depth'
      include 'common/pii'
      include 'common/laydata'
      include 'common/svp'
      include 'common/paths'
      include 'common/discpt' 
      include 'common/pathway'
      include 'common/pathchr'
      include 'common/charcom'
      include 'common/bdcom'
      include 'common/weakcom'
      include 'common/caustix'
      character*32 labths,labthr,labr,labb
      character*64 matfile(10)
      character*3 c1
      integer*4 mathead(5,3),npmax,nfunm
      parameter(npmax=200,nfunm=1000)
      real*8 facax(nrth),rtmat(3,nrth),avals(npmax,nfunm),
     .   rvals(npmax,nfunm),tvals(npmax,nfunm),mvals(npmax,nfunm),
     .   raypaths(5,nfunm)

cxx   data labths/'SOURCE GRAZING ANGLE (DEG)'/
cxx   data labthr/'RECEIVER GRAZING ANGLE (DEG)'/
cxx   data labr/'RANGE (M)'/,labb/' '/
c
      do 940 jcan=1,nfunm
         do 942 jcan2=1,npmax
            avals(jcan2,jcan)=0.
            rvals(jcan2,jcan)=0.
            tvals(jcan2,jcan)=0.
            mvals(jcan2,jcan)=0.
942      continue
940   continue
      thmn=max(1.e-6,thmn)
      thmx=min(90-1.e-6,thmx)
      if(iirth .ge. 4) then
         thmn=1.e-3
         thmx=90-thmn
         if(iirth .eq. 4) then
            nval=100
         else
            nval=iirth
         endif
         iirth=4
      endif
      call axlabel(thmn,thmx,thlo,thhi,thdelt)
      call axlabel(rmnn,rmxx,rlo,rhi,rdelt)
      dimh=10.75
      dimv=8.0
      icxx=-1
c     call pltlfn(l"raypic")
c     call pltdim(11.75,9.0,-1)
c     call pltorg(.5,.5)
c     if(iirth .eq. 1) then
c        call pltaxis(0.,0.,dimh,0.,thhi,thlo,thdelt/10.,labths,-27,10)
c     else
c        call pltaxis(0.,0.,dimh,0.,thhi,thlo,thdelt/10.,labthr,-29,10)
c     endif
c     call pltaxis(0.,dimv,dimh,0.,thhi,thlo,thdelt/10.,labb,1,0)
c     call pltaxis(0.,0.,dimv,90.,rlo,rhi,rdelt/10.,labr,9,10)
c     call pltaxis(dimh,0.,dimv,90.,rlo,rhi,rdelt/10.,labb,-1,0)
      hhfac=dimh/(thhi-thlo)
      vvfac=dimv/(rhi-rlo)
      thmin=thmn*pierad
      thmax=thmx*pierad
      thint=thmax-thmin
c
      kbdf=nfeig
      cref=cs
      if(iirth .eq. 2) cref=cr
      if(iirth .eq. 3) then
c: For iirth=3, output is a matlab file of range vs time for the ncp
c: ocean paths and only the first bottom path:
         open(52,file=outn(1:loutn)//'_rt',status='unknown',
     .        form='unformatted')
c set MATLAB header values:
         lmat=12
         mathead(1,1) = 1000       ! MATLAB type flag
         mathead(2,1) = nrth       ! Number of rows in the matrix
         mathead(3,1) = 3          ! Number of columns in the matrix
         mathead(4,1) = 0          ! Imaginary flag (=false)
         mathead(5,1) = lmat       ! Length in bytes of variable name +1
         write(c1(3:3),'(a1)') char(0)
         do 145 jpt=1,nrth
            facax(jpt)=(float(jpt-1)/float(nrth-1))**(.3333)
145      continue
      endif
      nbwk=0
c: loop for the ncp number of ocean paths
      cscr=max(cs,cr)
      nrtfile=0
      nfun=0
c: EKW: 12-12-91
      mbp=0
      npointmax=0
      do 10 icpth=1,ncp
         ntop=ncpth(icpth,1)/2
         nbotx=ncpth(icpth,3)/2
         nttot(0)=2*nbotx
         call setop(kst,kend,kdir,czprop(icpth))
         mmr0b=1
 44      continue
cpln         write(6,*)
cpln         write(6,*)'RTHPLOT NLAYM: ',NLAYM
         call setbp(kst,kstt,kend,kendd,norefl,czprop(icpth),
     .      kdone,mmr0b,ntlay,NLAYM)
cxx   print *,'after setbp: kstt,kendd,norefl,kdone=',kstt,kendd,
cxx  .   norefl,kdone
         if(kdone .eq. 1) goto 10
         a1=aaxis(1)
         a2=aaxis(kdisc(kendd)/10 - 1)
         thint=asin(a2*cref) - asin(a1*cref)
cxx   print *,'icpth = ',icpth,ncp,nrth,nbotx,ntop
         if(iirth .eq. 3) then
            sinth1=sqrt(1. - (a1*cscr)**2)
            sinth2=sqrt(1. - (a2*cscr)**2)
            delsin=sinth2 - sinth1
cxx         dela=(a2 - a1)/max0(1,nrth-1)
cxx   print *,'icpth = ',icpth,ncp,nrth,nbotx,ntop
            do 150 jpt=1,nrth
               sinth=sinth1 + facax(jpt)*delsin
cxx            apt=a1 + (jpt-1)*dela
               apt=max(1.e-10,sqrt(1. - sinth**2)/cscr)
               call rdrcalc(apt,rtmat(1,jpt),dr,ibd,rbd,drbd)
               call timatt(apt,attpt,rtmat(2,jpt),ibd)
               rtmat(3,jpt)=apt
150         continue
            nrtfile=nrtfile + 1
            c1(1:1)=char(48 + nrtfile/10)
            c1(2:2)=char(48 + mod(nrtfile,10))
            matfile(1)='rng_time_'//c1(1:3)
            write(52) (mathead(jcan,1),jcan=1,5),matfile(1)(1:lmat)
            write(52) ((rtmat(jcan2,jcan),jcan=1,nrth),jcan2=1,3)
            goto 44
         endif
         thpl=0.
         rpl=-.20
         kaxlast=kdisc(kstt)/10
         thlast=thmax
         thnext=max(acos(aaxis(kaxlast-1)*cref),thmin)
         if(nfun .ge. nfunm) then
            print *,'# ray paths for iirth=4 > nfunm, skipping: ',nfunm
            goto 10
         endif
         nfun=nfun + 1
         if(nfun .gt. nfunm) then
            print *,'Increase array sizes in rthplot!! nfun = ',nfun
            stop
         endif
         raysign=1.
         if(ldir(icpth) .eq. 'u') raysign=-1.
         raypaths(1,nfun)=1.
         if(fdir(icpth) .eq. 'd') raypaths(1,nfun)=-1.
         raypaths(2,nfun)=raysign
         raypaths(3,nfun)=ntop
         raypaths(4,nfun)=nbotx
         raypaths(5,nfun)=nttot(1)
cxx   print *,'nfun = ',nfun
         npoint=0
      if(iidiag .eq. 1) print *,'rthplot: th = ',thlast,thnext
         if((czprop(icpth) .eq. 'r') .or. (norefl .eq. 1) .or.
     .      (thlast .le. thnext)) goto 58
         thdif=thlast - thnext
         npt=max0(6,int(nval*thdif/thint))
         if(iirth .eq. 4) then
            do 85 n=0,npt
               th=thlast - thdif*(1.-cos(float(n)*pie/float(npt)))/2.
               apt=cos(th)/cref 
               call rdrcalc(apt,rth,dr,ibd,rbd,drbd)
               npoint=npoint + 1
               if(npoint .gt. npmax) then
                  print *,'Increase array sizes in rthplot!! ',
     .               'npoint = ',npoint
                  stop
               endif
               avals(npoint,nfun)=raysign*apt
               rvals(npoint,nfun)=rth
               call timatt(apt,e,tvals(npoint,nfun),ibd)
c: EKW 5-3-93: rtc_calc handles attenuation correctly:
ccc            call trcalc(apt,trm,trp,nps)
               call rtc_multi(apt,trm,trp,nps)
               ps=1. - (apt*cs)**2
               pr=1. - (apt*cr)**2
               gsl=sqrt(apt*cs*cr/(sqrt(ps*pr)*abs(rth*dr)))
               mvals(npoint,nfun)=expo(freq,e)*gsl*trm
85          continue
            npointmax=max0(npointmax,npoint)
         else
            nline=3 
            do 55 n=0,npt
               th=thlast - thdif*(1.-cos(float(n)*pie/float(npt)))/2.
               apt=cos(th)/cref 
               call rdrcalc(apt,rth,dr,ibd,rbd,drbd)
               if((rth .le. rhi) .and. (rth .ge. rlo)) then 
cxx               thpl=(thhi - th*piedeg)*hhfac
cxx               rpl=(rth-rlo)*vvfac
cxx               call plt(thpl,rpl,nline)
c: temporary write to file while plotting not available:
                  write(50,100) th*piedeg,rth
100               format(f8.4,1x,f10.2,1x,i1)
                  nline=2
               else 
                  nline=3
               endif
55          continue
         endif
58       continue
         do 50 kint=kstt+1,kendd
            kdirpt=mod(kdisc(kint),10)
cxx   print *,'kdirpt,kdir = ',kdir,kdirpt
            if(mod(kdir,kdirpt) .ne. 0) goto 50
            thlast=min(acos(aaxis(kaxlast)*cref),thmax)
            kaxnext=kdisc(kint)/10
            thnext=max(acos(aaxis(kaxnext-1)*cref),thmin)
            thdif=thlast - thnext
      if(iidiag .eq. 1) print *,'rthplot: kint = ',kint,thdif
            if(thdif .le. 0.) goto 50
            npt=max0(6,int(nval*thdif/thint))
            nline=3 
            if(iirth .eq. 4) then
               do 64 n=0,npt
                  th=thlast-thdif*(1.-cos(float(n)*pie/float(npt)))/2. 
                  apt=cos(th)/cref 
                  call rdrcalc(apt,rth,dr,ibd,rbd,drbd)
                  npoint=npoint + 1
                  if(npoint .gt. npmax) then
                     print *,'Increase array sizes in rthplot!! ',
     .                  'npoint = ',npoint
                     stop
                  endif
                  avals(npoint,nfun)=raysign*apt
                  rvals(npoint,nfun)=rth
                  call timatt(apt,e,tvals(npoint,nfun),ibd)
c: EKW 5-3-93: rtc_calc handles attenuation correctly:
ccc               call trcalc(apt,trm,trp,nps)
                  call rtc_multi(apt,trm,trp,nps)
                  ps=1. - (apt*cs)**2
                  pr=1. - (apt*cr)**2
                  gsl=sqrt(apt*cs*cr/(sqrt(ps*pr)*abs(rth*dr)))
                  mvals(npoint,nfun)=expo(freq,e)*gsl*trm
64             continue
               npointmax=max0(npointmax,npoint)
            else
               do 60 n=0,npt
                  th=thlast-thdif*(1.-cos(float(n)*pie/float(npt)))/2. 
                  apt=cos(th)/cref 
                  call rdrcalc(apt,rth,dr,ibd,rbd,drbd)
                  if((rth .le. rhi) .and. (rth .ge. rlo)) then 
                     thpl=(thhi - th*piedeg)*hhfac
                     rpl=(rth-rlo)*vvfac
c                 call plt(thpl,rpl,nline)
c: temporary write to file while plotting not available:
                     write(50,100) th*piedeg,rth
                     nline=2
                  else 
                     nline=3
                  endif
60             continue
            endif
            kaxlast=kaxnext
50       continue
c        call pltline(thpl,rpl,-.16)
c        write(8,200) (nttot(j),j=1,min0(10,ndp))
200      format(10(i2,1x))
         goto 44
c        call pltline(thpl,rpl-.20,-.16)
c        write(8,210) fdir(icpth),ntop,nbotx
210      format(a1,1x,i2,1x,i2)
c
10    continue
c     call pltline(dimh/2.,dimv+.15,-.16)
c     write(8,300) ktitle(1:lktit),zs,zr,cs,cr
300   format(a,'; zs = ',f7.2,'; zr = ',f7.2,'; cs = ',f8.3,
     .   '; cr = ',f8.3)
c     call pltend(0.0)
      if(iirth .eq. 4) then
c: For iirth=4, output is a matlab file of a vs. r and t:
         lenrec=8*(4*(npointmax*nfun) + 4*(1) + 1*(5*nfun)) + 
     .      1*(4*6 + 4*5 + 1*9) + 4*(4*5 + 4*5 + 1*5)
         print *,'lenrec = ',lenrec
         open(62,file=outn(1:loutn)//'_art.mat',access='direct',
     .      recl=lenrec,status='unknown')
c set MATLAB header values:
         mathead(1,1) = 1000       ! MATLAB type flag
         mathead(2,1) = npointmax  ! Number of rows in the matrix
         mathead(3,1) = nfun       ! Number of columns in the matrix
         mathead(4,1) = 0          ! Imaginary flag (=false)
         mathead(5,1) = 6          ! Length in bytes of variable name +1
         write(matfile(1)(6:6),'(a1)') char(0)
         matfile(1)(1:5)='avals'
         write(matfile(2)(6:6),'(a1)') char(0)
         matfile(2)(1:5)='rvals'
         write(matfile(3)(6:6),'(a1)') char(0)
         matfile(3)(1:5)='tvals'
         write(matfile(4)(6:6),'(a1)') char(0)
         matfile(4)(1:5)='mvals'
         mathead(1,2) = 1000       ! MATLAB type flag
         mathead(2,2) = 1          ! Number of rows in the matrix
         mathead(3,2) = 1          ! Number of columns in the matrix
         mathead(4,2) = 0          ! Imaginary flag (=false)
         mathead(5,2) = 5          ! Length in bytes of variable name +1
         write(matfile(5)(5:5),'(a1)') char(0)
         matfile(5)(1:4)='csrc'
         write(matfile(6)(5:5),'(a1)') char(0)
         matfile(6)(1:4)='crec'
         write(matfile(7)(5:5),'(a1)') char(0)
         matfile(7)(1:4)='zrec'
         write(matfile(8)(5:5),'(a1)') char(0)
         matfile(8)(1:4)='zsrc'
         mathead(1,3) = 1000       ! MATLAB type flag
         mathead(2,3) = 5          ! Number of rows in the matrix
         mathead(3,3) = nfun       ! Number of columns in the matrix
         mathead(4,3) = 0          ! Imaginary flag (=false)
         mathead(5,3) = 9          ! Length in bytes of variable name +1
         write(matfile(9)(9:9),'(a1)') char(0)
         matfile(9)(1:8)='raypaths'
         write(62,rec=1) (mathead(jcan,1),jcan=1,5),matfile(1)(1:6),
     .             ((avals(jcan,jcan2),jcan=1,npointmax),jcan2=1,nfun),
     .             (mathead(jcan,1),jcan=1,5),matfile(2)(1:6),
     .             ((rvals(jcan,jcan2),jcan=1,npointmax),jcan2=1,nfun),
     .             (mathead(jcan,1),jcan=1,5),matfile(3)(1:6),
     .             ((tvals(jcan,jcan2),jcan=1,npointmax),jcan2=1,nfun),
     .             (mathead(jcan,1),jcan=1,5),matfile(4)(1:6),
     .             ((mvals(jcan,jcan2),jcan=1,npointmax),jcan2=1,nfun),
     .             (mathead(jcan,2),jcan=1,5),matfile(5)(1:5),cs,
     .             (mathead(jcan,2),jcan=1,5),matfile(6)(1:5),cr,
     .             (mathead(jcan,2),jcan=1,5),matfile(7)(1:5),zr,
     .             (mathead(jcan,2),jcan=1,5),matfile(8)(1:5),zs,
     .             (mathead(jcan,3),jcan=1,5),matfile(9)(1:9),
     .      ((raypaths(jcan,jcan2),jcan=1,5),jcan2=1,nfun)
         close(62)
      endif
c
      return
      end 
