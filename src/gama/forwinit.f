c******************************************************************
       subroutine forwinit
c: for a given range-invariant ocean/layered bottom environment, this 
c: program calculates the acoustic pressure field using ray theory.
c: two input files are required: (1) the svp file and (2) the opt file.
c: the svp file contains the svp in the ocean and the geoacoustic
c: profile of the layered bottom.  the opt file contains the source/
c: receiver geometry, the ray paths to consider, and various plotting 
c: options.  graphical outputs include ray plots, propagation loss
c: plots, transfer function plots, and bar graphs of ray arrivals.
c: other outputs are eigenray lists and prop loss lists.

      implicit integer(i-n)
c From gamainit
      INCLUDE '../comopt.h'
      INCLUDE '../comforw.h'
c
      include './common/Parms_com'
      include './common/parms_conv'
c From GAMARAY
      include './common/timing' 
      include './common/saga_freqcom'
      include './common/srlox'
      include './common/gamaoptions'
      include './common/charcom'
      include './common/pii'
      include './common/depth'
      include 'common/bottom' 
      include './common/scatcom'
      include './common/scatarr'
      include './common/paths'
      include './common/causbad'
      include './common/laydata' 
      include './common/bdcom'
      include './common/tilt'
      include './common/saga_gama_fail'
      include './common/rthcom'
      include './common/saga_caustix'
cc
      complex*8 tf(NRTOTM*NFFTH1M),tfz(NFFTH1M),plcoh(NRNFMAX)
      complex*16 exparr(NFFTH1M),zzz(NRTOTM)
      real*4 ssf(2,NFFTH1M),ss2f(2,NFFTH1M),bsf(2,NFFTH1M),
     .   bs2f(2,NFFTH1M),plinc(NRNFMAX),phi(NRTOTM),
     .   angs(0:NASLM),frsl(0:NFSLM),dbsl(NAFSLM),phsl(NAFSLM),
     .   angb(0:NASLM),frbl(0:NFSLM),dbbl(NAFSLM),phbl(NAFSLM)
      real*8 frf(NFFTH1M),rangx(NRANGXM),tyme(NRANGXM),tmin(NRTOTM),
     .   dbmax(NRTOTM),range(NRTOTM),
     .   xxx(NRTOTM),yyy(NRTOTM),xlist(7*50*NRTOTM),r0(3,KAXMAX),
     .   dr0(3,KAXMAX),rbd0(3,KAXMAX),drbd0(3,KAXMAX),
     .   rlay(NLAYM,KAXMAX),drlay(NLAYM,KAXMAX),
     .   rbdlay(NLAYM,KAXMAX),drbdlay(NLAYM,KAXMAX),
     .   rc(KAXMAX),drc(KAXMAX),rbdc(KAXMAX),drbdc(KAXMAX),
     .   facax(NRTHMAX),rtmat(3*NRTHMAX)
      external time
      integer*4 time,jjnr(NRTOTM),kkk(NRTOTM),nda(NRTOTM,0:50),
     .   klist(2*50*NRTOTM),nlist(2*NRTOTM),ibd0(3,KAXMAX),
     .   ibdlay(NLAYM,KAXMAX),ibdc(KAXMAX),ntlay(-1:NLAYM+2,MBPMAXX)
      logical wk(NRTOTM),ap(NRTOTM),bintrf
c      real*4 etime,dtime,cptim,cpsec
      character*64 eline
      data eline/'INVALID INPUT IN OPT FILE: '/
c
      real*4 sngzs,c60
      character*2 chzs
      data c60/60./
cpln      INCLUDE 'gama_saga.h'

      integer*4 isectr,intrp(150),isubb
      real*4 rdstep
      integer ioer, izero, iilay
      integer noff,ndt
      common /info/ienv_info,lu_env,iiwrite
      common/iterpar/iter,iforwpop
      integer ienv_info,lu_env,iiwrite,iter,iforwpop
      integer ipos,j
      common /div/ipos
      integer*4 isec1
      real*8 rtmpold,zr1old
c
      if(iiwrite.gt.0)
     .  print *,'Running GAMARAY ...'
c
c
      iiwrite=-3
      if (iopt(6).eq.1) then
         iiwrite=10
         ienv_info=1         ! write the environment for each fm to file 90
         lu_env=90
      endif
C***********************************************
c
c
      entry forw2

c
c     Initialize Ocean ray paths
      ncp=0
      do ii=1,195
         do jj=1,3
            ncpth(ii,jj)=0
         end do
      end do
c
c---- If writing out trf function
c
        if (iWriteTrf.eq.1 .and. 
     &     dflast.gt.0 .and. fmaxlast.gt.fminlast) then
         dfrq=dflast
         fminsaga=fminlast
         fmaxsaga=fmaxlast
         nf1=max(1,nint(fminsaga/dfrq))
         nf2=nint(fmaxsaga/dfrq)
         fminsaga=nf1*dfrq
         fmaxsaga=nf2*dfrq
         nfbb=nf2 - nf1 + 1
         nfrq=nfbb
         nfftbb=2**(int(log10(float(nf2))/log10(2.) + .9))*2
         fs= dfloat(nfftbb)*dfrq
         write(*,*)' For optimised env Fmin, Fmax, DF, nfrq: ', 
     .        Fminsaga, Fmaxsaga, DFrq,nfrq,nf2,nfftbb
         if (nfrq.gt.mfreq) stop 'increase mfreq! '
         
         do j=1,nfrq
            frqsaga(j)=fminsaga + (j-1)*dfrq
         enddo           
      endif

c      write(*,*)'iiwrite=',iiwrite
c
c

      ctab=char(9)
cpln      isec1=time()
cpln      call idate(kdat(1),kdat(2),kdat(3))
cpln      call itime(ktim)
c
cpln      cptim=dtime(cpsec)
      nrangx=-nrangx
      flo=flotmp
      fhi=fhitmp
      nraytot=0
c
      if(ienv_info .NE. 1)
     &      iiwrite= iiwrite-1
      lwascomp=1
      iifail=0
C      print *,' Enter  PROSIM ...', lwascomp,iopt(17)
c
c---- Use of EOF
c
      if (iopt(17).eq.1) then
         call eofval
         if (lwascomp.eq.-1) then
            return
         endif
      endif
      nsctor=1
c
      if (ipos.eq.1) then
         do isectr=1,nsctor
c the impedance contrast should be positive
            if(R_H1(isectr) .ne. 0.0) then
               if (r_r1(1,isectr).gt.r_r2(isectr)) then
                  lwascomp=-1
                  write(*,*)' Model',iforwpop,' rejected. dens-contr'
                  return
               elseif (R_c1(R_ND1(ISECTr),isectr).gt.
     .                (R_c2(isectr))) then
                  lwascomp=-1
                  write(*,*)' Model',iforwpop,' rejected. Csed > Cbot'
                  return
               endif
c     for sediment profiles with negative slope of ssp
               do iilay=2,R_ND1(ISECTR)
                  if(R_c1(iilay-1,isectr).gt. 
     .                 R_c1(iilay,isectr)) then
                     lwascomp=-1
                     write(*,*)' Model',iforwpop,' rejected. neg slope'
                     write(*,*)' between layer',iilay-1,' and ',iilay
                     return
                  end if
               end do
            endif
         enddo
      endif
      if (ipos.eq.2) then
         do isectr=1,nsctor
            if(R_c2(isectr).lt.1530)then
               write(*,*)' Model',iforwpop,' rejected. Cbot<1530'
               lwascomp=-1
               return
            endif
         enddo
      endif
c
      isubb = 1
c: Set iiwrite=1 so that output files and messages are sent:
c      
      do isectr=1,nsctor
         R_z0(R_nd0(isectr),isectr)=R_h0(isectr) ! last point is the water depth
         if (R_nd1(isectr).gt.0) then
            R_z1(R_nd1(isectr),isectr)=R_h1(isectr)
c            write(6,*)'R_z1,R_h1,R_nd1: ',R_z1(R_nd1(isectr),isectr),
c     .           R_h1(isectr),R_nd1(isectr)
c            pause
         end if
      enddo
      
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
      if(iifail.eq.1) then
         lwascomp=-1
         return
      end if

      nlayer=nlay
      call mem_lim(nasl,NASLM,eline,27,'NASLM',5,iibad)
      call mem_lim(nfsl,NFSLM,eline,27,'NFSLM',5,iibad)
      call mem_lim(nasl*nfsl,NAFSLM,eline,27,'NAFSLM',6,iibad)
      call mem_lim(nrth,NRTHMAX,eline,27,'NRTHMAX',7,iibad)
      if(iibad .eq. 1) stop 'Stop in forwinit #1'

      IF (iopt(6).eq.1)then         ! option Q.
         call write_info
         call write_gama_svp
         call write_gama_opt
      ENDIF
c
      if(iiwrite.gt.0)
     .   print *,'Running GAMARAY ...'
c
cpln      call flush(6)
c
      do isectr=1,nsctor 

c INSERT GAMA HERE
         outn='gamaout'
         loutn=7
         nrayt=0
c     
c     Initialize parameters
c
         do jcan=1,NRTOTM
            xxx(jcan)=0.
            yyy(jcan)=0.
            zzz(jcan)=dcmplx(0.d0,0.d0)
            phi(jcan)=0.
            jjnr(jcan)=0
            kkk(jcan)=0
         end do
c
         if(iiffi .ne. 0) then
            if(iisl .eq. 0) then
               do 940 jcan=1,nffth1
                  ssf(1,jcan)=1.
                  ssf(2,jcan)=0.
                  ss2f(1,jcan)=1.
                  ss2f(2,jcan)=0.
 940           continue
            endif
            if(iibl .eq. 0) then
               do 942 jcan=1,nffth1
                  bsf(1,jcan)=1.
                  bsf(2,jcan)=0.
                  bs2f(1,jcan)=1.
                  bs2f(2,jcan)=0.
 942           continue
            endif
            do 10 j=1,nffth1
               frf(j)=float(j-1)*dfrq
 10         continue
         endif
         if(iisl .eq. 0) then
            do 944 jcan=1,nfr2
               ss(1,jcan)=1.
               ss(2,jcan)=0.
               ss2(1,jcan)=1.
               ss2(2,jcan)=0.
 944        continue
         endif
         if(iibl .eq. 0) then
            do 946 jcan=1,nfr2
               bs(1,jcan)=1.
               bs(2,jcan)=0.
               bs2(1,jcan)=1.
               bs2(2,jcan)=0.
 946        continue
         endif
c
         call geom(rangx,tyme,range,jjnr,kkk,xxx,yyy,phi)
c: EKW: 5-3-93
cxx      call pastri
         call pasc_make
cpln      call saga_getdat3(nline1)
c
         ii12=0
         if((iipl .gt. 0) .or. (iitf .gt. 0)) then
            ii12=1
            open(12,file=outn(1:loutn)//'.ls',status='unknown')
            rewind(12)
cpln            call prhead(12)
         endif
c
c
      rtmpold=rtmp(1)
      zr1old=zr1
      do 20 jzs=1,nzs
c: construct different names for ir and fft files for each zs:
         if(iiffi .gt. 0) then
            if(nzs .gt. 1) then
               chzs(1:1)=char(48+jzs/10)
               chzs(2:2)=char(48+mod(jzs,10))
               lchzs=2
            else 
               lchzs=0
            endif
            tfmax=0.
            tfmin=1.e13
         endif
         if(iiir .eq. 1) then
            if(jzs .gt. 1) close(50)
            open(50,file=outn(1:loutn)//chzs(1:lchzs)//
     .         '_ir',status='unknown',form='unformatted')
            rewind(50)
         endif
c         if(iifft .ne. 0) then
c            if(jzs .gt. 1) close(16)
cpln            open(16,file=outn(1:loutn)//chzs(1:lchzs)//
cpln     .         '_fft',access='direct',recl=NDAW*(22+nfft),
cpln     .           status='unknown')
c            open(16,file=outn(1:loutn)//chzs(1:lchzs)//
c     .         '_fft',status='unknown',
c     .           form='formatted')
c            open(16,file=outn(1:loutn)//
c     .         '_fft',status='unknown',
c     .           form='formatted')
c            rewind(16)
c         endif
         if(iibmp .eq. 1) then
            sngzs=zs
            write(18) sngzs
         endif         
         call newsdep(nlist,plcoh,plinc,nda,tmin,dbmax,range,NDAW)
         ncbad=0
         kph=0
         ntcut=0
         nbdbad=0
c: plot ray picture axes:
         if(iipic .ge. 1) call plotpic(range,jjnr)
         do 30 jzr=1,nzr
c
            if(tilth.or.tiltv) then
               nrangx=-1
c               rtmp(1)=rangx(jzr)*0.001
c               rtmp(1)=dble(jzr-1)*drtmp+rangx(1)*0.001
               rtmp(1)=dble(jzr-1)*drtmp+rtmpold
               call getdat2tilt
               call geomtilt(rangx,tyme,range,jjnr,kkk,xxx,yyy,phi)
               call pasc_make
               call getdat3tilt
c               pause
            end if
            if (iiwrite.gt.0) then
               write(6,*)'zr0,zr1,dz0: ',zr0,zr1,dz0(1)
               write(6,*)'drtmp,drtmp0: ',drtmp,drtmp0
               write(6,*)'nzr,nrangx: ',nzr,nrangx
               write(6,*)'jzr, rtmp(1): ',jzr,rtmp(1)
               write(6,*)'jzr,rangx,range: ',jzr,rangx(jzr),range(jzr)
               write(6,*)'jzr,rangx,range: ',jzr,rangx(1),range(1)
               write(6,*)'xyz(1,3),dz(1): ',xyz(1,3),dz(1)
               write(6,*)'Horizontal tilt: ',tilth,dtilth
               write(6,*)'Vertical tilt: ',tiltv,dtiltv
               write(6,*)
            endif
c
            call newrdep(range)
c: determine sample points in the snell invariant (si) axis.
            call discpts(kax)
            call mem_lim(kax,KAXMAX,eline,27,'KAXMAX',6,iibad)
            if(iibad .eq. 1) stop 'Stop in forwinit #2'
c: calculate r(a) and r'(a) functions in ocean and bottom layers.
            call ocean(r0,dr0,ibd0,rbd0,drbd0)
            call layers(rlay,drlay,ibdlay,rbdlay,drbdlay)
            if(iipic+iibet+iipltf+iiffi .eq. 0) goto 29
c
c: find eigenrays:
            call rays(xlist,klist,nlist,range,wk,ap,kkk,xxx,yyy,nda,
     .         tmin,dbmax,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl,
     .         r0,dr0,ibd0,rbd0,drbd0,rlay,drlay,ibdlay,rbdlay,drbdlay,
     .         ntlay,NLAYM,rc,drc,ibdc,rbdc,drbdc,KAXMAX)
            call prwarn
c
29          continue
            if(iirth .ne. 0) call rthplot(iirth,facax,rtmat)
c: plot geoacoustic parameters on ray picture.
            if(iipic .ge. 1) call acouspl
30       continue
c
         cptim=etime(cpsec)
         write(8,300) zs,int(cptim/c60),int(mod(cptim,c60))
300   format(/'FINISHED EIGENRAY SEARCH FOR SOURCE DEPTH ',f7.2,
     .   ';  CP MINUTES = ',i4,':',i2.2)
c: calculate field at receivers and output desired ray lists or fft's: 
         call srcstep(xlist,klist,nlist,tf,tfz,frf,ssf,ss2f,bsf,bs2f,
     .      plcoh,plinc,rangx,tyme,range,phi,jjnr,nda,tmin,dbmax,
     .      angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl,exparr,
     .      ntlay,NLAYM)
c
cpln         write(8,310) zs,int(cptim/c60),int(mod(cptim,c60))
cpln 310   format(/'FINISHED EIGENRAY PROCESSING FOR SOURCE DEPTH ',f7.2,
cpln     .   ';  CP MINUTES = ',i4,':',i2.2)
c: make prop loss or transfer function plots:
c         call plcalls(range,tyme,plcoh,plinc,xxx,yyy,zzz,jjnr,
c     .      angs,frsq,dbsl)
cpln         cptim=etime(cpsec)
         if(iifft .ne. 0) then
            if(tfmax .ne. 0.) tfmax=20.*log10(tfmax)
            if(tfmin .ne. 0.) tfmin=20.*log10(tfmin)
            write(8,400) tfmin,tfmax
400         format(/'MIN,MAX VALUES IN TRANSFER FUNCTIONS = ',f7.2,2x,
     .         f7.2,' dB')
         endif
         close(54)
         if(iiir .ne. 0) close(50)
         if(iifft .ne. 0) then
           close(16)
           close(51)
         end if
 20   continue
c
      if(iieig .ne. 0) close(59)
      if(ii12 .eq. 1) close(12)
      if(iibmp .ne. 0) close(18)
      if(iidat .ne. 0) close(17)
      if(iirth .eq. 1 .or. iirth .eq. 2) close(50)
      if(iipic .ne. 0) then
         close(55)
         close(51)
      endif
      close(56)
      close(57)
      close(8)
      nraytot=nrayt
      rtmp(1)=rtmpold
      zr1=zr1old
c
      if(iifail .eq. 1) then
         print *,'FAILURE OCCURRED FOR THIS GAMARAY RUN.'
         lwascomp=-1
cpln         stop
         return
      endif

      call saga_out(tf)

      enddo                     !FINISH SUBBAND CALCULATIONS
c
cpln      cptim=dtime(cpsec)
cpln      isec2=time() - isec1
cpln      cptim=etime(cpsec)
cpln      print *,'isec1,isec2 = ',isec1,isec2,time()
cpln      write(8,330) int(cptim/c60),int(mod(cptim,c60)),
cpln     .   int(isec2/60),int(mod(isec2,60)),nraytot
cpln330   format(//'TOTAL CP MINUTES TAKEN FOR RUN = ',i4,':',i2.2,
cpln     .   '; WALL CLOCK TIME = ',i4,':',i2.2,/
cpln     .   'TOTAL NUMBER OF RAYS PROCESSED = ',i6)
cpln      close(8)


      if ((iWriteTrf.gt.0)) then

         write(6,*)'Calling bb_fft_out for best model'
c         open(16,file=outn(1:loutn)//
c     .        '.trf',status='unknown',
c     .        form='formatted')
         rewind(16)
         call bb_fft_out(rangx,tf)
         close(16)
      endif
cpln      rewind(16)
cpln      call bb_fft_out(rangx,tf)
cpln      close(16)
cpln      pause
c
      isubb=1
               
      return
      
      end
c
