      subroutine srcstep(xlist,klist,nlist,tf,tfz,frf,ssf,ss2f,bsf,bs2f,
     .   plcoh,plinc,rangx,tyme,range,phi,jjnr,nda,tmin,dbmax,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl,exparr,ntlay,NLAYM)
c
c: this subroutine steps through the source positions and calls
c: raycalc for each receiver in the array.
      implicit integer*4(i-n)
c
      include 'common/gamaoptions'
      include 'common/srlox'
      include 'common/depth'
      include 'common/freqcom'
      include 'common/svp'
      include 'common/charcom'
      include 'common/paths'
      complex*8 plcoh(nrtot,nfr2),tf(nzr*nrangx*nfbb),tfz(nffth1)
      real*4 xhead(50),xh2(20),plinc(nrtot,nfr2),phi(nrtot),sngzr,sngrng
      real xlist(7,mray,nrtot)
      integer*4 klist(2,mray,nrtot),nlist(2,nrtot)
      real tmin(nrtot),dbmax(nrtot),frf(nffth1)
c: count on less than 2500 rays for each s/r configuration: 
      real rangx(nrangx),tyme(nrangx),range(nrtot),xl(7,2500)
      integer*4 jjnr(nrtot),nda(nrtot,0:50),kl(2,2500),khead(50)
      equivalence(khead,xhead)
      data xh2/2.,0.,0.,2.,2.,2.,0.,0.,0.,0.,
     .         2.,1.,1.,0.,1.,0.,0.,0.,0.,0./
c
c: eliminate the 3-dB margin for error introduced in rays:
c     cutmag=10**(cutdb/-20.)
c
c     print *,'# direct access writes: ',nnda,' # kbytes ',3.201*nnda
c: initialize header for fft or ir files: 
      if(iiffi .ne. 0) call csfft1d(nfft,tfz,0)
      if(iiffi .ne. 0 .or. iidat .ne. 0) then
c: xhead,khead is header for impulse response file; xh2 is header for
c: fft file:
         xhead(3)=float(nfft)
         xhead(4)=fs
         xhead(6)=float(nrec)
         xh2(5)=float(nfft) 
         xh2(6)=fs
         xh2(8)=float(nzr)
         xh2(10)=nffth1
         xh2(11)=float(nrec)
         xh2(12)=1.
c: EKW uses this for image, xcorr program:
         xh2(16)=csarr(jzs)
         xh2(19)=zs
         xh2(20)=float(nrangx)
         xhead(29)=zs
         jrec=0
         do 110 jzr=1,nzr
            do 120 jxy=1,nxy(jzr)
               jrec=jrec + 1
               if(jrec .gt. 21) goto 112
               xhead(6+jrec)=xyz(kxy(jzr)+jxy,3)
               xhead(29+jrec)=xyz(kxy(jzr)+jxy,1)
120         continue
110      continue
112      continue
      endif
      irec=0
      do 10 ksrc=1,nrangx
         ksrc1=ksrc-1
         rsrc=rangx(ksrc)
         tsrc=tyme(ksrc)
         tmall=0.
         if(iiffi .ne. 0 .or. iidat .ne. 0) then
            if(iifft .le. 3) then
               tmall=1.e13
               do 50 jzr=1,nzr
                  kref=mr1(jzr) + nxy(jzr)*ksrc1 - 1 
                  do 55 jxy=1,nxy(jzr)
                     tmall=min(tmall,tmin(jjnr(kref+jxy)))
55                continue
50             continue
               if(tmall .eq. 1.e13) then
                  tmall=0.
               else
                  tmall=tmall - .04*nfft/fs
               endif
            endif
            xhead(2)=rangx(ksrc)
            xh2(3)=float(ksrc)
            if(iiir .eq. 1) then
               write(50) ktitle(1:8)
               write(50) (khead(jcan),jcan=2,50)
            endif
         endif
         jrec=0
         do 20 jzr=1,nzr
            zr=xyz(kxy(jzr)+1,3)
            cr=crarr(jzr)
c: EKW uses this for image, xcorr program:
            xh2(17)=cr
            kref=mr1(jzr) + nxy(jzr)*ksrc1 - 1
            do 30 jxy=1,nxy(jzr)
               mr=jjnr(kref + jxy)
               rng=range(mr)
               tminsr=tmin(mr)
c: For iifft=4,5,6 set different reference times for each receiver:
               if(iifft .gt. 3) tmall=tmin(mr) - .04*nfft/fs
               xh2(1)=tmall
               xhead(5)=tmall
               dbsr=dbmax(mr) 
               nray=0
c     print *,'mr,nda(mr,0),nda = ',mr,nda(mr,0),nda(mr,1:nda(mr,0))
               do 35 kda=1,nda(mr,0)
                  read(54,rec=nda(mr,kda)) nrr,((xl(jcan,jcan2),
     .               jcan=1,7),jcan2=nray+1,nray+nrr),
cxx  .               ((kl(jcan,jcan2),jcan=1,2),jcan2=nray+1,nray+nrr)
c: Bug found by Cederberg:
     .                ((kl(jcan,jcan2),jcan2=nray+1,nray+nrr),jcan=1,2)
                  nray=nray + nrr
35             continue
               nrnow=nlist(1,mr)
               do 940 jcan=1,nrnow
                  do 942 jcan2=1,7
                     xl(jcan2,nray+jcan)=xlist(jcan2,jcan,mr)
942               continue
                  kl(1,nray+jcan)=klist(1,jcan,mr)
                  kl(2,nray+jcan)=klist(2,jcan,mr)
940            continue
               nray=nray + nrnow
               jrec=jrec + 1
               xh2(4)=float(jrec)
c: EKW uses this for image, xcorr program:
               xh2(13)=xyz(kxy(jzr)+jxy,1)
               xh2(14)=xyz(kxy(jzr)+jxy,2)
               xh2(15)=xyz(kxy(jzr)+jxy,3)
               xh2(18)=rng
c: write to beam pattern file:
               if(iibmp .eq. 1) then
                  sngzr=zr
                  sngrng=rng
                  write(18) sngzr,sngrng,phi(mr),nray
               endif
c: EKW (3-25-92): Need reinitialize SVP stuff if plotting rays:
               if(iipic .ne. 0) call newrdep(range)
               call raycalc(mr,ksrc,xl,kl,nray,tf,tfz,frf,ssf,ss2f,bsf,
     .            bs2f,plcoh,plinc,range,xh2,jrec,
     .            angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl,exparr,
     .            ntlay,NLAYM,irec)
               if(iidat .ne. 0) then
c: note that xh2 is real*4; other reals are real*8:
                  write(17) nray,(xh2(jcan),jcan=1,20),((xl(jcan,jcan2),
     .               jcan=1,7),jcan2=1,nray),
     .               ((kl(jcan,jcan2),jcan=1,2),jcan2=1,nray)
               endif
30          continue
20       continue
10    continue
c
      return
      end 
