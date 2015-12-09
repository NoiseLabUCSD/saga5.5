      subroutine rays(xlist,klist,nlist,range,wk,ap,kkk,xxx,yyy,nda,
     .   tmin,dbmax,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl,
     .   r0,dr0,ibd0,rbd0,drbd0,rlay,drlay,ibdlay,rbdlay,drbdlay,
     .   ntlay,NLAYM,rc,drc,ibdc,rbdc,drbdc,KAXMAX)
c
c: this subroutine finds the eigenrays between the source and receivers
c
      implicit integer*4(i-n)
      include 'common/dublcom'
      include 'common/srlox'
      include 'common/laydata'
      include 'common/paths'
      include 'common/pathway'
      include 'common/depth'
      include 'common/discpt' 
      include 'common/pathchr'
      include 'common/weakcom'
      include 'common/discint'
      include 'common/hccom'
      include 'common/gamaoptions'
      include 'common/bdcom'
      include './common/saga_gama_fail'
      integer*4 nlist(2,nrtot),kkk(nrtot),kf(2),KAXMAX
      real range(nrtot),xxx(nrnf),yyy(nrnf)
      real amono(3),amonox(3),aprev(3),anext(3),a1(3,2),a2(3,2)
      real r0(3,inx),dr0(3,inx),rbd0(3,inx),drbd0(3,inx)
      real rc(KAXMAX),drc(KAXMAX),rbdc(KAXMAX),drbdc(KAXMAX)
      integer*4 ibd0(3,inx),ibdc(KAXMAX)
      logical wk(nrtot),ap(nrtot)
      integer*4 ntlay(-1:NLAYM+2,50000)
      common /info/ienv_info,lu_env,iiwrite
      integer ienv_info,lu_env,iiwrite
c
c: allow a 3-dB margin for error in cutting rays due to weakness:
c     cutmag=10**((cutdb+3)/-20.)
      nbwk=0
c: initialize mrmin(nbotx) [for r<range(mrmin(nbotx)), rays with
c: nbot>nbotx are too weak to include.] 
      do 5 nb=0,nbmax+1
         mrmin(nb)=m22
5     continue
c
      icxx=-1
c: BUG discovered 4-2-90: Initialize mbp=0 in newsdep only:
      mbp=0

c: loop for the ncp number of ocean paths
      do 10 icpth=1,ncp
c
         ntop=ncpth(icpth,1)/2
         nbotx=ncpth(icpth,3)/2
         nttot(0)=2*nbotx
         nbotx1=nbotx+1
c: if rays have been found for previous bottom bounce: 
         if(mrmin(nbotx) .ne. m22) then
            mmr0b=mrmin(nbotx)
         else
            mmr0b=mmr1 
         endif
         if(mmr0b .gt. mmr2) then
            if(iiwrite.gt.0)
     .       write(6,110) nbotx
110   format('SKIPPING OCEAN PATHS WITH NBOT >= ',i3) 
            goto 99
         endif
c
c: add up the r(a) and r'(a) functions for the ocean path.
         do 20 k=1,inx 
            rc(k)=0.
            drc(k)=0.
            rbdc(k)=0.
            drbdc(k)=0.
            ibdc(k)=0
            do 30 icleg=1,3
               rc(k)=rc(k) + r0(icleg,k)*ncpth(icpth,icleg)
               drc(k)=drc(k) + dr0(icleg,k)*ncpth(icpth,icleg)
               if(ibd0(icleg,k) .eq. 1) then
                  ibdc(k)=1
                  rbdc(k)=rbdc(k) + rbd0(icleg,k)*ncpth(icpth,icleg)
                  drbdc(k)=drbdc(k) + drbd0(icleg,k)*ncpth(icpth,icleg)
               endif
30          continue
20       continue
c
         call setop(kst,kend,kdir,czprop(icpth))
c
44       continue
c: initialize mmr0 to value for current ocean path: 
cpln         write(6,*)
cpln         write(6,*)'RAYS NLAYM: ',NLAYM
         call setbp(kst,kstt,kend,kendd,norefl,czprop(icpth),kdone,
     .      mmr0b,ntlay,NLAYM)
         if(czprop(icpth) .eq. 'p' .and. nttot(1) .eq. 0) goto 44
         if(kdone .eq. 1) goto 10
         if(iidiag .ne. 0) print *,'new ray path: icpth,mbp,norefl = ',
     .      icpth,mbp,norefl
c
c: krmin is the min kr for which a strong ray is found for this path: 
         krmin=m22
c: krweak is the max kr for which a weak ray is found for this path:
         krweak=0
c: check the ranges from a=0 to the lowest discontinuity point.
c: the graph is increasing in this region.
         kaxlast=kdisc(kstt)/10
         do 48 kr=mmr1,mmr2
            nlist(2,kr)=nlist(1,kr)
48       continue
c: for refracting ocean paths (which do not touch the ocean bottom),
c: do not check for eigenrays which interact with the bottom: 
c: also, do not check for rays reflecting off interface w/ no mismatch: 
         if((czprop(icpth) .ne. 'r') .and. (norefl .eq. 0)) then
            kbdf=1
            aprev(1)=0.
            aprev(2)=0.
            aprev(3)=1.
            anext(1)=aaxis(kaxlast-1)
            call rdrsum(kaxlast-1,anext(2),anext(3),rc,drc,ibdc,
     .         rbdc,drbdc,rlay,drlay,ibdlay,rbdlay,drbdlay)
            nlom=0
cpln            write(6,*)'rays: rcheck #1'
            call rcheck(aprev,anext,nlom,xlist,klist,nlist,
     .         range,wk,ap,kkk,xxx,yyy,nda,tmin,dbmax,
     .         angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
            if(iifail.eq.1) return
            call maxup(xlist,klist,nlist,tmin,dbmax)
         endif
         do 940 jcan=1,3
            adlo(jcan)=aprev(jcan)
            adhi(jcan)=anext(jcan)
940      continue
c
         iihc=0
c: check the remaining intervals between discontinuity points.
         do 50 kint=kstt+1,kendd
            kdirpt=mod(kdisc(kint),10)
c: check if the next disc pt is a disc pt for the current ocean path. 
            if(mod(kdir,kdirpt) .ne. 0) goto 50
c
c: update frequency for bd calcs:
            kbdf=1
            mmbd=0
150         if(iibd .eq. 1) then
               bdfreq=bdf(kbdf)
               bdomega=bdom(kbdf)
               bdfac=bdf(1)/bdfreq - 1.
            endif
            rmn=max(0.,range(mmr0) - 10.*bdlam(kbdf))
            rmx=range(mmr2) + 10.*bdlam(kbdf)
      if(iidiag .ne. 0) print *,'disc int: kint,kbdf = ',kint,kbdf
c
c: initialize variables for checking for caustic doublet pairs: 
            nlom=0
c: kc is the pointer for the cautics found, alternates between 1 and 2:
            kc=2
            do 52 kr=mmr1,mmr2
               nlist(2,kr)=nlist(1,kr)
52          continue
            amono(1)=aaxis(kaxlast)
            call rdrsum(kaxlast,amono(2),amono(3),
     .         rc,drc,ibdc,rbdc,drbdc,rlay,drlay,ibdlay,rbdlay,drbdlay)
            do 942 jcan=1,3
               aprev(jcan)=amono(jcan)
               adlop(jcan)=adlo(jcan)
               adhip(jcan)=adhi(jcan)
942         continue
            kaxnext=kdisc(kint)/10
c: go through the aaxis points in the interval and check for sign
c: changes in dr.
            do 60 kax=kaxlast+1,kaxnext-1
               anext(1)=aaxis(kax)
               call rdrsum(kax,anext(2),anext(3),
     .      rc,drc,ibdc,rbdc,drbdc,rlay,drlay,ibdlay,rbdlay,drbdlay)
c: when a sign change in dr occurs, loma will find the pt where dr=0. 
               if(aprev(3)*anext(3) .le. 0.) then
                  kd=kc
                  kc=3 - kc
cpln                  write(6,*)'LOMA #1'
                  call loma(aprev,anext,a1(1,kc),a2(1,kc),kf(kc),nlom)
cpln                  write(6,*)'rays: rcheck #2'
                  call rcheck(amono,a1(1,kc),nlom,
     .               xlist,klist,nlist,range,wk,ap,kkk,xxx,yyy,nda,
     .               tmin,dbmax,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
                  if(iifail.eq.1) return
c: check for double extremum discovered: 
                  if(kbad .eq. 1) then
                     nlom=nlom-1
                     call dubloma(amono,a1(1,kc),amonox,
     .   a1(1,kd),a2(1,kd),kf(kd),nlom,xlist,klist,nlist,range,wk,ap,
     .   kkk,xxx,yyy,nda,tmin,dbmax,angs,frsl,dbsl,phsl,
     .   angb,frbl,dbbl,phbl)
                     if(iifail.eq.1) return
cpln            write(6,*)'rays: rcheck #3'
                     call rcheck(a2(1,kd),a1(1,kc),nlom,
     .   xlist,klist,nlist,range,wk,ap,kkk,xxx,yyy,nda,tmin,dbmax,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
                     nlom=nlom+1
                     if(iifail.eq.1) return
                  endif
                  if(nlom .gt. 1) then 
         call shadow(a1(1,kd),a2(1,kd),kf(kd),
     .   amonox,a1(1,kc),iifc,xlist,klist,nlist,range,wk,ap,nda,
     .   tmin,dbmax,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
                  end if
                  do 944 jcan=1,3
                     amonox(jcan)=amono(jcan)
                     amono(jcan)=a2(jcan,kc)
944               continue
c: check for double kink, even though r's appear monotonic: 
               elseif((anext(2)-aprev(2))*anext(3) .lt. 0.) then
                  call drswich(aprev,anext)
                  if(iifail.eq.1) return
                  call dubloma(amono,anext,amonox,a1(1,kc),a2(1,kc),
     .   kf(kc),nlom,xlist,klist,nlist,range,wk,ap,kkk,xxx,yyy,nda,
     .   tmin,dbmax,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
                  if(iifail.eq.1) return
               endif
               aprev(1)=anext(1)
               aprev(2)=anext(2)
               aprev(3)=anext(3)
60          continue
c: check the final range interval.
cpln            write(6,*)'rays: rcheck #4'
            call rcheck(amono,anext,nlom,xlist,klist,nlist,range,wk,ap,
     .         kkk,xxx,yyy,nda,tmin,dbmax,
     .         angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
            if(iifail.eq.1) return
c: check for double extremum discovered: 
            if(kbad .eq. 1) then
               call dubloma(amono,anext,amonox,a1(1,kc),
     .   a2(1,kc),kf(kc),nlom,xlist,klist,nlist,range,wk,ap,kkk,xxx,yyy,
     .   nda,tmin,dbmax,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
               if(iifail.eq.1) return
cpln            write(6,*)'rays: rcheck #5'
               call rcheck(a2(1,kc),anext,nlom,xlist,klist,nlist,
     .            range,wk,ap,kkk,xxx,yyy,nda,tmin,dbmax,
     .            angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
               if(iifail.eq.1) return
            endif
cxx         if((kint .eq. kendd) .and. (ndp .eq. 0) .and. 
cxx  .            (nldif .ne. 0)) then
cxx            iihc=1
cxx   print *,'iihc set to 1: kint,kendd,nldif = ',kint,kendd,nldif
cxx   print *,'amono,anext = ',1./amono(1),1./anext(1)
cxx            kpair=(3+ncpth(icpth,nldif)-ncpth(icpth,2))/2
cxx            kp0=(3+ncpth(icpth,nlsam)-ncpth(icpth,2))/2
cxx            iifc=0
cxx         endif
            if(nlom .gt. 0) then
         call shadow(a1(1,kc),a2(1,kc),kf(kc),
     .   amonox,anext,iifc,xlist,klist,nlist,range,
     .   wk,ap,nda,tmin,dbmax,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
         end if
c
c: NEW 4-18-90: check for horizontal caustic pairs:
cxx         if(iihc .eq. 1) then
c     print *,'iihc=1: kp0,kpair,mbp,iifc = ',kp0,kpair,mbp,iifc
cxx            iifcs(kp0,kpair)=iifc
cxx            if(kpair .eq. 1) then
cxx               mbphc(kp0)=mbp
cxx            else
cxx               call horcaus(mbp,xlist,klist,nlist)
cxx            endif
cxx         endif
            call maxup(xlist,klist,nlist,tmin,dbmax)
c
c: go back and find eigenrays with bd at new frequency:
            if((mmbd .ne. 0) .and. (kbdf .lt. nbdf)) then
               kbdf=kbdf + 1
               goto 150
            endif
c
            do 946 jcan=1,3
               adlo(jcan)=amono(jcan)
               adhi(jcan)=anext(jcan)
946         continue
            kaxlast=kaxnext
50       continue
c
c: enter nttot bottom path in mweak if weak rays found: 
c: krmw points to the shortest range for which bottom path needs to be
c: considered:
         if(krmin .ne. m22) then
c: we know some strong rays were found:
            if(krweak .ne. 0) then
c: some weak rays were found:
               krmw=min0(krmin,krweak)
            else
c: no weak rays were found:
               krmw=krmin
            endif
         else
c: no strong rays were found:
c: if no weak rays either, goto next bottom path:
            if(krweak .eq. 0) goto 144
c: set krmw to the max range for which weak ray found:
            krmw=krweak
         endif
c
         if(krmw .le. mmr0) goto 144
c: don't enter path if all Mj are equal (because such paths have only
c: one multipath and can be unexpectedly weaker than similar paths):
         if(ndp .eq. 0) goto 69
         do 68 j=1,ndp
            if(nttot(j) .ne. nttot(0)) goto 69
68       continue
         goto 144
69       continue
c
c: check if nttot already in list: 
         do 70 jwk=nbwk,1,-1
c: skip rest of checks if nbotx already higher:
            if(mweak(jwk,3) .lt. nbotx) goto 75
            mw=mweak(jwk,1)
c: skip this weak bottom path if ndp or Mj not same:
            if(ndp .ne. ntlay(-1,mw)) goto 70
            do 72 j=1,ndp
               if(nttot(j) .ne. ntlay(j,mw)) goto 70
72          continue
c: nttot found in list as entry # jwk. update min kr to consider: 
c     print *,'jwk,mweak,krmw = ',jwk,mweak(jwk,2),krmw
            mweak(jwk,2)=max0(krmw,mweak(jwk,2))
            goto 144
70       continue
75       nbwk=nbwk+1
         if(nbwk .gt. 600) then
      print *,'nbwk exceeded 600 in rays. increase array dim.'
            do 82 j=1,99
               j99=99 + j
               do 948 jcan=1,3
                  mweak(j,jcan)=mweak(j99,jcan)
948            continue
82          continue
            nbwk=100
         endif
c     print *,'mbp,krmw = ',mbp,krmw
         mweak(nbwk,1)=mbp
         mweak(nbwk,2)=krmw
         mweak(nbwk,3)=nbotx
144      goto 44
c
10    continue
99    continue
c
      return
      end 
