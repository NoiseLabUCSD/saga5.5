      subroutine rayplot(a,r,jcaus,ntlay,NLAYM)
c
c: this subroutine plots the eigenray on the ray picture.
c
      implicit integer*4(i-n)
      include 'common/svp'
      include 'common/bottom' 
      include 'common/pic'
      include 'common/srlox'
      include 'common/depth'
      include 'common/paths'
      include 'common/pathchr'
      include 'common/pathway'
      include 'common/gamaoptions'
      include 'common/charcom'
      include 'common/plotcom'
c
      real dro(240),dzo(240),raydr(240,3),raydz(240,3),drb(240,42),
     .   dzb(240,42)
      integer*4 line(0:3),nbounce(0:43,200),nbleg(200),jup(3),iup(3)
      integer*4 npt(3),nptb(0:42)
      integer*4 ntlay(-1:NLAYM+2,mbp)
      data line/2,4,9,6/,jup/1,2,2/,iup/1,-1,-1/
      
c
      ndp=ntlay(-1,mbp)
c
      ainv=1./a
c: nleg is the total number of traversals of the ocean sections;
c: iiup=1 means upward, iiup=-1 means downward; ldir is the last
c: direction of the ray (u or d), leg is the first ocean section 
c: from the receiver.
      nleg=ncpth(icpth,1) + ncpth(icpth,2) + ncpth(icpth,3)
      iiup=1
      if(ldir(icpth) .eq. 'u') iiup=-1
      if((iiup .eq. 1) .and. (zr .lt. zs)) then
         leg=1
      elseif((iiup .eq. -1) .and. (zr .ge. zs)) then
         leg=3
      else
         leg=2
      endif
c
c: figure out an ordering of the bottom paths: ntlay holds the total
c: # of legs in each layer to be divided among nbotx bottom bounces:
c: first, divide up legs evenly among bottom bounces:
      do 940 jcan=1,nbotx
         do 942 jcan2=0,ndp+1
            nbounce(jcan2,jcan)=0
942      continue
940   continue
      do 110 jl=1,ndp
         nptb(jl)=0
         jb=0
         do 120 j=1,ntlay(jl,mbp)/2
            jb=jb+1
            if(jb .gt. nbotx) jb=1
            nbounce(jl,jb)=nbounce(jl,jb) + 2
120      continue
110   continue
c: second, check for any bottom bounces for which ray can't get down
c: to layers where bounces are:
      do 130 jb=nbotx,2,-1
         do 140 jl=1,ndp-1
            if(nbounce(jl,jb) .eq. 0) then
               do 944 jcan=jl,ndp
                  nbounce(jcan,jb-1)=nbounce(jcan,jb-1) + 
     .               nbounce(jcan,jb)
                  nbounce(jcan,jb)=0
944            continue
               goto 130
            endif
140      continue
130   continue
      do 150 jb=1,nbotx
         nbleg(jb)=0
         do 160 jl=1,ndp
            nbleg(jb)=nbleg(jb) + nbounce(jl,jb)
160      continue
150   continue
c
c: compute r vs. z for ocean sections and bottom layers:
      do 10 j=1,3
         rst=0.
         zst=0.
         npt(j)=0
         jjup=jup(j)
         radd=0.
         zadd=0.
         if(ncpth(icpth,j) .eq. 0) goto 10
         do 20 k=kseq(jjup,j,1),kseq(jjup,j,2),kseq(jjup,j,3)
            kk=k+jjup-1
            iilast=0
            if(k .eq. kseq(jjup,j,2)) iilast=1
            call rzlay(kseg,a,ainv,zsvp(kk)-zsvp(kk-1),
     .         csvp(k),csvp(k-iup(j)),g0(kk,jjup),g1(kk,jjup),
     .         g2(kk,jjup),r,dro,dzo,nd,nptmax,kturn,iilast,radd,zadd)
            do 22 jpt=1,nd
               raydr(npt(j)+jpt,j)=rst + dro(jpt)
               raydz(npt(j)+jpt,j)=zst + dzo(jpt)
22          continue
            npt(j)=npt(j) + nd
            rst=raydr(npt(j),j)
            zst=raydz(npt(j),j)
            if(kturn .eq. 1) goto 10
20       continue
         if((iibd .eq. 1) .and. (j .ne. 2)) then
            jj=j-3
            if(ainv .lt. cp1(jj+1)) then
               call delta(1,a,rho2(jj),cp2(jj),cs2(jj),
     .            rho1(jj+1),cp1(jj+1),cs1(jj+1),rb,drb,dr2b,ibd)
               npt(j)=npt(j) + 1
               raydr(npt(j),j)=raydr(npt(j)-1,j) + rb/2.
               raydz(npt(j),j)=raydz(npt(j)-1,j)
c     print *,'bd for j,rb = ',j,rb,ainv,cp2(jj),cp1(jj+1)
            endif
         endif
10    continue
      iilast=0
      do 30 j=1,ndp
         call rzlay(kprof(j),a,ainv,z(j),cp1(j),cp2(j),bp(j),bet(j),0.,
     .      r,drb(1,j),dzb(1,j),nptb(j),nptmax,kturn,iilast,radd,zadd)
cxx      if(kturn .eq. 1 .and. j .lt. ndp) then
         if(kturn .eq. 1) goto 32
30    continue
      if(iibd .eq. 1 .and. ndp .gt. 0) then
         if(ainv .lt. cp1(ndp+1)) then
            call delta(1,a,rho2(ndp),cp2(ndp),cs2(ndp),rho1(ndp+1),
     .         cp1(ndp+1),cs1(ndp+1),rb,drb,dr2b,ibd)
            nptb(ndp)=nptb(ndp) + 1
            drb(nptb(ndp),ndp)=drb(nptb(ndp)-1,ndp) + rb/2.
            dzb(nptb(ndp),ndp)=dzb(nptb(ndp)-1,ndp)
c     print *,'bd at ndp: ndp,rb = ',ndp,rb
         endif
      endif
32    continue
c     print *,'npt = ',npt(1:3)
c     print *,'nptb = ',nptb(1:ndp)
c
      jbot=0
      rst=0.
      zst=zr
c: start plotting at receiver:
c     call plt(zr*zfac,0.,3)
c: temporary write to file while plotting unavailable:
      write(55,200) 0.,ctab,-zr,ctab,2
200   format(f9.4,a1,f8.2,a1,i1)
c: plot the legs of the ray through the ocean sections.
      do 40 ileg=1,nleg
c: plot path in ocean section given by current value of leg:
c        print *,'ileg,leg,iiup = ',ileg,leg,iiup
         if(npt(leg) .eq. 0) goto 54
         if(iiup .eq. iup(leg)) then
            do 50 j=1,npt(leg)
cxx            rpt=rfac*(rst + raydr(j,leg))
cxx            zpt=zfac*(zst - iiup*raydz(j,leg))
c: put range into km:
               rpt=.001*(rst + raydr(j,leg))
               zpt=(zst - iiup*raydz(j,leg))
c              call plt(zpt,rpt,line(jcaus))
c: temporary write to file while plotting unavailable:
               write(55,200) rpt,ctab,-zpt,ctab,0
c     print *,'rpt,zpt = ',j,rpt,zpt
50          continue
            rst=rst + raydr(npt(leg),leg)
            zst=zst - iiup*raydz(npt(leg),leg)
c     print *,'end pt: ',rst,zst
         else
            rst=rst + raydr(npt(leg),leg)
            zst=zst - iiup*raydz(npt(leg),leg)
            do 52 j=npt(leg)-1,1,-1
cxx            rpt=rfac*(rst - raydr(j,leg))
cxx            zpt=zfac*(zst + iiup*raydz(j,leg))
               rpt=.001*(rst - raydr(j,leg))
               zpt=(zst + iiup*raydz(j,leg))
c              call plt(zpt,rpt,line(jcaus))
c: temporary write to file while plotting unavailable:
               write(55,200) rpt,ctab,-zpt,ctab,0
c     print *,'ne rpt,zpt = ',j,rpt,zpt
52          continue
c           call plt(zst*zfac,rst*rfac,line(jcaus))
c: temporary write to file while plotting unavailable:
            write(55,200) .001*rst,ctab,-zst,ctab,0
c     print *,'final: rpt,zpt = ',rst,zst
         endif
54       if((leg .ne. 3) .or. (iiup .ne. -1)) goto 25
c: plot bottom path given in nbounce:
c: (jbot is the current bottom path, jl is the bottom layer being
c: plotted, ibup indicates the direction of the ray)
         jbot=jbot + 1
         jl=0
         ibup=-1
c        print *,'bottom path: jbot,nbleg = ',jbot,nbleg(jbot)
         do 55 jbleg=1,nbleg(jbot)
            jl=jl - ibup
            if(nptb(jl) .eq. 0) goto 62
c     print *,'layer,ibup = ',jl,ibup
            if(ibup .eq. -1) then
               do 57 j=1,nptb(jl)
cxx               rpt=rfac*(rst + drb(j,jl))
cxx               zpt=zfac*(zst - ibup*dzb(j,jl))
                  rpt=.001*(rst + drb(j,jl))
                  zpt=(zst - ibup*dzb(j,jl))
c                 call plt(zpt,rpt,line(jcaus))
c: temporary write to file while plotting unavailable:
                  write(55,200) rpt,ctab,-zpt,ctab,0
c     print *,'rpt,zpt = ',rpt,zpt
57             continue
               rst=rst + drb(nptb(jl),jl)
               zst=zst - ibup*dzb(nptb(jl),jl)
c     print *,'end pt: ',rst,zst,nptb(jl)
            else
               rst=rst + drb(nptb(jl),jl)
               zst=zst - ibup*dzb(nptb(jl),jl)
               do 59 j=nptb(jl)-1,1,-1
cxx               rpt=rfac*(rst - drb(j,jl))
cxx               zpt=zfac*(zst + ibup*dzb(j,jl))
                  rpt=.001*(rst - drb(j,jl))
                  zpt=(zst + ibup*dzb(j,jl))
c                 call plt(zpt,rpt,line(jcaus))
c: temporary write to file while plotting unavailable:
                  write(55,200) rpt,ctab,-zpt,ctab,0
c     print *,'rpt,zpt = ',rpt,zpt
59             continue
c              call plt(zst*zfac,rst*rfac,line(jcaus))
c: temporary write to file while plotting unavailable:
               write(55,200) .001*rst,ctab,-zst,ctab,0
c     print *,'final: ',rst,zst
            endif
62          nbounce(jl,jbot)=nbounce(jl,jbot) - 1
c: check for a reflection rather than a transmission:
            if(((ibup .eq. -1) .and. (nbounce(jl+1,jbot) .eq. 0)) .or.
     .            (ibup .eq. 1) .and. (nbounce(jl-1,jbot) .le. 
     .            nbounce(jl,jbot))) then
               ibup=-1*ibup
               jl=jl + ibup
            endif
55       continue
c: update leg (ocean section) and iiup (direction) for next section:
25       continue
         if(leg .eq. 2) then
            leg=leg - iiup
         elseif(leg+iiup .eq. 2) then
            iiup=-1*iiup
         else
            leg=2
         endif
40    continue
c: temp: return plot cursor horizontally to range 0, then it will go
c: vertically to zr again on next ray:
      write(55,200) 0.,ctab,-zs,ctab,1
c
      return
      end 
