      subroutine plcalls(range,tyme,plcoh,plinc,xxx,yyy,zzz,jjnr,
     .   ang,blph,blmag)
c
c: this subroutine calls plplot to do the prop loss and transfer
c: function plots.
c
      implicit integer*4(i-n)
      include 'common/timing' 
      include 'common/gamaoptions'
      include 'common/srlox'
      include 'common/freqcom'
      include 'common/depth'
      complex*8 plcoh(nrtot,nfr2)
      real*4 plinc(nrtot,nfr2)
      complex zzz(nrnf)
      real tyme(nrangx),range(nrtot),xxx(nrnf),yyy(nrnf)
      integer*4 jjnr(nrtot)
      character*10 labr,labf,labt,labb
      data labr/'RANGE (KM)'/,labf/'FREQ  (HZ)'/,labt/'SOURCE T/S'/
      data labb/'          '/ 
c
c: plot propagation loss curves: 
      if(iipl .gt. 0) then
         cptim=etime(cpsec)
         do 25 kfr=1,nfr
            do 26 jzr=1,nzr
               zr=xyz(kxy(jzr)+1,3)
               nnxy=nxy(jzr)
               mmr1=mr1(jzr)
               do 27 jxy=1,nnxy
                  mm=mmr1 + jxy - 1
                  do 28 kr=1,nrangx
                     mr=jjnr(mm + nnxy*(kr-1))
                     xxx(kr)=range(mr)
                     zzz(kr)=plcoh(mr,kfr)
                     yyy(kr)=plinc(mr,kfr)
28                continue
                  call plplot(tyme,xxx,nrangx,zzz,yyy,iipl,
     .               labt,10,labr,frq(kfr),labf,frq(kfr),labf,cptim)
27             continue
26          continue
25       continue
      endif
c: output bottom loss table if desired:
      if(iiblt .eq. 1) then
         do 60 kr=1,nrangx
            mr=jjnr(mmr1 + nnxy*(kr-1))
            xxx(kr)=range(mr)
60       continue
         call botlos(range,plcoh,jjnr,ang,blph,blmag)
         close(14)
      endif
c: plot transfer functions: 
      if(iitf .gt. 0) then
         do 29 jzr=1,nzr
            zr=xyz(kxy(jzr)+1,3)
            nnxy=nxy(jzr)
            mmr1=mr1(jzr)
            do 30 kr=1,nrangx 
               do 31 jxy=1,nnxy
                  mr=jjnr(mmr1 + nnxy*(kr-1) + jxy - 1)
                  do 32 kfr=1,nfr
                     zzz(kfr)=plcoh(mr,kfr)
                     yyy(kfr)=plinc(mr,kfr)
32                continue
                  call plplot(frq,frq,nfr,zzz,yyy,iitf,labf,
     .               10,labf,tyme(kr),labt,range(mr),labr,cptim) 
31             continue
30          continue
29       continue
      endif
c
      return
      end 
