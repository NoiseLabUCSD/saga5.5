      subroutine loma(al,ah,a1,a2,kf,nlom)
c
c: this subroutine finds the extremal point which lies
c: between al and ah, for which drl is of opposite sign as drh.
c: the interval is bisected until the a values are close, the r
c: values are close, and the dr values still change sign.
c: a1 and a2 are the si values that incrementally bracket the extremum.
      implicit integer*4(i-n)
      include 'common/pathway'
      include 'common/caustix'
      include 'common/srlox'
      include 'common/gamaoptions'
      real al(3),ah(3),a1(3),a2(3)
      common /info/ienv_info,lu_env,iiwrite
      integer ienv_info,lu_env,iiwrite
c
      nlom=nlom + 1 
      kflip=1
      kbad=0
      kloop=0
      kf=1
      do 940 jcan=1,3
         a1(jcan)=al(jcan)
         a2(jcan)=ah(jcan)
940   continue
      rtoll=.1*rtol 
c     print *,'start of loma: al,ah = ',al,ah
10    continue
c: return if a's get so close that they are equal on the computer: 
      if(a1(1) .eq. a2(2)) then
         print *,'a1=a2 in loma.'
         return
      endif
c: check to see if the extremum is clearly above or below the ranges
c: being looked for.
      if(((a1(3) .gt. 0.) .and. (a1(2) .gt. rmx) .and. (a2(2) .gt. rmx))
     .      .or. ((a1(3) .lt. 0.) .and. (a1(2) .lt. rmn) .and. 
     .      (a2(2) .lt. rmn))) then
         kf=0
         return
      endif
c
      if((abs(a2(1)-a1(1)) .gt. atol).or.(abs(a2(2)-a1(2)) .gt. rtoll))
     .      then
c: alternate between fitting a cubic and bisecting: 
         if((kflip .le. 5) .and. (kbad .eq. 0)) then
            kflip=kflip + 1
            call polfit(a1(1),a2(1),a1(2),a2(2),a1(3),a2(3))
            call cubext(a1(1),a2(1),amid,a2(3),kbad)
         else
            amid=(a1(1) + a2(1))/2. 
         endif
         call rdrcalc(amid,rmid,drmid,ibd,rbd,drbd)
         if(sign(1.,drmid) .eq. sign(1.,a2(3))) then
            a2(1)=amid 
            a2(2)=rmid 
            a2(3)=drmid
         else
            a1(1)=amid 
            a1(2)=rmid 
            a1(3)=drmid
         endif
         kloop=kloop+1
         if(kloop .gt. 100) then
            if(iiwrite.gt.0) then
            print *,'kloop exceeded 100 in loma: ',1./a1(1),1./a2(1),
     .         1./amid,rmid,drmid
            print *,'original: ',1./al(1),al(2),al(3),1./ah(1),ah(2),
     .         ah(3)
            end if
            return
         endif
         goto 10
      endif
c     print *,'loma kloop = ',kloop
c
      if(iidiag .ne. 0) print *,'loma done: 1/a,r,dr,kf,nlom = ',
     .      1./amid,rmid,drmid,kf,nlom
c
      return
      end 
