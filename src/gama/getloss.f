      subroutine getloss(lfile,na,nf,ang,fr,db,ph,fmaxx)
c
c: this subroutine reads in a surface loss or bottom loss table file
c: and places the na angles in array ang(), the nf frequencies in the
c: array fr(), the positive dB losses in db(:,:), and the phases in
c: degrees in ph(:,:)
c
      implicit integer*4(i-n)
      real*4 ang(0:na),fr(0:nf),db(0:nf,0:na),ph(0:nf,0:na)
      character*64 lfile
c
      open(53,file=lfile,status='unknown')
      rewind(53)
      nline=0
      call star(53,nline)
      call star(53,nline)
      read(53,*,end=500,err=500) (fr(jcan),jcan=1,nf-1)
      fr(0)=-1.
      fr(nf)=max(fmaxx,fr(nf-1) + 1.)
      do 5 jf=1,nf
         if((fr(jf) .le. fr(jf-1)) .or. (fr(jf) .lt. 0.)) then
            print *,'non-increasing frequency read from loss table.'
            stop
         endif
5     continue
      ang(0)=-1.
      ang(na)=91.
      call star(53,nline)
      do 10 ja=1,na-1
         read(53,*,end=500,err=500) ang(ja),(db(k,ja),ph(k,ja),k=1,nf-1)
c     print *,'bl read: ',ja,ang(ja),db(1,ja),ph(1,ja)
         if((ang(ja) .le. ang(ja-1)) .or. (ang(ja) .lt. 0.)) then
            print *,'non-increasing angle read from loss table.'
            stop
         endif
         db(0,ja)=db(1,ja)
         db(nf,ja)=db(nf-1,ja)
         ph(0,ja)=ph(1,ja)
         ph(nf,ja)=ph(nf-1,ja)
10    continue
      close(53)
      do 20 jf=0,nf
         db(jf,0)=db(jf,1)
         db(jf,na)=db(jf,na-1)
         ph(jf,0)=ph(jf,1)
         ph(jf,na)=ph(jf,na-1)
20    continue
      do 30 ja=0,na-1
         do 40 jf=1,nf
38          dph=ph(jf,ja) - ph(jf-1,ja)
            if(abs(dph) .gt. 180.) then
               do 940 jcan=jf,nf
                  ph(jcan,ja)=ph(jcan,ja) - sign(360.,dph)
940            continue
c     print *,'adjusting ph for ja,jf = ',ja,jf
               goto 38
            endif
40       continue
         if(ja .eq. na) goto 42
39       dph=ph(0,ja+1) - ph(0,ja)
         if(abs(dph) .gt. 180.) then
            do 942 jcan=0,nf
               do 944 jcan2=ja+1,na
                  ph(jcan,jcan2)=ph(jcan,jcan2) - sign(360.,dph)
944            continue
942         continue
cmay        ph(0:nf,ja+1:na)=ph(0:nf,ja+1:na) - sign(360.,dph)
c     print *,'adjust#2  ph for ja,jf = ',ja,jf
            goto 39
         endif
30    continue
42    continue
c
c     print *,'done getloss. db: '
c     do 60 ja=0,na
c        write(*,100) ja,(db(jcan,ja),jcan=0,nf)
100   format(i1,2x,12(f6.2,1x))
60    continue
c     print *,'ph = '
c     do 70 ja=0,na
c        write(*,110) ja,(ph(jcan,ja),jcan=0,nf)
110   format(i1,2x,12(f6.0,1x))
70    continue
c
      return
500   print *,'error reading loss table file ',lfile
      stop
      end
