      subroutine ckdup
c
c: this subroutine checks the ray paths specified for duplicates.
c
      implicit integer*4(i-n)
      include 'common/paths'
      include 'common/pathchr'
c
c: check for duplicate ocean paths:
      iibad=0
      do 10 k=1,ncp-1
         do 20 j=k+1,ncp
            if((ncpth(j,1) .eq. ncpth(k,1)) .and. (ncpth(j,2) .eq.
     .         ncpth(k,2)) .and. (ncpth(k,3) .eq. ncpth(j,3))) then
      print *,'duplicate ocean path specified: path = ',fdir(j),
     .   '  ',ncpth(j,1)/2,ncpth(j,3)/2
               iibad=1
            endif
20       continue
10    continue
      if(iibad .eq. 1) stop 'Stop in ckdup'
c
      return
      end 
