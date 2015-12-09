      subroutine pasc_make
c
c: this subroutine fills a 76 by 76 pascal's triangle matrix whose
c: upper left element, pasc(1,1), is the upper point of the pascal
c: triangle.
c
      implicit integer*4(i-n)
      include 'common/pascal' 
      npasc=76
      do 10 j=1,npasc
         p_mat(j,1)=1 
10    continue
      do 20 j=2,npasc
         do 30 k=2,j
            p_mat(j,k)=p_mat(j-1,k-1) + p_mat(j-1,k)
30       continue
20    continue
c
      return
      end 
