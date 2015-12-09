      subroutine star(nf,nline)
c
c: this subroutine searches for the next line in an input file that
c: begins with key symbol '*', which indicates that the next line
c: has data that is to be read.
c
      implicit integer*4(i-n)
      character*3 nfile(3)
      character*1 chfs,chstar 
      data chstar/'*'/,nfile/'opt','svp','tab'/
c
      nline=nline + 1
10    read(nf,100,end=500,err=500) chfs 
100   format(a1)
      if(chfs .ne. chstar) goto 10
      return
500   print *,'program could not find starred line # ',nline,
     .   ' from file ',nfile(nf)
      stop       
c
      end 
