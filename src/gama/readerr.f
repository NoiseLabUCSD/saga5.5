      subroutine readerr(nf,nline)
c
      implicit integer*4(i-n)
      character*3 nfile(2)
      data nfile/'opt','svp'/ 
c
      print *,'error reading ',nfile(nf),' file after starred ',
     .   'line # ',nline,'.'
      stop
c
      end 
