      subroutine cperr(icpth,fd,ktop,kbot,cz)
c
      implicit integer*4(i-n)
      character*1 fd,cz
c
      print *,'illegal ocean path specification: ',icpth,'  ',fd,
     .      '  ',ktop,'  ',kbot,'  ',cz
      print *,'initial direction must be u or d.  the number of times'
      print *,'the ray enters the top and bottom sections must be '
      print *,'compatible.  the ray type must be a, p, w, or r.'
      print *,'(also, character data must be between single quotes.)' 
c
      stop
      end 
