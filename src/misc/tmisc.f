      program tmisc
      parameter (iargmax = 20)

      integer noa,larg(iargmax)
      character*10 args(iargmax)

      call getcmds(noa,args,larg)
      write (*,*) 'NOA = ',noa
      write (*,*) 'ARGS ARRAY : Length '
      do i= 1,noa
         write (*,*) args(i),larg(i)
      end do

      stop
      end
