      SUBROUTINE plotvelprof
c     reads and interpretates the input file to the genetic algorithm
c     optimization program
c     PETER GERSTOFT, 1992
c
      INCLUDE 'comopt.h'
      include 'comforw.h'
      INTEGER m

      include 'comsnap.h'
      OPEN(unit=8,file='vel',access='append',
     &           status='unknown')
c        WRITE(8,*)' res=['
        do m=1,nd0
          write(8,*)z0(m),c0(m)
        enddo
c        do m=1,nd1
c          write(8,*)z1(m)+z0(nd0),c1(m)
c        enddo
c          write(8,*)z1(nd1)+z0(nd0)+10,c2,'];'
         write(8,*) 999.999, 999.999
        close(8)
      end
  
