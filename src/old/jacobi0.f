      SUBROUTINE jacobi0(theta)  
c*** JACOBI0 DOES NOT COMPUTE GRADIENTS
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comgrad.h'
      REAL theta(Mpar)		 ! the current parameter vector
      INTEGER j,k,index
c---  set the model and call the modelling routine
      REAL*8 rnorm
      CALL setmodelrealjac(theta)
      CALL forw2
c
c---- if we have to use weighting
      IF (iopt(3).ge.3) THEN
        DO k=1,ndep
          index=(k-1)*nx
          DO j=1,nx
           resp(index+j)=resp(index+j)*weight(j)
         ENDDO
        ENDDO
       ENDIF
c
c---  we are only using amplitude information
c
       IF (iopt(25).eq.1) THEN     ! normalizing to real
c--- for response
        DO k=1,ndep
          index=(k-1)*nx
          DO j=1,nx
           resp(index+j)=abs(resp(index+j))
          ENDDO
        ENDDO
       ENDIF                             ! end transformation to real
c
c---- for unit-vector optimization
c
       IF ((iopt(13).eq.1).or.(iopt(27).eq.1)) THEN
c---   for response
        rnorm=0.
        DO k=1,ndep
          index=(k-1)*nx
          DO j=1,nx
           rnorm = rnorm 
     1          +real(resp(index+j))**2 +aimag(resp(index+j))**2
          ENDDO
        ENDDO
c 
        rnorm=1./sqrt(rnorm)
c----   normalize gradient and resp
        DO k=1,ndep
          index=(k-1)*nx
          DO j=1,nx
            resp(index+j)=resp(index+j)*rnorm
         ENDDO
        ENDDO
       ENDIF   ! unit normalization

      END
