      SUBROUTINE gnmin(my,ny,theta,maxiter)
c-----DESCRIPTION
c     Least square minimization using the damped Gauss-Newton algorithm.
c     by using SVD 
c     my=nx, ny=ndep

      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comgrad.h'
      INCLUDE 'comforw.h'
c-----dimensions of DATA
      INTEGER my,ny
c-----PARAMETER vector
      REAL theta(mpar)
      INTEGER maxiter           !    maximum number of iterations     
c     REAL tol/0.00001/           !    minimum norm of step-vector
      INTEGER maxbisect      !    maximum number of bisections
      data maxbisect/3/
c-----LOCAL PARAMETERS
c-----gradient vector
      REAL G(mPAR)
c-----pseudo-Hessian
      REAL zx(mx,mpar)          ! matrix to find SVD of 
c     EQUIVALENCE (grad(1,1),zx(1,1))
      REAL ssvd(mpar+1),esvd(mpar),usvd(mx,mpar ),vsvd(mpar,mpar)
      REAL work2(mx),wsvd(mpar)             
      INTEGER nsvd              ! number of observations
c-----old PARAMETER vector
      REAL theta0(mpar)
      REAL V,V0                 ! current and old value of square error
c-----counters
      INTEGER i,j,k,n,m,index
      LOGICAL stop2             ! loop control
      LOGICAL linfo             ! info control
      INTEGER iterloc,iterforw,ieps !  number of iterations
c-----number of bisections
      INTEGER bisect, ibisect_loop
c-----information returned
      INTEGER info
c-----norm of step-vector
c     REAL normparm
c-----step size
c     REAL mu,scaleeps
      REAL eps
c-----check dimensions
      nsvd=ny*my
      IF(nsvd.GT.mx) STOP 'gnmin: Nsvd > mx; increase mx'
      IF (nthet.GT.mpar) STOP ' gnmin: increase mpar  '
c-----calculate gradient
      CALL jacobi(theta)
      iterforw=1
      ibisect_loop=0
      ieps=0
      MAXBISECT=3
      IF (iopt(4).NE.2) THEN
         linfo=.FALSE.
      ELSE
         linfo=.TRUE.
      ENDIF

      V=0
      DO n=1,ny
         index=(n-1)*my
         DO m=1,my
            resp(m+index)=(resp(m+index))-(DATA(m+index))
            V=V+(resp(m+index))**2
            WRITE(*,*)n,m,resp(m+index),(resp(m+index)),
     &           (DATA(m+index)),V
         ENDDO 
      ENDDO 
      IF (linfo) WRITE(*,*)'the start fitness is',V

      eps=1e-20
      DO 2000 iterloc=1,maxiter
         IF (linfo) WRITE(prtfil,*)' iteration',iterloc 
         IF (linfo) WRITE(*,*)' iteration',iterloc 
         stop2=.FALSE.          ! required accuracy
c     
c-----computation of SVD
c     
         DO j=1,nthet
            k=0
            DO n=1,ny
               index=(n-1)*my
               DO m=1,my
                  k=k+1
                  zx(k,j)=grad(m+index,j)             
               ENDDO 
            ENDDO 
         ENDDO 
         
         CALL ssvdc(zx,mx,nsvd,nthet,ssvd,esvd,usvd,mx,vsvd,mpar,
     &        work2,21,info)
         IF (info.NE.0) THEN
            WRITE(*,*)'********* info from zsvdc ************',info
            STOP
         ENDIF
c     
c---  WRITE out the singular vectors
c     
         IF (linfo) THEN 
            DO j=1,nthet
               WRITE(prtfil,*)'svd',ssvd(j)
            ENDDO  
            WRITE(prtfil,*)' v-svd'
            DO i=1,nthet
               WRITE(prtfil,'(12f8.4)')(vsvd(i,j),j=1,nthet)
            ENDDO  
            WRITE(prtfil,*)' u-svd'
            DO i=1,nthet
               WRITE(prtfil,'(12f9.5)')(usvd(i,j),j=1,nthet)
            ENDDO  
         ENDIF
c     
c     
c     
         DO j=1,nthet
            wsvd(j)=(0.,0.)
            k=0
            DO n=1,ny
               index=(n-1)*my
               DO m=1,my
                  k=k+1
                  wsvd(j)= wsvd(j)+(usvd(k,j))*resp(m+index)
               ENDDO 
            ENDDO 
         ENDDO  

         IF (ieps.EQ.0) THEN
            eps=ABS(ssvd(1))**2 *0.00001
         ENDIF
         ibisect_loop=0
 100     CONTINUE
c---  reset 
         DO j=1,nthet
            g(j)=(0.,0.) 
         ENDDO  
c     
c---  NEW step size
c     
         DO i=1,nthet
            DO j=1,nthet
               g(j)=g(j)+wsvd(i)*(ssvd(i)/(ssvd(i)**2+eps))*vsvd(j,i) 
            ENDDO  
         ENDDO  

c-----calculate norm
c     normparm=0
c     DO  i=1,nthet
c     normparm=normparm+G(i)**2
c     ENDDO
c     normparm=SQRT(normparm)
c     IF (linfo)WRITE(*,*) ' norm of parameters', normparm
c-------Find the jump in physical numbers.
         DO i=1,nthet
            g(i)=g(i)*jacscale(i)
         ENDDO
c-----SAVE old values
 1099    DO  i=1,nthet
            theta0(i)=theta(i)
         ENDDO
         V0=V
c-----initialization of loop-values
         bisect=0
 1100    CONTINUE

         IF (.NOT.stop2) THEN
c-------calculate NEW guess (consistency check)
            DO 1140 i=1,nthet
               theta(i)=(theta0(i)-G(i))
 1140       CONTINUE
c-------perform cosistency check on the NEW guess
            DO  i=1,nthet
c     temporary disabled            IF (theta(i).LT.fmin(forw2opt(i)) ) THEN
c     theta(i)=fmin(forw2opt(i))
c     ELSEIF (theta(i).GT.fmax(forw2opt(i)) ) THEN
c     theta(i)=fmax(forw2opt(i))
c     ENDIF
            ENDDO
c-------calculate DATA and gradient
            iterforw=iterforw+1
            IF (iterforw.GE.maxiter) THEN
c     CALL jacobi0(theta) w/o gradient disabled
               CALL jacobi(theta)
            ELSE
               CALL jacobi(theta)
            ENDIF
c--------compute error
            V=0
            DO n=1,ny
               index=(n-1)*my
               DO m=1,my
                  resp(m+index)=(resp(m+index))-(DATA(m+index))
                  V=V+(resp(m+index))**2
               ENDDO 
            ENDDO 

            IF (linfo) WRITE(*,*)'new, old error**2',V,V0
c     IF (normparm.LT.tol) THEN
c     stop2=.TRUE.
c     IF (linfo) WRITE(*,*)'stopping due to norm of parameters'
c     GOTO 111
c     ENDIF
            DO i=1,nthet
               IF (ABS(g(i)).GT.df(i)/5) THEN
                  GOTO 111
               ENDIF
            ENDDO
            IF (linfo) WRITE(*,*)' Stopping, crit. for each parameter'
            stop2=.TRUE.
 111        CONTINUE

            IF (V.LE.V0) THEN
c---------accept NEW guess
               eps=eps/2
               IF(linfo)WRITE(*,*)' new theta',iterloc,
     &              (theta(i),i=1,nthet)
               IF (iterforw.GE.maxiter) GOTO 2222
            ELSE
               V=V0
c---------reestablish old values
               DO 1130 i=1,nthet
                  theta(i)=theta0(i)
 1130          CONTINUE
               eps=eps*10
               IF (ieps.EQ.0) eps=eps*1000 
               ieps=1
               IF (linfo)WRITE(*,*)'new eps',eps
               bisect=bisect+1
               IF (iterforw.GE.maxiter) GOTO 2222
               IF (stop2) THEN
               ELSEIF (bisect.GT.maxbisect) THEN
                  stop2=.TRUE.
                  IF (linfo) WRITE(*,*)' Stopping due to max bisections'
               ELSE
                  ibisect_loop=1
                  GOTO 100
               ENDIF
            ENDIF
         ENDIF
         IF (stop2) GOTO 2222
 2000 CONTINUE
      IF (linfo) WRITE(*,*)' Stopping due to maxiterations'
 2222 CONTINUE
      IF (linfo) WRITE(*,*)' no of forward runs',iterforw
      
c---  For matlab
      IF (iopt(4).EQ.2) THEN
c---- computation of covariance
c---- computation of resolution
c---- computation of SVD
         CALL jacobi(theta)

         DO j=1,nthet
            k=0
            DO n=1,ny
               index=(n-1)*my
               DO m=1,my
                  k=k+1
                  zx(k,j)=grad(m+index,j)             
               ENDDO 
            ENDDO 
         ENDDO 
         OPEN(13,status='unknown')
         WRITE(13,*)' T= ['  
         DO m=1,nsvd
            WRITE(13,'(20f8.3)')(zx(m,j),j=1,nthet)
         ENDDO
         WRITE(13,*)' ];'
         WRITE(13,*)' layers=',nthet
         WRITE(13,*)' thiclness=-(1:1:layers)'';'
         WRITE(13,*)' lambda=',eps
         WRITE(13,*)' h=eye(layers,layers);'
         WRITE(13,*)' sol=['
         WRITE(13,'(1(15f8.2))')(theta(i),i=1,nthet)
         WRITE(13,*)']'''
c     CLOSE(13)
      ENDIF

      END







