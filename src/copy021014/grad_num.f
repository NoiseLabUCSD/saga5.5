      SUBROUTINE gaunew(theta)
c****************************************************************************
c     PROGRAM for inversion of seismo-acoustic DATA using genetic algorithms.
c     The RESULT is presented as the most likely model, the mean model and the 
c     standard variation.
c     
c     PETER GERSTOFT, 1992 
c****************************************************************************
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
C     
C**   local variables
      INTEGER i,J,jn
      REAL ran2
      REAL fitn
c-----PARAMETER vector
      REAL theta(mpar)
      REAL hstep(mpar)
      INTEGER ifirst
      COMMON /deriv/ifirst,hstep
c---  initialize the starting point
      ifirst=0
c     DO j=1,nparm
c     WRITE(*,*)' Random staring point'
c     theta(j)=fmin(j)+ (fmax(j)-fmin(j))*ran2(1)
c     theta(j)=fmin(j)+ (fmax(j)-fmin(j))*0.5
c     ENDDO


c     WRITE(*,*) 'starting'
      
      theta(7)=45
      theta(8)=45
      theta(1)=90               !1650
      theta(2)=1750
      theta(3)=2500
      theta(4)=550
      theta(5)=1050
      theta(6)=1500


c     theta(1)=1600
c     theta(2)=1600       
c     theta(3)=1600       
c     theta(4)=500       
c     theta(5)=500       
c     theta(6)=500       
      WRITE(*,*)'The starting point  is..'
      DO j=1,nparm
         WRITE(*,*)' ',j,theta(j)
      ENDDO
      WRITE(prtfil,*)' The starting point  is..'
      DO j=1,nparm
         WRITE(prtfil,*)' ',j,theta(j)
      ENDDO
      WRITE(*,*)'nx,ny',nx,ndep
c     
      IF (iopt(22).EQ.1) THEN 
c     CALL uncer(nx,ndep,theta,nparm)
c     CALL plhess(nx,ndep,nparm)
      ENDIF

      CALL gnmin(nx,ndep,theta,nparm)

      CALL setmodelreal(theta)
      CALL forw2
      IF (iopt(3).GE.1) CALL norm ! normalize the response 
      CALL cost(fitn)

      WRITE(*,*)'The best  energy',fitn
      WRITE(prtfil,*)'The best  energy',fitn

      DO j=1,nparm
         WRITE(prtfil,*)' ',j,theta(j)
         WRITE(*,*)' ',j,theta(j)
      ENDDO

      END
c**************************************************************
      SUBROUTINE jacobi(grad,theta)  
c     -    This computes the forward dirivatives of the response.
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
c     
      REAL theta(Mpar)          ! the current parameter vector
      COMPLEX respx(Mx,mdep)    ! the computed data points
      COMPLEX grad(Mx,mdep,Mpar) ! the gradient 
      INTEGER ixlay,ixpar,i,j,k
      REAL thetapar,h,eps,xsum,fitn
c     DATA eps/3.e-1/  ! trial for range
c     DATA eps/1.e-3/  ! trial for range
c     DATA eps/1.e-4/  ! thial for velocity (sspmis)
c     DATA eps/8.e-5/  ! thial for velocity
      DATA eps/1.5e-3/          ! thial for velocity
c     DATA eps/3.e-3/  ! thial for velocity (collins.dat)
c     DATA eps/3.e-3/  ! this should correspond to sqrt(maceps)
c---  set the model and CALL the modelling routine
      REAL hstep(mpar)
      INTEGER ifirst
      COMMON /deriv/ifirst,hstep
c     WRITE(prtfil,*)'*************************************************'
c     WRITE(prtfil,*)'theta',(theta(i),i=1,nparm)
      DO i=1,ndep
         DO j=1,nx
            respx(j,i)=ABS(resp(j,i)) ! move the computed response
         ENDDO
      ENDDO

      IF (ifirst.NE.1) THEN
         ifirst=1
c---- find step value
         DO i=1,nparm
            hstep(i)=theta(i)*eps/100 !*100
c     CALL dfridr(respx,i,theta(i),hstep(i))       
            hstep(i)=2.5
c*******************************
            WRITE(*,*)' hstep',i,hstep(i)
         ENDDO

      ENDIF
c---  now compute the  jacobian
      DO i=1,nparm
         xsum=0
         ixlay= par2lay(i)
         ixpar= par2phy(i)
c---- the stepsize is chosen as SQRT(maxeps)*MAX(xc,xref)
c     h=eps*MAX(theta(i),1.)
c     IF (h.EQ.0) h=eps
         h=hstep(i)
         thetapar=theta(i)+h
         CALL setmodelx(ixlay,ixpar,thetapar)
         CALL forw2
         IF (iopt(3).GE.1) CALL norm ! normalize the response 
         CALL cost(fitn)
         DO k=1,ndep
            DO j=1,nx
c******this is important *************
c---  note gradient is normalized
c     grad(j,k,i)=(resp(j,k)-respx(j,k))/h*theta(i)
               grad(j,k,i)=(ABS(resp(j,k))-ABS(respx(j,k)))/h*theta(i)
c     grad(j,k,i)=(ABS(resp(j,k))-ABS(respx(j,k)))/h
               xsum=xsum+ABS(grad(j,k,i)) 
            ENDDO
c     WRITE(*,*)j,k,grad(1,k,i),ABS(resp(1,k)),ABS(respx(1,k))
         ENDDO
         WRITE(prtfil,*)'sum of deriv',i,xsum
c     WRITE(*,*)'sum of deriv',i,xsum
         
         thetapar=theta(i)
         CALL setmodelx(ixlay,ixpar,thetapar)

      ENDDO        
c     WRITE(*,*) 'the gradient are',(grad(1,1,i),i=1,nparm)

      END


c**************************************************************
      SUBROUTINE dfridr(respx,iparm,theta,h)  
c     -    This dirivatives of one PARAMETER. this SUBROUTINE is based on 
c     -    numerical recipiers p 183, Peter gerstoft 9307
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'

      REAL theta                ! the current parameter vector
      COMPLEX respx(Mx,mdep)    ! the computed data points
      INTEGER ixlay,ixpar,i,j,k,imax,iparm
      REAL thetapar,h,eps,xsum,fitn
      REAL con,con2,big,safe,fac,hh
      INTEGER ntab
      PARAMETER (con=4, con2=con*con,big=1e30,ntab=50,safe=2)
      REAL a(ntab,ntab),errt,err,xdiff
c---  set the model and CALL the modelling routine
c     WRITE(prtfil,*)'*************************************************'
c     WRITE(prtfil,*)'theta',theta
c---  now compute the  jacobian
      ixlay= par2lay(iparm)
      ixpar= par2phy(iparm)
c---- the stepsize is chosen as SQRT(maxeps)*MAX(xc,xref)
c     h=eps*MAX(theta(i),1.)
c     IF (h.EQ.0) h=eps
      hh=h*con
      err=big
      DO 999 i=1,ntab
         hh=hh/con
         thetapar=theta+hh
         CALL setmodelx(ixlay,ixpar,thetapar)
         CALL forw2
         IF (iopt(3).GE.1) CALL norm ! normalize the response 
         CALL cost(fitn)
         xsum=0
         DO k=1,ndep
            DO j=1,nx
c---  note gradient is normalized
               xsum=xsum+ABS((ABS(resp(j,k))-ABS(respx(j,k)))/hh*theta)
            ENDDO
         ENDDO
         a(1,i)=xsum
         fac=con2
         thetapar=theta
         CALL setmodelx(ixlay,ixpar,thetapar) ! return to old vqalue
         IF (i.EQ.1) GOTO 999
c---  
         DO j=2,i,1
c     WRITE(*,*)'i,j',i,j
            a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.)
c     WRITE(*,*)j,i,a(j,i),  a(j-1,i),fac,a(j-1,i-1),(fac-1.)
            fac=fac*con2
            errt=MAX(ABS(a(j,i)-a(j-1,i)),ABS(a(j,i)-a(j-1,i-1)))
            IF (errt.LE.err) THEN
               err=errt
               xdiff=a(j,i)
            ENDIF
         ENDDO	  
         IF (ABS(a(i,i)-a(i-1,i-1)).GE.safe*err) THEN
            GOTO 3000   
         ENDIF
 999  CONTINUE
 3000 imax=i
      h=hh*con
      WRITE(*,*)'best derivative obtained',imax,xdiff
      WRITE(prtfil,*)'derivatives'
      DO i=1,imax
         WRITE(prtfil,'(10e11.4)')(a(j,i),j=1,i)
      ENDDO
      END


