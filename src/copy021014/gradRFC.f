c****************************************************************************
c     Subroutines for performing gauss-Newton inversion
c     
c     PETER GERSTOFT, 1992 
c****************************************************************************
      SUBROUTINE gaunew(theta,maxiter)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comgrad.h'
C     
C**   local variables
      INTEGER J
      REAL fitn
c-----PARAMETER vector
      REAL theta(mpar)
      REAL thetain(mpar)
      INTEGER maxiter
      COMPLEX a1(mpar*mpar),a2(mpar*mpar)
      INTEGER jdep,ifrq,index,i
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'
c     
      CALL cost(fitn)
      WRITE(*,*)' The initial fitness',fitn       
      WRITE(13,*) 'gradorg=['
      DO jdep=1,ndep 
         DO ifrq=1,nfrq
            j=(ifrq-1)*ndep+jdep
            index=(j-1)*nx
            DO i=1,nx
               WRITE(13,*)REAL(DATA(i+index)),REAL(resp(i+index))
               DATA(i+index)=DATA(i+index)
     1              -resp(i+index)
            ENDDO
         ENDDO
      ENDDO
      WRITE(13,*) '];'

c---  initialize the starting point
c     IF (iopt(4).EQ.2) THEN
      WRITE(*,*)'The starting point  is..'
      nthet=nno
      DO j=1,nthet
         df(j)=0.1
         theta(j)=ccs(j)
         thetain(j)=theta(j) 
         WRITE(*,*)' ',j,theta(j)
      ENDDO
c     ENDIF
C     

      IF (iopt(22).EQ.1) THEN 
         IF (iopt(13).EQ.1) THEN
            IF (nparm.GE.2) THEN
               CALL crrao(thetain,nparm,a1) ! nparm=2 for hess. pl.
            ENDIF

            WRITE(*,*)' calling crline'
            DO j=1,nparm 
               thetain(j)=theta(j)
               WRITE(*,*)' ',j,theta(j)
            ENDDO
            CALL setmodelreal(theta)
            CALL crline(thetain,1,a1,a2) ! nparm=2 for hess. pl.
         ELSE
            WRITE(*,*)' routines found in the directory  old'
c     CALL unceir(nx,ndep,theta,nparm)
c     CALL plhess(nx,ndep,theta,2) ! nparm=2 for hess. pl.
         ENDIF
         RETURN
      ENDIF


      IF (iopt(25).EQ.1) THEN
c     WRITE(*,*) 'gnmin w/o phase'
         CALL gnmin(nx,ndep,theta,maxiter)
      ELSE
         WRITE(*,*) 'gnmin w phase is not developped'
      ENDIF

      IF (iopt(4).EQ.2) THEN
         CALL cost(fitn)
         WRITE(*,*)'The best  energy',fitn          
         DO j=1,nthet
            WRITE(prtfil,*)' ',j,theta(j)
            WRITE(*,*)' ',j,theta(j)
         ENDDO
      ENDIF 

      WRITE(13,*) 'gradmatch=['
      DO jdep=1,ndep 
         DO ifrq=1,nfrq
            j=(ifrq-1)*ndep+jdep
            index=(j-1)*nx
            DO i=1,nx
               WRITE(13,*)REAL(DATA(i+index)),REAL(resp(i+index))
            ENDDO
         ENDDO
      ENDDO
      WRITE(13,*) ']'

      END
c***************************************
      SUBROUTINE jacobi(theta)  
c-    This computes the forward dirivatives of the response.
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comgrad.h'
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'


       REAL theta(Mpar)		 ! the current parameter vector
      INTEGER i,j,index,ix,ifrq,jdep
      REAL   xtemp,xtemp2,xtemp3
c---  set the model and CALL the modelling routine

c Frequecy loop:
        WRITE(*,*)'nthet:',  nthet

       DO 200 ifrq=1,nfrq      
          DO ix=1,nx            
             DO jdep=1,ndep
                j=(ifrq-1)*ndep+jdep
                index=(j-1)*nx
                resp(ix+index)=0.
                DO i=1,nthet
                   xtemp=(xranges(ix)-elenod(i-1))
     1                  *(xranges(ix)-elenod(i+1))
c        WRITE(*,*)'temp1:',  xranges(ix),elemid(i),
c     1                  elelen,xtemp,theta(i)
                   IF (xtemp.LT.0) THEN
                      xtemp3=xranges(ix)-elenod(i)
                      IF (xtemp3.LT.0) THEN  !left  element
                         xtemp2=1-xtemp3/(elenod(i-1)-elenod(i))
                      ELSE                   !right element
                         xtemp2=1-xtemp3/(elenod(i+1)-elenod(i))
                      ENDIF
                      resp(ix+index) =resp(ix+index)+theta(i)*xtemp2
                      grad(ix+index,i)=xtemp2
                   ELSE
                      grad(ix+index,i)=0.
                   ENDIF
                ENDDO
c                WRITE(*,*)' from jacobi: resp', resp(ix+index)
            ENDDO
         ENDDO
 200  CONTINUE
c
      DO i=1,nthet
         jacscale(i)=1
      ENDDO
c
      END



