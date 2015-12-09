      SUBROUTINE plotstddev(expfit,totobs)
c     computes standard deviation... genetic algorithm
c     optimization PROGRAM
c     PETER GERSTOFT, 1992
c
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      REAL fval              ! function for computation of real values
      INTEGER i,j,iq
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'
      REAL xsed(10,mq_post),xmeansed(10),svarsed(10),xeoftemp(30)
      REAL*8 xsumsed(10),xsumsqrsed(10),zsedi(10) ! param for statis.
      REAL*8 expfit(mq_post),totobs,xhelp
c     REAL   xbestbest(meofvar),xmeanmean(meofvar)
      INTEGER npoint,jmax
      REAL base_par(4)

      WRITE(*,*)' std dev for mean model...'

      nbase_par=3
      npoint=3
      DO i=1,npoint
         xsumsed(i)=0.
         xsumsqrsed(i)=0.
      ENDDO
 
      OPEN(unit=8,file='velppd.m',access='append',
     &     status='unknown')
      WRITE(8,*)'velsed=['
      
      DO iq=1,q 
         DO i=1,nparm
            CALL seofgen(cap_coef,Ncap, rf.hmsl(1,1),
     1           rf.refmsl(1,1),lvlep_start,i_valid,base_par)
         IF (i_valid.EQ.0) THEN ! this is not a valid profile
            STOP 'a non valid profile'
         ENDIF
         jmax=0
         DO i=1,nbase_par
            j=i
            xsed(i,iq)=base_par(i)
            xhelp=xsed(j,iq)*expfit(iq)
            xsumsed(j)=xsumsed(j)+xhelp
            xsumsqrsed(j)=xsumsqrsed(j)+xhelp*xsed(j,iq)
         ENDDO
c     WRITE(*,*)'jmax,iq,xsed(1,iq)',jmax,iq,xsed(1,iq)
         WRITE(8,'(i3,f8.4,10f8.2)')iq, expfit(iq), 
     .        (xsed(j,iq),j=1,jmax)
         
      ENDDO                     !number of populations, iq
      WRITE(8,*)'];'
      CLOSE(8)
c     c
c---- for standard deviation
c     
      DO i=1,npoint         
         xmeansed(i)=xsumsed(i)/totobs ! /2 
         svarsed(i)
     &        =SQRT(xsumsqrsed(i)/totobs-(xsumsed(i)/totobs)**2)
c     &     =SQRT(xsumsqrsed(i)/2/totobs-(xsumsed(i)/2/totobs)**2)
c     IF (xmean(i).NE.0)svar(i)=svar(i)/ABS(xmean(i))
         WRITE(*,*)'mean,std.dev', xmeansed(i), svarsed(i)
      ENDDO
      zsedi(1)=0
      zsedi(2)=5
      zsedi(3)=10
      zsedi(4)=20
      OPEN(unit=8,file='../vel.m',access='append',
     &     status='unknown')
      
      WRITE(8,*)' vel=['
      DO i=1,npoint
c     WRITE(8,*)zsedi(i),xbestbest(i)-svarsed(i),xbestbest(i),
c     &              xbestbest(i)+svarsed(i)
         WRITE(8,*)zsedi(i),xmeansed(i)-svarsed(i),xmeansed(i),
     &        xmeansed(i)+svarsed(i)
      ENDDO
      WRITE(8,*)' ];'
c     WRITE(8,*) 999.999, 999.999
c     DO i=1,npoint
c     WRITE(8,*)zsedi(i),xmeansed(i)-svarsed(i)
c     ENDDO
c     WRITE(8,*) 999.999, 999.999
      
      END
  
