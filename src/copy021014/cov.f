c     234567
      SUBROUTINE covmat(mvec,m1,m2,nf,nmodel,evec,M)
      INTEGER m1,m2,nf,nmodel 
      INTEGER m,i,j,imodel
c      PARAMETER(M=50) 
c     REAL fit(nf)
      INTEGER mvec(m1,m2)   
      REAL mbar(M),total(M),cov(M,M),eval(M),evec(M,M),trace
c     reset
      DO i=1,nf
         DO j=1,nf 
            cov(i,j)=0
         ENDDO
         total(i)=0
      ENDDO
      
      DO imodel=1,nmodel
         DO i=1,nf
            DO j=1,nf
               cov(i,j)=cov(i,j)+mvec(i,imodel)*mvec(j,imodel)
            ENDDO
            total(i)=total(i)+mvec(i,imodel)
         ENDDO
      ENDDO
      
      DO i=1,nf
         mbar(i)=total(i)/nmodel
      ENDDO
      
      DO i=1,nf
         DO j=1,nf
            cov(i,j)=cov(i,j)/nmodel-mbar(i)*mbar(j)
         ENDDO
      ENDDO

c      WRITE(2,*)nf,nmodel
c      DO  j=1,nf
c         WRITE(2,*)'cova',j,(cov(j,i),i=1,nf)
c      ENDDO

      trace=0.0
      DO j=1,nf
         trace=trace+cov(j,j)
      ENDDO

      DO i=1,nf
         DO j=1,nf
            cov(i,j)=cov(i,j)/trace
         ENDDO
      ENDDO

      CALL sa_eig(cov,M,nf,eval,evec)
      END
c     
c     this code finds the eigenvalues of a covariance
c     matrix
c     
      SUBROUTINE sa_eig(cov,M,nf,eval,evec)

c     INCLUDE 'common.h'
      INTEGER M,nf,i,j,nrot
      REAL    enrm
      REAL cov(M,M)
      REAL eval(M),evec(M,M)
c     
      CALL jacobirec(cov,nf,M,eval,evec,nrot)
      CALL eigsrt(eval,evec,nf,M)
c     
      DO 3 i=1,nf
c     WRITE(*,*)'   ev # ',i,' = ',SQRT(eval(i)),eval(i)
         enrm=0.0
         DO 1 j=1,nf
            enrm=amax1(enrm,ABS(evec(j,i)))
 1       CONTINUE
c     
         DO 2 j=1,nf
            evec(j,i)=evec(j,i)/enrm
 2       CONTINUE
 3    CONTINUE
c     
c     DO 5 i=1,nf
c     WRITE(2,*)'lamda',1,(SQRT(eval(i)),i=1,nf)
c     DO 4 j=1,nf
c     WRITE(2,*)'ev',j,(evec(j,i),i=1,nf)
c 4    CONTINUE
c     5 CONTINUE
c     
      RETURN
      END
c     
c     from numerical recipes.
c     
      SUBROUTINE jacobirec(a,n,np,d,v,nrot)
      INTEGER nmax,ip,iq,nrot,n,np,i,j
      PARAMETER (nmax=100)
      REAL a(np,np),d(np),v(np,np),b(nmax),z(nmax)
      REAL tresh,sm,g,h,t,theta,c,s,tau
      DO 12 ip=1,n
         DO 11 iq=1,n
            v(ip,iq)=0.
 11      CONTINUE
         v(ip,ip)=1. 
 12   CONTINUE
      DO 13 ip=1,n
         b(ip)=a(ip,ip)
         d(ip)=b(ip)
         z(ip)=0.
 13   CONTINUE
      nrot=0
      DO 24 i=1,50
         sm=0.
         DO 15 ip=1,n-1
            DO 14 iq=ip+1,n
               sm=sm+ABS(a(ip,iq))
 14         CONTINUE
 15      CONTINUE
         IF(sm.EQ.0.)RETURN
         IF(i.LT.4)THEN
            tresh=0.2*sm/n**2
         ELSE
            tresh=0.
         ENDIF
         DO 22 ip=1,n-1
            DO 21 iq=ip+1,n
               g=100.*ABS(a(ip,iq))
               IF((i.GT.4).AND.(ABS(d(ip))+g.EQ.ABS(d(ip)))
     *              .AND.(ABS(d(iq))+g.EQ.ABS(d(iq))))THEN
                  a(ip,iq)=0.
               ELSE IF(ABS(a(ip,iq)).GT.tresh)THEN
                  h=d(iq)-d(ip)
                  IF(ABS(h)+g.EQ.ABS(h))THEN
                     t=a(ip,iq)/h
                  ELSE
                     theta=0.5*h/a(ip,iq)
                     t=1./(ABS(theta)+SQRT(1.+theta**2))
                     IF(theta.LT.0.)t=-t
                  ENDIF
                  c=1./SQRT(1+t**2)
                  s=t*c
                  tau=s/(1.+c)
                  h=t*a(ip,iq)
                  z(ip)=z(ip)-h
                  z(iq)=z(iq)+h
                  d(ip)=d(ip)-h
                  d(iq)=d(iq)+h
                  a(ip,iq)=0.
                  DO 16 j=1,ip-1
                     g=a(j,ip)
                     h=a(j,iq)
                     a(j,ip)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
 16               CONTINUE
                  DO 17 j=ip+1,iq-1
                     g=a(ip,j)
                     h=a(j,iq)
                     a(ip,j)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
 17               CONTINUE
                  DO 18 j=iq+1,n
                     g=a(ip,j)
                     h=a(iq,j)
                     a(ip,j)=g-s*(h+g*tau)
                     a(iq,j)=h+s*(g-h*tau)
 18               CONTINUE
                  DO 19 j=1,n
                     g=v(j,ip)
                     h=v(j,iq)
                     v(j,ip)=g-s*(h+g*tau)
                     v(j,iq)=h+s*(g-h*tau)
 19               CONTINUE
                  nrot=nrot+1
               ENDIF
 21         CONTINUE
 22      CONTINUE
         DO 23 ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0.
 23      CONTINUE
 24   CONTINUE
      PAUSE '50 iterations should never happen'
      RETURN
      END
c     
c     from numerical recipes.
c     
      SUBROUTINE eigsrt(d,v,n,np)
      INTEGER n,np,i,j,k    
      REAL d(np),v(np,np),p
      DO 13 i=1,n-1
         k=i
         p=d(i)
         DO 11 j=i+1,n
            IF(d(j).GE.p)THEN 
               k=j
               p=d(j)
            ENDIF
 11      CONTINUE
         IF(k.NE.i)THEN
            d(k)=d(i)
            d(i)=p
            DO 12 j=1,n
               p=v(j,i)
               v(j,i)=v(j,k)
               v(j,k)=p
 12         CONTINUE
         ENDIF
 13   CONTINUE
      RETURN
      END
