      SUBROUTINE plotstddev(expfit,totobs)
c     computes standard deviation... genetic algorithm
c     optimization program
c     PETER GERSTOFT, 1992
c
      INCLUDE 'comopt.h'
      include 'comforw.h'
      real fval              ! function for computation of real values
      INTEGER i,j,iq
      include 'comsnap.h'
      real xsed(10,mq_post),xmeansed(10),svarsed(10),xeoftemp(30)
      REAL*8 xsumsed(10),xsumsqrsed(10),zsedi(10) ! param for statis.
      real*8 expfit(mq_post),totobs,xhelp
c      real   xbestbest(meofvar),xmeanmean(meofvar)
      integer npoint,jmax

      write(*,*)' std dev for mean model...'
      npoint=4
      DO i=1,npoint
        xsumsed(i)=0.
        xsumsqrsed(i)=0.
      enddo

 
      OPEN(unit=8,file='velppd.m',access='append',
     &           status='unknown')
      write(8,*)'velsed=['
      
      do iq=1,q 
      DO i=1,nparm
	IF (par2phy(i).EQ.11) THEN
          aeof(par2lay(i))=fval(model(i,iq),i)
	ENDIF
      ENDDO
      CALL eofvalpoint(xeoftemp)
      jmax=0
      DO i=1,neofvar
         j=0
         if (par2phy_eof(i).eq.3) then
           xsed(par2lay_eof(i),iq)=xeoftemp(i)
           j=par2lay_eof(i)
         endif
         if (j.ge.1) then
          xhelp=xsed(j,iq)*expfit(iq)
          xsumsed(j)=xsumsed(j)+xhelp
          xsumsqrsed(j)=xsumsqrsed(j)+xhelp*xsed(j,iq)
         endif
         jmax=max(j,jmax)
      enddo
c      write(*,*)'jmax,iq,xsed(1,iq)',jmax,iq,xsed(1,iq)
      write(8,'(i3,f8.4,10f8.2)')iq, expfit(iq), (xsed(j,iq),j=1,jmax)

      enddo !number of populations, iq
      write(8,*)'];'
      close(8)
cc
c---- for standard deviation
c
      DO i=1,npoint         
        xmeansed(i)=xsumsed(i)/totobs      ! /2 
        svarsed(i)
     &     =sqrt(xsumsqrsed(i)/totobs-(xsumsed(i)/totobs)**2)
c     &     =sqrt(xsumsqrsed(i)/2/totobs-(xsumsed(i)/2/totobs)**2)
c        IF (xmean(i).ne.0)svar(i)=svar(i)/abs(xmean(i))
        write(*,*)'mean,std.dev', xmeansed(i), svarsed(i)
      ENDDO
        zsedi(1)=0
        zsedi(2)=5
        zsedi(3)=10
        zsedi(4)=20
      OPEN(unit=8,file='../vel.m',access='append',
     &           status='unknown')
 
        write(8,*)' vel=['
       do i=1,npoint
c         write(8,*)zsedi(i),xbestbest(i)-svarsed(i),xbestbest(i),
c     &              xbestbest(i)+svarsed(i)
         write(8,*)zsedi(i),xmeansed(i)-svarsed(i),xmeansed(i),
     &              xmeansed(i)+svarsed(i)
       enddo
        write(8,*)' ];'
c       write(8,*) 999.999, 999.999
c       do i=1,npoint
c         write(8,*)zsedi(i),xmeansed(i)-svarsed(i)
c       enddo
c       write(8,*) 999.999, 999.999

      end
  
