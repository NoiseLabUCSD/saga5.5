      SUBROUTINE PLhess(my,ny,thetain,ntheta)
      INCLUDE 'comopt.h'
      INCLUDE 'comgrad.h'
      INCLUDE './oases/complo.f'
c *** Node holograms   
      integer nipt,niptp3,maxpoint  
      real pi
      parameter (nipt=2,pi=3.14159, niptp3=nipt+3,maxpoint=1000)
c-----dimensions of data
      integer my,ny
c-----number of parameters
      integer ntheta
c-----parameter vector
      real theta(mpar),thetain(mpar)
      integer n1,n2,ilog,nlab,nc,ic,j1,i1,k,j,n,m,i
      integer ixlay,iylay,ixpar,iypar,ix3,iy3
      real XLEFT,XRIGHT,XSCALE,XINC,ydown,yup,
     &     dx,dy
      real zmaxx,zmaxy
      real zmin,zmax,zstep,zmax10,zdiv,xmax,ymax,xleft0,ydown0
      real xle,xri,xin,xd,ytop,ylow,yin,yd,ylen,xlen,dxa,dya
      integer nxdif,nydif,info,igrid
c-----pseudo-Hessian
      REAL zx(mx,mpar)   ! matrix to find SVD of 
      REAL ssvd(mpar+1),esvd(mpar),usvd(mx,mpar ),vsvd(mpar,mpar)
      REAL work2(mx)             
      INTEGER index   
c-----
      character*80 titlex,titley,TITLEH
      real svdvector(mpar,mpar,maxpoint)
      real v
      integer nsvd   ! number of observations
      real ca,sa
      real  xval(niptp3),yval(niptp3)
      CHARACTER*6 OPTION(4),OPTTL(2)
      DATA OPTTL /'UHNODE','THNODE'/
      do i=1,ntheta
        theta(i) = thetain(i)
      enddo
      nthet=ntheta
      OPTION(1)=PROGNA
      OPTION(2)='DDDDDD'
      OPTION(3)='      '
      OPTION(4)='      '
c      IF (DEBUG) WRITE(prtfil,*) 'ENTERING PLTLOB'
       ptit='ERROR ANALYSIS'
      NLAB=0
      zmax =-10e20
      zmaxx=-10e20
      zmaxy=-10e20
      nsvd=ny*my
      ixlay= par2lay(1)
      ixpar= par2phy(1) 
      iylay= par2lay(2)
      iypar= par2phy(2)
      if (iopt(12).eq.0) then
        ix3=par3(1)
        iy3=par3(2)
      endif 
       titlex=phystxt(par2phy(1))
       titley=phystxt(par2phy(2))
       xleft  = fmin(1)
       xleft0 = fmin(1)
       xright = fmax(1)
       
       ydown  = fmin(2)
       ydown0 = fmin(2)
       yup    = fmax(2)
      CALL AUTOAX(Xleft,Xright,XLE,XRI,XIN,XD,NXDIF)
      XLE=Xleft
      XRI=Xright
      nxdif=0
      xd=1
c >>> z-axis
c *** write plp file
      CALL AUTOAX(Ydown,Yup,Ylow,Ytop,YIN,YD,NYDIF)
      Ylow = Ydown
      Yup  = Ytop
      nydif=0
      yd=1
       n1=10
       n2=10
c      n1=min(51,ndigit(1))
c      n2=min(51,ndigit(2))
       nc=n2*n1*2 
       dx=(xright-xleft)/(n1+1)    !n1-1
       dy=(yup-ydown)/(n2+1)       !n2-1

c       if (par2phy(1).eq.9) then
c          xLEFT=xLEFT/1000
c          xright=xright/1000
c       endif
c       if (par2phy(2).eq.9) then
c          ydown=ydown/1000
c          yup=yup/1000
c       endif

      ic=0
c--   y-loop
      do j1=1,n2
c---  x-loop 
         do i1=1,n1
c change of x axis          
           theta(1)=(i1)*dx+xleft0        
c change of y axis          
           theta(2)=(j1)*dy+ydown0
           if (iopt(12).eq.0) then
              call setmodelx(ixlay,ixpar,theta(1))
              call setmodelx(iylay,iypar,theta(2))
           else
              call setmodelx(ixlay,ixpar,ix3,theta(1))
              call setmodelx(iylay,iypar,iy3,theta(2))
           endif
c change of x axis
         ic=ic+1          

c-----calculate gradient
      call jacobi(theta)
mkdir c
c-----  computation of SVD
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
    
        call ssvdc(zx,mx,nsvd,nthet,ssvd,esvd,usvd,mx,vsvd,mpar,
     &           work2,21,info)
        if (info.ne.0) THEN
          WRITE(*,*)'********* info from zsvdc ************',info
          stop
        ENDIF
  
      do j=1,ntheta
        if (vsvd(1,j).lt.0) then
          do i=1,ntheta
            vsvd(i,j)=-vsvd(i,j)
          enddo
        endif    
        WRITE(prtfil,'(a,2i4,3e14.6)')
     &            'svd',i1,j1,ssvd(j),(vsvd(i,j),i=1,2)
      enddo  

          if ((1./ssvd(2)).gt.zmax) then
              xmax=theta(1)
              ymax=theta(2)
              zmax = 1./ssvd(2)
            endif

        do i=1,ntheta
          do j=1,ntheta
            svdvector(i,j,ic)= (1./ssvd(j))*vsvd(i,j)     
          enddo !x 
        enddo !x 

          do j=1,ntheta
        if (abs(svdvector(1,j,ic)).gt.zmaxx) then
           zmaxx=abs(svdvector(1,j,ic))
        endif
        if (abs(svdvector(2,j,ic)).gt.zmaxy) then
           zmaxy=abs(svdvector(2,j,ic))
        endif
        enddo

        enddo !x 
      enddo ! y

      write(*,*)' zmax=',zmax
      write(*,*)' zmaxx=',zmaxx
      write(*,*)' zmaxy=',zmaxy
      XTYP='LIN'
      YTYP='LIN'
      IGRID=0
c >>> axes lengths
      xlen=10.0
      ylen=10
      CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,XLEN,YLEN,
     &                  IGRID,XLE,XRI,XIN,XD,titlex,XTYP,
     &                  YLOW,YTOP,YIN,YD,titley,YTYP,NC)
c >>> cos and sin for arrow heads
      ca=cos(PI*0.166666)
      sa=sin(pi*0.166666)
      ic=0
      do 100 i1=1,n1
        do 100 j1=1,n2
            ic=ic+1
            xval(1)=(i1)*dx+xleft0        
c change of y axis          
            yval(1)=(j1)*dy+ydown0
           do 100 k=1,2 
           xval(2)=xval(1)+svdvector(1,k,ic)*dx/zmax
           yval(2)=yval(1)+svdvector(2,k,ic)*dy/zmax
c           xval(2)=xval(1)+svdvector(1,k,ic)*dx/zmax
c           yval(2)=yval(1)+svdvector(2,k,ic)*dy/zmax
c >>> arrow head
       dxa=(xval(nipt-1)-xval(nipt))*ca - 
     &    (yval(nipt-1)-yval(nipt))*sa
       dya=(xval(nipt-1)-xval(nipt))*sa + 
     &    (yval(nipt-1)-yval(nipt))*ca
       xval(nipt+1)=xval(nipt)+dxa
       yval(nipt+1)=yval(nipt)+dya
       xval(nipt+2)=xval(nipt)
       yval(nipt+2)=yval(nipt)
       dxa=(xval(nipt-1)-xval(nipt))*ca + 
     &    (yval(nipt-1)-yval(nipt))*sa
       dya=- (xval(nipt-1)-xval(nipt))*sa + 
     &    (yval(nipt-1)-yval(nipt))*ca
       xval(nipt+3)=xval(nipt)+dxa
       yval(nipt+3)=yval(nipt)+dya
       CALL PLTWRI(nipt,0.,0.,0.,0.,xval(1),1,yval(1),1)
c       CALL PLTWRI(niptp3,0.,0.,0.,0.,xval(1),1,yval(1),1)
 100  continue
c *** formats
      stop                         !      RETURN
      END
