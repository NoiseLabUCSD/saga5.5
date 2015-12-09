
      SUBROUTINE crrao(theta,ntheta,crb)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comgrad.h'
      INCLUDE './oases/complo.f'
c-----number of parameters
      INTEGER ntheta
      COMPLEX a1(100000),a2(100000),a3(100000),a4(100000)
      COMPLEX a5(100000),a6(100000),a7(100000),a8(100000)
      
      COMPLEX crb(ntheta,ntheta)
c     *** Node holograms   
      INTEGER nipt,maxpoint  
      PARAMETER (nipt=6,maxpoint=10000)
c-----PARAMETER vector
      REAL theta(mpar),t1
      INTEGER n1,n2,nlab,nc,ic,j1,i1,j,i,idep 
      INTEGER ixlay,iylay,ixpar,iypar,ix3,iy3
      REAL XLEFT,XRIGHT,ydown,yup,dx,dy
      REAL zmaxx,zmaxy,scale
      REAL zmax,xmax,ymax
      REAL xle,xri,xin,xd,ytop,ylow,yin,yd
      REAL ylen,xlen    ! lenght of axis
      data zmaxx/-10e20/,zmaxy/-10e20/, zmax/-10e20/,ylen/10/,xlen/10/ 
      INTEGER nxdif,nydif,info,igrid
      data nlab/0/,igrid/0/
      INTEGER jgrad,igrad       !,ix,j1start,ibart,ierr
c-----for svd analysis of a quadratic matrix
      REAL hess(mpar,mpar)      ! matrix to find SVD of 
      REAL ssvd(mpar+1),esvd(mpar),usvd(mpar,mpar ),vsvd(mpar,mpar)
      REAL work2(mpar)
c-----------
      REAL f_crb(210,210,10)
      CHARACTER*80 titlex,titley
      REAL svdvector(mpar,mpar,maxpoint)
      COMPLEX cx            
      REAL  xval(nipt),yval(nipt)
      CHARACTER*6 OPTION(6)
      data  OPTION/'      ','DDDDDD',4*'      '/

c---- For contour plot:
      CHARACTER*80 conopt 
      REAL XSCALE,XINC,ySCALE,yINC,
     &     ypointup,ypointdown,xpointleft,xpointright,
     &     zmin,zstep,v
      conopt='COL,NCL,SCS                '
      OPTION(1)=PROGNA
      ptit='ERROR ANALYSIS'
      ixlay= par2lay(1)
      ixpar= par2phy(1) 
      iylay= par2lay(2)
      iypar= par2phy(2)
      IF (iopt(12).EQ.0) THEN
         ix3=par3(1)
         iy3=par3(2)
      ENDIF 
      titlex=phystxt(par2phy(1))
      titley=phystxt(par2phy(2))
      xleft  = fmin(1)
      xright = fmax(1)
      
      ydown  = fmin(2)
      yup    = fmax(2)
      CALL AUTOAX(Xleft,Xright,XLE,XRI,XIN,XD,NXDIF)
      XLE=Xleft
      XRI=Xright
      nxdif=0
      xd=1
c     >>> z-axis
c     *** WRITE plp file
      CALL AUTOAX(Ydown,Yup,Ylow,Ytop,YIN,YD,NYDIF)
      Ylow = Ydown
      Ytop = Yup
      nydif=0
      yd=1

      n1=MIN(201,ndigit(1))
      n2=MIN(201,ndigit(2))
      nc=n2*n1
      IF (nc.GT.maxpoint)   STOP 'dimension of maxpoint'
      dx=(xright-xleft)/(n1+1)  
      dy=(yup   -ydown)/(n2+1)  
c     
c---- CRB for reference environment
c     
c-----calculate gradient
      CALL jacobi(theta)
c     
c-----computation of simple fisher
      CALL  crb_sto(ntheta,ndep,crb,a1,a2,a3,a4,a5,a6,a7,a8)
c---- WRITE to matlab file
      WRITE(13,*)' crb= ['  
      DO i=1,ntheta
         WRITE(13,'(20g13.3)')(REAL(crb(i,j)),j=1,ntheta)
      ENDDO
      WRITE(13,*)' ];'  
      WRITE(13,'(a,20f10.3)')'theta=[',(theta(j),j=1,ntheta)
      WRITE(13,*)' ];'  
      

      ic=0
c--   y-loop
      DO j1=1,n2
         IF (j1.EQ.2) THEN
            CALL rdtime(t1)
            WRITE(*,320) t1*n2
 320        FORMAT(' >>> Estimated contour plot time:',F12.3,' secs.')
         ENDIF
c---  x-loop 
         DO i1=1,n1
c     change of x & y  axis          
            theta(1)=(i1)*dx+xleft       
            theta(2)=(j1)*dy+ydown
            IF (iopt(12).EQ.0) THEN
               CALL setmodelx(ixlay,ixpar,theta(1))
               CALL setmodelx(iylay,iypar,theta(2))
            ELSE
               CALL setmodelx(ixlay,ixpar,ix3,theta(1))
               CALL setmodelx(iylay,iypar,iy3,theta(2))
            ENDIF
            ic=ic+1          

c-----calculate gradient
            CALL jacobi(theta)
c     
c-----computation of simple fisher
c     
            CALL  crb_sto(ntheta,ndep,crb,a1,a2,a3,a4,a5,a6,a7,a8)

            DO jgrad=1,ntheta
               DO igrad=1,ntheta
                  cx=0
                  DO idep=1,ndep
                     cx=cx +CONJG(grad(idep,igrad))*grad(idep,jgrad)
                  ENDDO
                  hess(jgrad,igrad)=REAL(cx)
                  hess(jgrad,igrad)=REAL(crb(jgrad,igrad))
               ENDDO 
               f_crb(i1,n2+1-j1,jgrad+1)=SQRT(REAL(crb(jgrad,jgrad)))
            ENDDO 

            CALL norm
            CALL cost(v)

c     WRITE(*,*)i,j,v
            IF (iopt(13).EQ.1 .AND. iopt(20).EQ.0) THEN
               IF (v.GE.0) THEN
                  v=10.*LOG10(  ABS(v)  )
               ELSE
                  v=-300
               ENDIF
            ELSE
               IF (v.NE.1) THEN
                  v=10.*LOG10(  ABS(1.-v)  )
               ELSE
                  v=-300
               ENDIF
            ENDIF
            f_crb(i1,n2+1-j1,1)=v       

c     
            CALL ssvdc(hess,mpar,2,2,ssvd,esvd,usvd,mpar,vsvd,mpar,
     &           work2,01,info)
            IF (info.NE.0) THEN
               WRITE(*,*)'********* info from zsvdc ************',info
               STOP
            ENDIF

            IF (ABS(ssvd(1)).EQ.0) THEN
               WRITE(*,*)'First eigenvalue is zero!, j,i1,theta='
               WRITE(*,*)1,i1,theta(1)
               STOP
            ENDIF
            DO j=1,ntheta
c     ssvd(j)=1./SQRT(ssvd(j))
               IF ((vsvd(1,j)).LT.0) THEN
                  DO i=1,ntheta
                     vsvd(i,j)=-vsvd(i,j)
                  ENDDO
               ENDIF    
               WRITE(prtfil,'(a,2i4,6e14.6)')
     &              'svd',i1,j1,ssvd(j),(vsvd(i,j),i=1,2)
            ENDDO  
            
            IF ((ABS(ssvd(1))).GT.zmax) THEN
               xmax=theta(1)
               ymax=theta(2)
               zmax = ABS(ssvd(1))
            ENDIF
            
            DO i=1,ntheta
               DO j=1,ntheta
                  svdvector(i,j,ic)= (ssvd(j))*vsvd(i,j)     
c     svdvector(i,j,ic)= vsvd(i,j)     
               ENDDO              
            ENDDO                
c     
            DO j=1,ntheta
               IF (ABS(svdvector(1,j,ic)).GT.zmaxx) THEN
                  zmaxx=ABS(svdvector(1,j,ic))
               ENDIF
               IF (ABS(svdvector(2,j,ic)).GT.zmaxy) THEN
                  zmaxy=ABS(svdvector(2,j,ic))
               ENDIF
            ENDDO
         ENDDO                  !x 
      ENDDO                     ! y

c     WRITE(*,*)' zmax=',zmax
      WRITE(*,*)' zmaxx=',zmaxx
      WRITE(*,*)' zmaxy=',zmaxy
      XTYP='LIN'
      YTYP='LIN'
      CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,XLEN,YLEN,
     &     IGRID,XLE ,XRI ,XIN,XD,titlex,XTYP,
     &     YLOW,YTOP,YIN,YD,titley,YTYP,NC)

      ic=0
      scale=MIN(dx/zmaxx,dy/zmaxy)*0.45
      WRITE(*,*)'scale=',scale
      DO 100 j1=1,n2
         DO 100 i1=1,n1
            ic=ic+1
            xval(2)=(i1)*dx+xleft        
            yval(2)=(j1)*dy+ydown
            xval(1)=xval(2)+svdvector(1,1,ic)*scale
            yval(1)=yval(2)+svdvector(2,1,ic)*scale
            xval(3)=xval(2)+svdvector(1,2,ic)*scale
            yval(3)=yval(2)+svdvector(2,2,ic)*scale
            xval(4)=xval(2)-svdvector(1,2,ic)*scale
            yval(4)=yval(2)-svdvector(2,2,ic)*scale
            xval(5)=xval(2)
            yval(5)=yval(2)
            xval(6)=xval(2)-svdvector(1,1,ic)*scale
            yval(6)=yval(2)-svdvector(2,1,ic)*scale
            CALL PLTWRI(nipt,0.,0.,0.,0.,xval(1),1,yval(1),1)
 100     CONTINUE

c----------Create a contour plot
         DO 200 igrad=1,ntheta+1
            v=0
            zmin=10e20
            zmax=-10e20
            DO j=1,n2
               DO i=1,n1
                  v= f_crb(i,j,igrad)
                  IF (v.GT.zmax) THEN
                     zmax = v
                  ENDIF
                  IF (v.LT.zmin) THEN
                     zmin = v
                  ENDIF
               ENDDO
            ENDDO
            yinc=yinc/2
            ypointup   = ytop-dy
            ypointdown  =ylow+dy
            xpointleft =xle+dx
            xpointright=xri-dx
            XSCALE= ABS(XRIGHT-XLEFT)/11. ! x - axis are 11cm.
            YSCALE= ABS(yup-ydown)/11. ! Y- AXIS IS 11 CM.
c     zmin=0
            zstep=(zmax-zmin)/10
            
            CALL CONDRW1(TITLE,N1,n2,n1,n2,XLE,XRI,XSCALE,XINC,
     1           ylow,ytop,YSCALE,yINC,ZMIN,ZMAX,
     2           ZSTEP,5.,1.,ypointup,ypointdown,xpointleft,
     &           xpointright,1,titlex,titley,conopt)
c--   y-loop
            DO j=1,n2
               DO i=1,n1
c     fobj(i,j)=fobj(i,j)*zdiv ! *zmax10/zmax
               ENDDO
               CALL CONDRB(1,1,n1,f_crb(1,j,igrad))
            ENDDO
 200     CONTINUE
         

         END

c***************************************************************
      SUBROUTINE crline(theta,ntheta,crb,crb2)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comgrad.h'
      INCLUDE './oases/complo.f'
c-----PARAMETER vector
      REAL theta(mpar)
      INTEGER ntheta
      INTEGER n1,ic,i1
      INTEGER ixlay,ixpar,ix3
      REAL XLEFT,XRIGHT,dx
c-----
      COMPLEX a1(100000),a2(100000),a3(100000),a4(100000)
      COMPLEX a5(100000),a6(100000),a7(100000),a8(100000)
      
      COMPLEX crb(ntheta,ntheta),crb2(ntheta,ntheta)

      nthet=ntheta
      ixlay= par2lay(1)
      ixpar= par2phy(1) 
      IF (iopt(12).EQ.0) THEN
         ix3=par3(1)
      ENDIF 
c     titlex=phystxt(par2phy(1))
c     titley=phystxt(par2phy(2))
      xleft  = fmin(1)
      xright = fmax(1)

c     n1=50
      n1=MIN(51,ndigit(1))

      dx=(xright-xleft)/(n1+1)  !n1-1

      ic=0
cf90: position='append', f77: access='append',
      OPEN(unit=80,file='crb.m',position='append',status='unknown')
      WRITE(80,*)'crb =[' 
c     OPEN(unit=81,file='grad.m',access='append',status='unknown')
c     WRITE(81,*)'r =[' 
c---  x-loop 
      DO i1=1,n1
c     change of y axis          
         theta(1)=(i1)*dx+xleft       
         IF (iopt(12).EQ.0) THEN
            CALL setmodelx(ixlay,ixpar,theta(1))
         ELSE
            CALL setmodelx(ixlay,ixpar,ix3,theta(1))
         ENDIF
         ic=ic+1          
c-----calculate gradient
         CALL jacobi(theta)
c     
c     ----stochastic CRB
         CALL  crb_sto(ntheta,ndep,crb,a1,a2,a3,a4,a5,a6,a7,a8)

c     
c-----deterministic
c     

         CALL  crb_det(ntheta,crb2,a1)
         WRITE(80,*)theta(1),SQRT(REAL(crb2(1,1))),
     &        SQRT(REAL(crb(1,1)))
c     DO ifrq=1,nfrq
c     i=2+(ifrq-1)*ndep
c     WRITE(81,*)
c     &           theta(1),ifrq,REAL(resp(i)),imag(resp(i)), 
c     &           REAL(grad(i,1)),imag(grad(i,1)) 
c     ENDDO 
      ENDDO                     !x 
c     *** formats
      WRITE(80,*)'];' 
c     WRITE(81,*)'];' 
      END

c***************************************************
      SUBROUTINE crb_sto(ntheta,ndep1,crb,cxmat,cm_inv,ctemp,ctot,
     &     cx_diff,cfish)
c     1       2     3    4    5      6
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comgrad.h'
      INCLUDE './oases/complo.f'
c-----PARAMETER vector
      INTEGER ntheta,ndep1
      INTEGER j,i,idep,jdep,index_i,index_j
      INTEGER index   
c-----
      COMPLEX csum
      COMPLEX cxmat(ndep1,ndep1),cm_inv(ndep1,ndep1)
      COMPLEX ctemp(ndep1,ndep1,ntheta),ctot(ndep1,ndep1)
      COMPLEX cx_diff(ndep1,ndep1)
      COMPLEX cfish(ntheta,ntheta),crb(ntheta,ntheta)
      INTEGER index2,index3,jgrad,igrad,ix,j1start,ibart,ierr
      COMPLEX czero
      data czero/(0.,0.)/
      REAL AVGSIGNAL,noise
      INTEGER ifrq
      DO igrad=1,ntheta
         DO jgrad=1,ntheta
            cfish(igrad,jgrad)= czero
         ENDDO
      ENDDO

      
      DO 200 ifrq=1,nfrq
         DO 200 ix=1,nx
            AVGSIGNAL = 0.0
            j1start=((ifrq-1))*ndep ! where the response starts
            ibart=ix+((ifrq-1))*nx ! which observation
            index3=(ibart-1)*ndep-1 
            DO idep=1,ndep
c     j=(ifrq-1)*ndep+idep
               index_i=(idep+j1start-1)*nx+ix
               DO jdep=1,ndep
                  index=ix+(jdep+j1start-1)*nx
                  index2=idep+(jdep+index3)*ndep 
                  cxmat(idep,jdep)=cov(index2)    
c     j=(ifrq-1)*ndep+jdep
                  index_j=(j-1)*nx+ix
                  index_j=(jdep+j1start-1)*nx+ix
                  cxmat(idep,jdep)=resp(index_I)*CONJG(resp(index_j))
               ENDDO            ! jdep
               AVGSIGNAL = AVGSIGNAL + ABS( resp(index_i ) ) ** 2
            ENDDO               ! idep
c     WRITE(*,*)' Adding white noise to the data in the *.obs file'
c     WRITE(*,'(a,f10.3,a)')'   with a  SNR =',snr_db,' dB'
c     WRITE(*,*)
            noise=avgsignal*10**(-snr_db/10)
            DO idep = 1, Ndep
               cxmat(idep,idep)=cxmat(idep,idep)+noise
            ENDDO

c---  Inverse coovariance
            CALL cminv(ndep,cxmat,cm_inv,ierr)
            IF (ierr.NE.0) THEN
               WRITE(*,*)' gnuncer: covariance matrix is singular'
               WRITE(*,*) 'ierr=',ierr
               STOP
            ENDIF
c     
c-----USE matrix multiplication to get CRB (baggeroer)
c     
            DO igrad=1,ntheta
               DO idep=1,ndep   ! trace
                  j=(ifrq-1)*ndep+idep
                  index_i=(j-1)*nx+ix
                  DO jdep=1,ndep ! matrix poduct
                     j=(ifrq-1)*ndep+jdep
                     index_j=(j-1)*nx+ix
                     cx_diff(idep,jdep)=
     &                    resp(index_i)*CONJG(grad(index_j,igrad))
     &                    +grad(index_i,igrad)*CONJG(resp(index_j))
                  ENDDO
               ENDDO
               CALL cmmul(cm_inv,ndep,ndep, 
     &              cx_diff,ndep,ctemp(1,1,igrad))
            ENDDO

            DO igrad=1,ntheta
               DO jgrad=1,ntheta
                  CALL cmmul(ctemp(1,1,igrad),ndep,ndep, 
     &                 ctemp(1,1,jgrad),ndep,ctot)
                  
                  csum=0     
                  DO i=1,ndep   ! trace
                     csum=csum+ctot(i,i)
                  ENDDO
c     WRITE(*,*)' matrix crb,', csum
                  cfish(igrad,jgrad)= cfish(igrad,jgrad)+csum
               ENDDO
            ENDDO
 200     CONTINUE
         
         CALL cminv(ntheta,cfish,crb,ierr)
         END

c***************************************************
      SUBROUTINE crb_det(ntheta,crb,cfish)
c     1       2     3    4    5      6
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comgrad.h'
      INCLUDE './oases/complo.f'
c-----PARAMETER vector
      INTEGER ntheta
      INTEGER j,idep,jdep,index_i,index_j
c-----
      COMPLEX csum
      COMPLEX cfish(ntheta,ntheta),crb(ntheta,ntheta)
      INTEGER index3,jgrad,igrad,ix,j1start,ibart,ierr
      COMPLEX czero
      data czero/(0.,0.)/
      REAL AVGSIGNAL,noise
      INTEGER ifrq
      DO igrad=1,ntheta
         DO jgrad=1,ntheta
            cfish(igrad,jgrad)= czero
         ENDDO
      ENDDO

      
      DO 200 ifrq=1,nfrq
         DO 200 ix=1,nx
            AVGSIGNAL = 0.0
            j1start=((ifrq-1))*ndep ! where the response starts
            ibart=ix+((ifrq-1))*nx ! which observation
            index3=(ibart-1)*ndep-1 
            DO idep=1,ndep
c     j=(ifrq-1)*ndep+idep
               index_i=(idep+j1start-1)*nx+ix
               AVGSIGNAL = AVGSIGNAL + ABS( resp(index_i ) ) ** 2
            ENDDO               ! idep
c     WRITE(*,*)' Adding white noise to the data in the *.obs file'
c     WRITE(*,'(a,f10.3,a)')'   with a  SNR =',snr_db,' dB'
c     WRITE(*,*)
            noise=avgsignal*10**(-snr_db/10)

c     
c-----USE matrix multiplication to get CRB (baggeroer)
c     
            DO igrad=1,ntheta
               DO jgrad=1,ntheta
                  csum=(0.,0.)
                  DO jdep=1,ndep ! matrix poduct
                     j=(ifrq-1)*ndep+jdep
                     index_j=(j-1)*nx+ix
                     csum=csum+
     &                    (grad(index_j,igrad))
     &                    *CONJG(grad(index_j,jgrad))
c     WRITE(*,*)'idep,gr,sum',csum,grad(index_j,igrad)
                  ENDDO
                  cfish(igrad,jgrad)= cfish(igrad,jgrad)
     &                 +REAL(csum)/noise
               ENDDO
            ENDDO


 200     CONTINUE
         
         CALL cminv(ntheta,cfish,crb,ierr)
         END
c***************************************************







