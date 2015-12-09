      SUBROUTINE plotsvd(vsvd,mpar1,ntheta)
C *** PLOT OF PDP DISTRIBUTIONS     
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './oases/complo.f'
      integer mpar1,ntheta
      real  vsvd(mpar1,mpar1)
      INTEGER nlab,igrid,nc,i,j
      REAL yval(2000),xoff
      REAL xdiv,ydiv,xdum,dx0,yoff
      REAL xmin,xmax,ymin,ymax
      REAL  XLEFT,XRIGHT,XAXIS,XINC ! data for the range axis
      REAL YUP,YDOWN,YAXIS,YINC
c      real x2,dx2,delpdp
c      integer j2,ndig,idum
c      integer np1,np2,ii
      INTEGER nxdif,nydif
      CHARACTER*6 OPTION(6),OPTTL(3)
      character*80 title1
      DATA OPTTL /' PSTCK','WTLRAN','UTLRAN'/
      data option /6*'      '/
c--- idum =1 is used for plotting figures using physical numbers
      xaxis=15
c     yaxis=10   
      yaxis=12

      OPTION(1)=PROGNA
      OPTION(2)=OPTTL(1)
      option(3)='    '
      option(4)='    '
c     IF (DEBUG) WRITE(prtfil,*) 'ENTERING PLTLOS'

c---  definition of xaxis 
      PTIT=' PDP-DISTRIBUTION'
      XMAX=ntheta
      XMIN=1
      CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
      xleft =xmin
      xright=xmax
      XDIV=1
c      xinc=0.2
      XTXT='  parameter$'
      PTIT='     '
c      WRITE(XTXT,820) NXDIF
c 820  FORMAT('range (10**',I3,')$')
c      NXDIF=1

c      I=MAX(INR,1)
      NLAB=0
      XTYP='LIN'
      YTYP='LIN'
      IGRID=0
c
c--- scale yaxis
c
      ymin=0
c      ymin=nc+1
      ymax=ntheta+1
      CALL AUTOAX(YMIN,YMAX,Ydown,YUP,YINC,YDIV,NYDIF)
      YTXT='Distribution Number$' 
      YTXT=' $'
      ydiv=1
      yup=ymax
      ydown=ymin
      yinc=1
      nc=ntheta
C *** WRITE PLP FILE
         title1=title
c      CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,Xaxis,Yaxis,
      CALL PLPWRI(OPTION,PTIT,TITLE1,NLAB,LAB,Xaxis,Yaxis,
     &                  IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &                  YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)

c
c**   write all the ppd's to file
c
      do i=1,ntheta
        yoff=i
        do j=1,ntheta
            yval(j)=vsvd(j,i)
        enddo
        xoff=1
        dx0=1
c       PLTWRI(N,XOFF,DX,YOFF,DY,X,IX,Y,IY)
        CALL PLTWRI(ntheta,xoff,dx0,yoff,0.,xdum,1,yval,1)
      enddo


C *** FORMATS
 811  FORMAT('SD:',F9.1,' m$')
 812  FORMAT('RD:',F9.1,' m$')
      RETURN
      END  
