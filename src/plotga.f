      SUBROUTINE PLT_tim()
C     *** plot of multi freq inversion     
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/complo.f'
      INTEGER nlab,igrid,nc,i,jj,idep,ifrq,nx_plotdim
      PARAMETER(nx_plotdim=5000)
      REAL xval(nx_plotdim),yvalresp(nx_plotdim),yvaldata(nx_plotdim)
      REAL xdiv,ydiv
      REAL xmin,xmax,ymin,ymax
      REAL  XLEFT,XRIGHT,XAXIS,XINC ! data for the range axis
      REAL YUP,YDOWN,YAXIS,YINC
      REAL dfrq
      INTEGER ipow,ii,IF,lfr,mfr,ntime,nfr
c     lfr and mfr correspond to HS lx and mx.
      INTEGER nxdif,nydif,index,jjplot
      CHARACTER*6 OPTION(6),OPTTL(3)
      DATA option /6*'      '/
      DATA OPTTL /'STLRAN','WTLRAN','UTLRAN'/
c     COMPLEX tc(5000)
c     REAL tr(10000)
c     EQUIVALENCE (tc(1),tr(1))


      xaxis=15                  ! 15cm recommended by finn
      yaxis=10
c     WRITE(*,*)'entering plt_tim..., iopt(1)',iopt(1)
      IF (nx.GT.nx_plotdim) STOP ' dimension in plotga'
      OPTION(1)=PROGNa
      
      OPTION(2)=OPTTL(1)
      option(3)=',DSD,U'
      option(4)='NI,   '
      
      dfrq=frq(2)-frq(1)
      xincre=1/(frq(nfrq)-frq(1))/2
      Lfr=NINT(Frq(1)/DFRQ+1)
      mfr=lfr+nfrq
      IPOW=0
 990  IPOW=IPOW+1
      II=2**IPOW
      IF (II.LT.mfr) GO TO 990
      Nfr=II
      ntime=nfr*2


c---  definition of xaxis for time

      XMAX=1/dfrq
      XMIN=0
      xminimum=0
      xmaxmum=1/dfrq
      xincre=1./(frq(nfrq)-frq(1))
c---  definition for freq
      ntime=nfrq
      xmin=frq(1)
      xmax=frq(nfrq)
      
      xmaxmum=xmax
      xminimum=xmin
      xincre=(frq(2)-frq(1))
      WRITE(*,*)'xminimum,xmaxmum,xincre'
      WRITE(*,*)xminimum,xmaxmum,xincre
      CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
      PTIT=' Multi freq'
      WRITE(XTXT,840) 
 840  FORMAT('Time (s) $')
      XTXT=' Frequency (Hz) $'
      DO i=1,ntime
         xval(i)=(xminimum + (xincre)*(i-1))
      ENDDO

      WRITE(13,'(a,i5,a)')'nplots=',ncurv,';'
c     WRITE(13,'(a,i5,a)')'nphone=', ndep,';'
      WRITE(13,'(a,i5,a)')'nphone=', ntime,';'
      IF (ntime.GT.nx_plotdim) STOP ' dimension in plotga'
      
      NLAB=0
      XTYP='LIN'
      YTYP='LIN'
      IGRID=0
      NC=2
c     
c---  for the computed response find yaxis
c     
      jjplot=0
      jj=0
      WRITE(*,*)'nfrq,ndep',nfrq,ndep
      DO 1000 i=1,nx
         DO 1000 idep=1,ndep
            jj=1+jj
            
c     I=MAX(INR,1)
            NLAB=2
            WRITE(LAB(1),810) jj
 810        FORMAT('Curve no:',i3)
c     WRITE(LAB(2),811) FRQ(ifrq)
c     811     FORMAT('Freq:',F7.0,' Hz$')
c     IF (iopt(1).NE.3) THEN
c     NLAB=3
            WRITE(LAB(3),812) rdep(idep)
 812        FORMAT('Depth:',F7.1,' m$')
c     ENDIF
            
            DO  IF=1,ntime
               yvalresp(IF)=0.
               yvaldata(IF)=0.
            ENDDO
c     
            DO ifrq=1,nfrq
               index= (idep +((ifrq-1))*ndep-1)*nx
               IF=2*(ifrq-1+lfr)
c     WRITE(*,*)'index,if,ifrq,ntime',index,IF,ifrq,ntime
               yvalresp(IF-1) = REAL(resp(i+index))
               yvalresp(IF  ) = AIMAG(resp(i+index))
               yvaldata(IF-1) = REAL(DATA(i+index))
               yvaldata(IF  ) = AIMAG(DATA(i+index))

               WRITE(99,*) REAL(resp(i+index)),AIMAG(resp(i+index))
               WRITE(98,*) REAL(DATA(i+index)),AIMAG(DATA(i+index))
            ENDDO
            DO ifrq=1,nfrq
               index= (idep +((ifrq-1))*ndep-1)*nx
               IF=2*(ifrq-1+lfr)
               yvalresp(ifrq) = 20 * LOG10(ABS(resp(i+index)))
               yvaldata(ifrq) = 20 * LOG10(ABS(DATA(i+index)))
            ENDDO
c     

c     CALL RFFT(yvalresp,Ntime,-1)
c     CALL RFFT(yvaldata,Ntime,-1)
c     CALL RFFT(tc,Ntime,-1)

c     
c---  scale yaxis
            
            ymin= 1e10
            ymax=-1e10
            DO IF=1,ntime
               IF (ymin.GT.yvalresp(IF)) ymin=yvalresp(IF)
               IF (ymax.LT.yvalresp(IF)) ymax=yvalresp(IF)
               IF (ymin.GT.yvaldata(IF)) ymin=yvaldata(IF)
               IF (ymax.LT.yvaldata(IF)) ymax=yvaldata(IF)
            ENDDO
            WRITE(*,*)'ymin,ymax',ymin,ymax
            
            CALL AUTOAX(YMIN,YMAX,Ydown,YUP,YINC,YDIV,NYDIF)
            
            YTXT='Y- value$'
            
C     *** WRITE PLP FILE
c     WRITE(*,*)'title',title
            CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,Xaxis,Yaxis,
     &           IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &           YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)
            
            
c**   WRITE the observed DATA to the plt file
            CALL PLTWRI(ntime,0.,0.,0.,0.,xval,1,yvaldata,1)
            
c**   WRITE synthetic DATA to plt file
            CALL PLTWRI(ntime,0.,0.,0.,0.,xval,1,yvalresp,1)
            
            jjplot=jjplot+1
            WRITE(13,'(a,i3,a,4f15.6,a)')
     &           ' plaxis(:,',jjplot,')=[',xleft,
     &           xright,ydown,yup,']'';'
            
 1000    CONTINUE               ! loop over ncurv

         RETURN
         END


      SUBROUTINE PLTLOS()
C     *** TRANSMISSION LOSS VS RANGE     
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/complo.f'
      INTEGER nlab,igrid,nc,i,jj,idep,ifrq,nx_plotdim
      PARAMETER(nx_plotdim=5000)
      REAL xval(nx_plotdim),yvalresp(nx_plotdim),yvaldata(nx_plotdim)
      REAL xdiv,ydiv,e1,e2,ahelp
      REAL xmin,xmax,ymin,ymax
      REAL  XLEFT,XRIGHT,XAXIS,XINC ! data for the range axis
      REAL YUP,YDOWN,YAXIS,YINC

      INTEGER nxdif,nydif,index,jjplot
      CHARACTER*6 OPTION(6),OPTTL(3)
      DATA option /6*'      '/
      DATA OPTTL /'STLRAN','WTLRAN','UTLRAN'/
      xaxis=15                  ! 15cm recommended by finn
      yaxis=10
      e1=0.
      e2=0.
      WRITE(*,*)'entering pltlos..., iopt(1)',iopt(1)
      IF (nx.GT.nx_plotdim) STOP ' dimension in plotga'
      OPTION(1)=PROGNa

      OPTION(2)=OPTTL(1)
      option(3)=',DSD,U'
      option(4)='NI,   '
      WRITE(13,'(a,i5,a)')'nplots=',ncurv,';'
c     WRITE(13,'(a,i5,a)')'nphone=', ndep,';'
      WRITE(13,'(a,i5,a)')'nphone=', nx,';'

      IF (iopt(1).EQ.1) THEN
c---  definition of xaxis for slowness
         PTIT=' Slowness Kernel'
         XMAX=1.0/cmin
         XMIN=1.0/cmax
         CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
         WRITE(XTXT,820) NXDIF
 820     FORMAT('Slowness (10**',I3,')$')
         DO i=1,nx
            xval(i)=(1/cmax + (1/cmin-1/cmax)/(nwave-1)*(i-1))
         ENDDO
      ELSE IF (iopt(1).EQ.2) THEN
c---  definition of xaxis for transmission loss
         PTIT='TRANSMISSION LOSS'
         XTXT='Range (km)$'
         XMAX=-1000
         XMIN=100000
         DO i=1,nx
            xval(i)=xranges(i)*0.001
            xmax=MAX(xmax,xval(i))
            xmin=MIN(xmin,xval(i))
         ENDDO
         WRITE(*,*)'xmin,xmax,xranges(1),xranges(2)'
         WRITE(*,*)xmin,xmax,xranges(1),xranges(2)
         CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
         xdiv=1

      ELSEIF (iopt(1).EQ.3) THEN
c---  definition of xaxis for Angle
         OPTION(2)='P-P'
         XMAX=angle2
         XMIN=angle1
         CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
         Xright=MIN(xright,90.)
         PTIT='Reflection coefficient'
         WRITE(XTXT,830) 
 830     FORMAT('Grazing angle $')
         DO i=1,nx
            xval(i)=(angle1 + (angle2-angle1)/(nang-1)*(i-1))
         ENDDO
      ELSEIF (iopt(1).EQ.4) THEN
c---  definition of xaxis for time
         OPTION(2)='REV'
         XMAX=xmaxmum
         XMIN=xminimum
         WRITE(*,*)'xminimum,xmaxmum,xincre'
         WRITE(*,*)xminimum,xmaxmum,xincre
         CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
         Xright=MIN(xright,90.)
         PTIT='Reverberation level'
         WRITE(XTXT,840) 
 840     FORMAT('Time (s) $')
         DO i=1,nx
            xval(i)=(xminimum + (xincre)*(i-1))
         ENDDO
      ENDIF

      NLAB=0
      XTYP='LIN'
      YTYP='LIN'
      IGRID=0
      NC=2
c     
c---  for the computed response find yaxis
c     
      jjplot=0
      jj=0
      WRITE(*,*)'nfrq,ndep',nfrq,ndep
      DO 1000 ifrq=1,nfrq
         DO 1000 idep=1,ndep
            jj=1+jj

c     I=MAX(INR,1)
            NLAB=2
            WRITE(LAB(1),810) jj
 810        FORMAT('Curve no:',i3)
            WRITE(LAB(2),811) FRQ(ifrq)
 811        FORMAT('Freq:',F7.0,' Hz$')
            IF (iopt(1).NE.3) THEN
               NLAB=3
               WRITE(LAB(3),812) rdep(idep)
 812           FORMAT('Depth:',F7.1,' m$')
            ENDIF

            index=(jj-1)*nx
            IF (iopt(1).EQ.2) THEN
               DO i=1,nx
c     yvalresp(i)=20*alog10(ABS(resp(i,jj)))
                  yvalresp(i)=REAL(resp(i+index))
                  e1=e1+yvalresp(i)
c     WRITE(*,*)'plotga:yvalresp',yvalresp(i)
               ENDDO
            ELSE
               DO i=1,nx
                  yvalresp(i)=REAL(resp(i+index))
c     WRITE(*,*)'plotga:yvalresp',yvalresp(i)
                  e1=e1+yvalresp(i)
               ENDDO
            ENDIF   
c     
c---  scale yaxis
c     
c-------------------------REAL or Synt DATA ------------------------
            IF ((iopt(1).EQ.2).AND.(itrans(3).NE.2)) THEN
c--   range-depth   and observed resp is NOT on dB scale
               DO i=1,nx
c     yvaldata(i)=20*alog10(ABS(DATA(i+index)))
                  yvaldata(i)=REAL(DATA(i+index))
                  e2=e2+yvaldata(i)
               ENDDO
            ELSE  
               DO i=1,nx
                  yvaldata(i)=REAL(DATA(i+index))
c     WRITE(*,*)'plotga:yvaldata',yvaldata(i)
                  e2=e2+yvaldata(i)
               ENDDO
            ENDIF 
            

c     234567
            IF ((iopt(1).EQ.3)
     *           .AND. (itrans(3).NE.2)) THEN
c     WRITE(*,*) ' am I here ?'
               DO i=1,nx
                  yvaldata(i)= -20. *alog10 (yvaldata(i))
                  yvalresp(i)= -20. *alog10 (yvalresp(i))
               ENDDO
            ENDIF
            ymin= 1e10
            ymax=-1e10
            DO i=1,nx
               IF (ymin.GT.yvalresp(i)) ymin=yvalresp(i)
               IF (ymax.LT.yvalresp(i)) ymax=yvalresp(i)
               IF (ymin.GT.yvaldata(i)) ymin=yvaldata(i)
               IF (ymax.LT.yvaldata(i)) ymax=yvaldata(i)
            ENDDO
            WRITE(*,*)'ymin,ymax',ymin,ymax

            CALL AUTOAX(YMIN,YMAX,Ydown,YUP,YINC,YDIV,NYDIF)
            IF (itrans(3).EQ.2 .AND.iopt(1).NE.4) THEN
               ahelp=ydown
               ydown=yup
               yup=ahelp
               yinc=-yinc
            ENDIF

            IF (iopt(1).EQ.3) THEN
c---  Reflection loss
               YTXT='Loss (dB)$'
c     YTXT='Loss (lin)$'
               ydiv=1
            ELSEIF (iopt(1).EQ.4) THEN
c---  Reflection loss
               YTXT='Signal level (dB)$'
c     YTXT='Loss (lin)$'
               ydiv=1
            ELSEIF (iopt(1).EQ.2) THEN
c---  transmission loss
               IF (iopt(7).EQ.1) THEN
                  YTXT='Y- value$'
               ELSEIF (iopt(7).EQ.2) THEN
                  YTXT='Vert. particle velocity (dB//1m/s)$'
               ELSEIF (iopt(7).EQ.3) THEN
                  YTXT='Hor. particle velocity (dB//1m/s)$'
               ENDIF
               Ydiv=1
            ELSE
               IF (iopt(7).EQ.1) THEN
                  WRITE(YTXT,'(''Normal stress (10**'',i3,'')$'')') 
     &                 NYDIF
               ELSEIF (iopt(7).EQ.2) THEN
                  WRITE(YTXT,'(''Vert velocity (10**'',i3,'')$'')')
     &                 NYDIF
               ELSEIF (iopt(7).EQ.3) THEN
                  WRITE(YTXT,'(''Hor  velocity (10**'',i3,'')$'')')
     &                 NYDIF
               ENDIF
            ENDIF

C     *** WRITE PLP FILE
c     WRITE(*,*)'title',title
            CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,Xaxis,Yaxis,
     &           IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &           YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)


c**   WRITE the observed DATA to the plt file
            CALL PLTWRI(nx,0.,0.,0.,0.,xval,1,yvaldata,1)

c**   WRITE synthetic DATA to plt file
            CALL PLTWRI(nx,0.,0.,0.,0.,xval,1,yvalresp,1)
            
            jjplot=jjplot+1
            WRITE(13,'(a,i3,a,4f15.6,a)')
     &           ' plaxis(:,',jjplot,')=[',xleft,xright,
     &           ydown,yup,']'';'

 1000    CONTINUE               ! loop over ncurv

         RETURN
         END
      SUBROUTINE PLTLOS_dep()
C     *** TRANSMISSION LOSS VS RANGE     
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/complo.f'
      INTEGER nlab,igrid,nc,i,jj,ifrq,nx_plotdim
      PARAMETER(nx_plotdim=5000)
      REAL xval(nx_plotdim),yvalresp(nx_plotdim),yvaldata(nx_plotdim)
      REAL xdiv,ydiv,e1,e2,ahelp
      REAL xmin,xmax,ymin,ymax
      REAL  XLEFT,XRIGHT,XAXIS,XINC ! data for the range axis
      REAL YUP,YDOWN,YAXIS,YINC
      INTEGER nxdif,nydif,index,jjplot,ibart,index3,j1start,ix
      CHARACTER*6 OPTION(6),OPTTL(3)
      DATA option /6*'      '/
      DATA OPTTL /'STLRAN','WTLRAN','UTLRAN'/
      xaxis=15                  ! 15cm recommended by finn
      yaxis=10
      e1=0.
      e2=0.
      WRITE(*,*)'entering pltlos_dep..., iopt(1)',iopt(1)
      IF (nx.GT.nx_plotdim) STOP ' dimension in plotga'
      OPTION(1)=PROGNa

      OPTION(2)=OPTTL(1)
      option(3)=',DSD,U'
      option(4)='NI,   '
      WRITE(13,'(a,i5,a)')'nplots=',ncurv,';'
      WRITE(13,'(a,i5,a)')'nphone=', ndep,';'
c     WRITE(13,'(a,i5,a)')'nphone=', nx,';'

      IF (iopt(1).EQ.1) THEN
c---  definition of xaxis for slowness
      ELSE IF (iopt(1).EQ.2) THEN
c---  definition of xaxis for transmission loss
         PTIT='TRANSMISSION LOSS'
         XTXT='Range (km)$'
         XMAX=-1000
         XMIN=100000
         DO i=1,ndep
            xval(i)=rdep(i)
            xval(i)=i
            xmax=MAX(xmax,xval(i))
            xmin=MIN(xmin,xval(i))
         ENDDO
         WRITE(*,*)'xmin,xmax,xranges(1),xranges(2)'
         WRITE(*,*)xmin,xmax,xranges(1),xranges(2)
         CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
         xdiv=1

      ELSEIF (iopt(1).EQ.3) THEN
c---  definition of xaxis for Angle
      ELSEIF (iopt(1).EQ.4) THEN
c---  definition of xaxis for time
      ENDIF

      NLAB=0
      XTYP='LIN'
      YTYP='LIN'
      IGRID=0
      NC=2
c     
c---  for the computed response find yaxis
c     
      jjplot=0
      jj=0
      WRITE(*,*)'nfrq,ndep,nx',nfrq,ndep,nx
      DO 1000 ifrq=1,nfrq
         DO 1000 ix=1,nx
            jj=1+jj
            ibart=ix+((ifrq-1))*nx ! which observation
            index3=(ibart-1)*ndep-1 
            j1start=((ifrq-1))*ndep ! where the response starts
c     I=MAX(INR,1)
            NLAB=2
            WRITE(LAB(1),810) jj
 810        FORMAT('Curve no:',i3)
            WRITE(LAB(2),811) FRQ(ifrq)
 811        FORMAT('Freq:',F7.0,' Hz$')
c     IF (iopt(1).NE.3) THEN
c     NLAB=3
c     WRITE(LAB(3),812) rdep(idep)
c     812    FORMAT('Depth:',F7.1,' m$')
c     ENDIF


            IF (iopt(1).EQ.2) THEN
               DO i=1,ndep
                  index=ix+(i+j1start-1)*nx
                  yvalresp(i)=REAL(resp(index))
                  e1=e1+yvalresp(i)
               ENDDO
            ELSE
               DO i=1,ndep
                  index=ix+(i+j1start-1)*nx
                  yvalresp(i)=REAL(resp(index))
c     WRITE(*,*)'plotga:yvalresp',yvalresp(i)
                  e1=e1+yvalresp(i)
               ENDDO
            ENDIF   
c     
c---  scale yaxis
c     
c-------------------------REAL or Synt DATA ------------------------
            IF ((iopt(1).EQ.2).AND.(itrans(3).NE.2)) THEN
c--   range-depth   and observed resp is NOT on dB scale
               DO i=1,ndep
                  index=ix+(i+j1start-1)*nx
                  yvaldata(i)=ABS(DATA(index))
                  e2=e2+yvaldata(i)
               ENDDO
            ELSE
               DO i=1,ndep
                  index=ix+(i+j1start-1)*nx
                  yvaldata(i)=REAL(DATA(index))
c     WRITE(*,*)'plotga:yvaldata',yvaldata(i)
                  e2=e2+yvaldata(i)
               ENDDO
            ENDIF 

c     234567
            IF ((iopt(1).EQ.3)
     *           .AND. (itrans(3).NE.2)) THEN
               DO i=1,ndep
                  yvaldata(i)= -20. *alog10 (yvaldata(i))
                  yvalresp(i)= -20. *alog10 (yvalresp(i))
               ENDDO
            ENDIF
            ymin= 1e10
            ymax=-1e10
            DO i=1,ndep
               IF (ymin.GT.yvalresp(i)) ymin=yvalresp(i)
               IF (ymax.LT.yvalresp(i)) ymax=yvalresp(i)
               IF (ymin.GT.yvaldata(i)) ymin=yvaldata(i)
               IF (ymax.LT.yvaldata(i)) ymax=yvaldata(i)
            ENDDO
            WRITE(*,*)'ymin,ymax',ymin,ymax

            CALL AUTOAX(YMIN,YMAX,Ydown,YUP,YINC,YDIV,NYDIF)
            IF (itrans(3).EQ.2 .AND.iopt(1).NE.4) THEN
               ahelp=ydown
               ydown=yup
               yup=ahelp
               yinc=-yinc
            ENDIF

            IF (iopt(1).EQ.3) THEN
c---  Reflection loss
               YTXT='Loss (dB)$'
c     YTXT='Loss (lin)$'
               ydiv=1
            ELSEIF (iopt(1).EQ.4) THEN
c---  Reflection loss
               YTXT='Signal level (dB)$'
c     YTXT='Loss (lin)$'
               ydiv=1
            ELSEIF (iopt(1).EQ.2) THEN
c---  transmission loss
               IF (iopt(7).EQ.1) THEN
                  YTXT='Y- value$'
               ELSEIF (iopt(7).EQ.2) THEN
                  YTXT='Vert. particle velocity (dB//1m/s)$'
               ELSEIF (iopt(7).EQ.3) THEN
                  YTXT='Hor. particle velocity (dB//1m/s)$'
               ENDIF
               Ydiv=1
            ELSE
               IF (iopt(7).EQ.1) THEN
                  WRITE(YTXT,'(''Normal stress (10**'',i3,'')$'')')
     &                 NYDIF
               ELSEIF (iopt(7).EQ.2) THEN
                  WRITE(YTXT,'(''Vert velocity (10**'',i3,'')$'')')
     &                 NYDIF
               ELSEIF (iopt(7).EQ.3) THEN
                  WRITE(YTXT,'(''Hor  velocity (10**'',i3,'')$'')')
     &                 NYDIF
               ENDIF
            ENDIF

C     *** WRITE PLP FILE
c     WRITE(*,*)'title',title
            CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,Xaxis,Yaxis,
     &           IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &           YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)


c**   WRITE the observed DATA to the plt file
            CALL PLTWRI(ndep,0.,0.,0.,0.,xval,1,yvaldata,1)

c**   WRITE synthetic DATA to plt file
            CALL PLTWRI(ndep,0.,0.,0.,0.,xval,1,yvalresp,1)
            
            jjplot=jjplot+1
            WRITE(13,'(a,i3,a,4f15.6,a)')
     &           ' plaxis(:,',jjplot,')=[',xleft,xright,
     &           ydown,yup,']'';'

 1000    CONTINUE               ! loop over ncurv

         RETURN
         END                    !pltlos_dep

c*******************************************************
      SUBROUTINE PLTPHASEF
C     *** PLOT OF PHASES ACROSS FREQUNCY BAND    
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/complo.f'
      INTEGER nlab,igrid,nc,i,jj,jdep,ifrq,id,nx_plotdim
      PARAMETER(nx_plotdim=5000)     
      REAL xvalresp(nx_plotdim),xvaldata(nx_plotdim)
      REAL xdiv,ydiv,e1,e2
      REAL xmin,xmax,ymin,ymax
      REAL  XLEFT,XRIGHT,XAXIS,XINC ! data for the range axis
      REAL YUP,YDOWN,YAXIS,YINC

      INTEGER nxdif,nydif,index
      CHARACTER*6 OPTION(6),OPTTL(3)
      INTEGER    j,jjplot
      INTEGER ix 
      DATA option /6*'      '/
      DATA OPTTL /'STLRAN','WTLRAN','UTLRAN'/
      WRITE(*,*)'entering into pltphaseF....'
      xaxis=15                  ! 15cm recommended by finn
      yaxis=10
      WRITE(13,'(a,i5,a)')'nplots=',
     1     ncurv,';'
      WRITE(13,'(a,i5,a)')'nphone=', nfrq,';'
      WRITE(13,'(a,i5,a)')'nfrq=', nfrq,';'
      IF (nx.GT.nx_plotdim) STOP ' dimension in plotga'
      OPTION(1)=PROGNa

      OPTION(2)=OPTTL(1)
      option(3)=',DSD,U'
      option(4)='NI,   '


      YTXT='Depth$'
      NLAB=0
      XTYP='LIN'
      YTYP='LIN'
      IGRID=0
c     
c---- scale xaxis
c     
      Ymin=Frq(1)
      ymax=Frq(nfrq)

      CALL AUTOAX(YMIN,YMAX,Ydown,YUP,YINC,YDIV,NYDIF)
      e1=ydown
      ydown=yup
      yup=e1  

c     
c---  for the computed response find yaxis
c     
      jj=0
      jjplot=0
      DO jdep=1,ndep
         DO  ix=1,nx
c     WRITE(*,*)'******************** range',ix
            NC=2
            jj=1+jj
            
            NLAB=2
            WRITE(LAB(1),810) jj
 810        FORMAT('Curve no:',i3)
c     WRITE(LAB(2),811) rdep(jdep)
c     811        FORMAT('Depth:',F7.1,' Hz$')

c*************************************************************
c     magnitude  plot.
c****************************************************
c---- calculated DATA
            e1=0
            e2=0
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               xvalresp(ifrq)=ABS(resp(ix+index))
               e1=e1+ABS(resp(ix+index))**2 ! resp
               e2=e2+ABS(DATA(ix+index))**2 ! data
               xvaldata(ifrq)=ABS(DATA(ix+index))
            ENDDO
            e1=SQRT(e1)
            e2=SQRT(e2)
            DO id=1,nfrq
               xvaldata(id)= xvaldata(id)/e2
               xvalresp(id)= xvalresp(id)/e1
            ENDDO
c     
            xmin= 1e10
            xmax=-1e10
            DO i=1,nfrq
               IF (xmin.GT.xvalresp(i)) xmin=xvalresp(i)
               IF (xmax.LT.xvalresp(i)) xmax=xvalresp(i)
               IF (xmin.GT.xvaldata(i)) xmin=xvaldata(i)
               IF (xmax.LT.xvaldata(i)) xmax=xvaldata(i)
            ENDDO
            CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
            WRITE(XTXT,820) NXDIF
 820        FORMAT('Magnitude (10**',I3,')$')
            PTIT='Magnitude'
c     
            
C     *** WRITE PLP FILE
            CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,Xaxis,Yaxis,
     &           IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &           YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)
            
c***  WRITE PLM FILe
            jjplot=jjplot+1
            WRITE(13,'(a,i3,a,4f15.6,a)')
     &           ' plaxis(:,',jjplot,')=[',
     &           xleft,xright,ydown,yup,']'';'
            
c**   WRITE the observed DATA to the plt file
            CALL PLTWRI(nfrq,0.,0.,0.,0.,xvaldata,1,frq,1)
            
c**   WRITE synthetic DATA to plt file
            CALL PLTWRI(nfrq,0.,0.,0.,0.,xvalresp,1,frq,1)
c     

         ENDDO 
      ENDDO
      RETURN
      END

C**************************************************************
c*******************************************************

c*******************************************************
      SUBROUTINE PLTPHASEk
C     *** PLOT OF PHASES ACROSS range of array   
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/complo.f'
      INTEGER nlab,igrid,nc,i,jj,jdep,ifrq,id,nx_plotdim
      PARAMETER(nx_plotdim=5000)     
      REAL xvalresp(nx_plotdim),xvaldata(nx_plotdim)
      REAL xdiv,ydiv,e1,e2
      REAL xmin,xmax,ymin,ymax
      REAL  XLEFT,XRIGHT,XAXIS,XINC ! data for the range axis
      REAL YUP,YDOWN,YAXIS,YINC

      INTEGER nxdif,nydif,index
      CHARACTER*6 OPTION(6),OPTTL(3)
      INTEGER    j,jjplot
      INTEGER ix 
      DATA option /6*'      '/
      DATA OPTTL /'STLRAN','WTLRAN','UTLRAN'/
      WRITE(*,*)'entering into pltphasek....'
      xaxis=15                  ! 15cm recommended by finn
      yaxis=10
      WRITE(13,'(a,i5,a)')'nplots=',
     1     ncurv,';'
      WRITE(13,'(a,i5,a)')'nphone=', nx,';'
      WRITE(13,'(a,i5,a)')'nfrq=', nfrq,';'
      IF (nx.GT.nx_plotdim) STOP ' dimension in plotga'
      OPTION(1)=PROGNa

      OPTION(2)=OPTTL(1)
      option(3)=',DSD,U'
      option(4)='NI,   '


      YTXT='Depth$'
      NLAB=0
      XTYP='LIN'
      YTYP='LIN'
      IGRID=0
c     
c---- scale xaxis
c     
      Ymin=xranges(1)
      ymax=xranges(nx)

      CALL AUTOAX(YMIN,YMAX,Ydown,YUP,YINC,YDIV,NYDIF)
      e1=ydown
      ydown=yup
      yup=e1  

c     
c---  for the computed response find yaxis
c     
      jj=0
      jjplot=0
      DO jdep=1,ndep
         DO  ifrq=1,nfrq
c     WRITE(*,*)'******************** range',ix
            NC=2
            jj=1+jj
            
            NLAB=2
            WRITE(LAB(1),810) jj
 810        FORMAT('Curve no:',i3)
c     WRITE(LAB(2),811) rdep(jdep)
c     811        FORMAT('Depth:',F7.1,' Hz$')

c*************************************************************
c     magnitude  plot.
c****************************************************
c---- calculated DATA
            e1=0
            e2=0
            DO ix=1,nx
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               xvalresp(ix)=ABS(resp(ix+index))
               e1=e1+ABS(resp(ix+index))**2 ! resp
               e2=e2+ABS(DATA(ix+index))**2 ! data
               xvaldata(ix)=ABS(DATA(ix+index))
            ENDDO
            e1=SQRT(e1)
            e2=SQRT(e2)
            DO id=1,nx
               xvaldata(id)= xvaldata(id)/e2
               xvalresp(id)= xvalresp(id)/e1
            ENDDO
c     
            xmin= 1e10
            xmax=-1e10
            DO i=1,nx
               IF (xmin.GT.xvalresp(i)) xmin=xvalresp(i)
               IF (xmax.LT.xvalresp(i)) xmax=xvalresp(i)
               IF (xmin.GT.xvaldata(i)) xmin=xvaldata(i)
               IF (xmax.LT.xvaldata(i)) xmax=xvaldata(i)
            ENDDO
            CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
            WRITE(XTXT,820) NXDIF
 820        FORMAT('Magnitude (10**',I3,')$')
            PTIT='Magnitude'
c     
            
C     *** WRITE PLP FILE
            CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,Xaxis,Yaxis,
     &           IGRID,
     &           YDOWN,YUP   ,YINC,YDIV,YTXT,YTYP,
     &           XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &           NC)
            
c***  WRITE PLM FILe
            jjplot=jjplot+1
            WRITE(13,'(a,i3,a,4f15.6,a)')
     &           ' plaxis(:,',jjplot,')=[',
     &           ydown,yup,xleft,xright,']'';'
c     &           xleft,xright,ydown,yup,']'';'
            
c**   WRITE the observed DATA to the plt file
c     CALL PLTWRI(nx,0.,0.,0.,0.,xvaldata,1,xranges ,1)
            CALL PLTWRI(nx,0.,0.,0.,0.,xranges ,1,xvaldata,1)
            
c**   WRITE synthetic DATA to plt file
c     CALL PLTWRI(nx,0.,0.,0.,0.,xvalresp,1,xranges ,1)
            CALL PLTWRI(nx,0.,0.,0.,0.,xranges ,1,xvalresp,1)
c     

         ENDDO 
      ENDDO
      RETURN
      END

C**************************************************************
c*******************************************************
      SUBROUTINE PLTPHASE
C     *** PLOT OF PHASES ACROSS THE ARRAY    
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/complo.f'
      INTEGER nlab,igrid,nc,i,jj,idep,ifrq,id,idum,nx_plotdim
      PARAMETER(nx_plotdim=5000)     
      REAL xvalresp(nx_plotdim),xvaldata(nx_plotdim)
      REAL xdiv,ydiv,e1,e2
      REAL xmin,xmax,ymin,ymax
      REAL  XLEFT,XRIGHT,XAXIS,XINC ! data for the range axis
      REAL YUP,YDOWN,YAXIS,YINC

      INTEGER nxdif,nydif,index
      CHARACTER*6 OPTION(6),OPTTL(3)
      INTEGER    info,j,index2,jjplot
      COMPLEX*16 zx(mdep,mdep), ssvd(mdep+1)  ! matrix to find SVD of 
      COMPLEX*16 esvd(mdep),usvd(mdep,mdep),vsvd(mdep,mdep)
      COMPLEX*16 work2(mdep)           
      COMPLEX ch,asum
      INTEGER ix, ibart, j1start,index3
      INTEGER iphase(3)
      DATA  iphase /1, 0, 0/ ! thus plotting magnitude, but not phase and power
      DATA option /6*'      '/
      DATA OPTTL /'STLRAN','WTLRAN','UTLRAN'/
      INTEGER xplane
      LOGICAL out_plane
      REAL arrayshape(10)
      COMMON /arrayparm/out_plane,arrayshape,xplane
      REAL yarray(mdep)
      IF (iopt(10).GE.2) iphase(2)=1 ! plot phase match
      IF (iopt(10).GE.3) iphase(3)=1 ! plot bartlett
      WRITE(*,*)'entering into pltphase....'
      xaxis=15                  ! 15cm recommended by finn
      yaxis=10
      WRITE(13,'(a,i5,a)')'nplots=',
     1     (iphase(1)+iphase(2)+iphase(3))*ncurv,';'
      WRITE(13,'(a,i5,a)')'nphone=', ndep,';'
      IF (nx.GT.nx_plotdim) STOP ' dimension in plotga'
      OPTION(1)=PROGNa

      OPTION(2)=OPTTL(1)
      option(3)=',DSD,U'
      option(4)='NI,   '

      YTXT='Depth$'
      NLAB=0
      XTYP='LIN'
      YTYP='LIN'
      IGRID=0
      IF ( out_plane) THEN
         DO i=1,ndep
            yarray(i)=i
         ENDDO
      ELSE
         DO i=1,ndep
            yarray(i)=rdep(i)
         ENDDO
      ENDIF
c     
c---- scale xaxis
c     
      Ymin=MIN(yarray(1),yarray(ndep))
      ymax=MAX(yarray(1),yarray(ndep))

      CALL AUTOAX(YMIN,YMAX,Ydown,YUP,YINC,YDIV,NYDIF)
      e1=ydown
      ydown=yup
      yup=e1  

c     
c---  for the computed response find yaxis
c     
      jj=0
      jjplot=0
      DO 1000 ix=1,nx
         WRITE(*,*)'******************** range',ix
         DO 1000 ifrq=1,nfrq
            ibart=ix+((ifrq-1))*nx ! which observation
            index3=(ibart-1)*ndep-1 
            j1start=((ifrq-1))*ndep ! where the response starts
            NC=2
            jj=1+jj
            IF (iopt(13).EQ.1) THEN      
c---  find the first eigenvector
               DO j=1,ndep
                  DO i=1,ndep
c     index2=i+ (j+(jj-1)*ndep-1)*ndep 
                     index2=i+(j+index3)*ndep 
                     zx(i,j)=cov(index2)             
                  ENDDO 
               ENDDO 
               
               CALL zsvdc(zx,mdep,ndep,ndep,ssvd,esvd,usvd,
     1              mdep,vsvd,mdep,work2,21,info)
               WRITE(*,*)' Largest SVD:',ssvd(1)
               IF (info.NE.0) THEN
                  WRITE(*,*)'********* info from zsvdc ************',
     &                 info
                  STOP
               ENDIF
            ENDIF
            NLAB=2
            WRITE(LAB(1),810) jj
 810        FORMAT('Curve no:',i3)
            WRITE(LAB(2),811) FRQ(ifrq)
 811        FORMAT('Freq:',F7.1,' Hz$')

c*************************************************************
c     magnitude plot.
c****************************************************
            IF (iphase(1).EQ.1) THEN
c---- calculated DATA
               e1=0
               e2=0
               DO j=1,ndep
c     index=(id +(ifrq-1)*ndep-1)*nx
                  index=(j+j1start-1)*nx
                  xvalresp(j)=ABS(resp(ix+index))
                  e1=e1+xvalresp(j)
               ENDDO
c---  scale DATA
               DO id=1,ndep
                  IF (iopt(13).EQ.1) THEN
                     xvaldata(id)=ABS(usvd(id,1))
                  ELSE
                     index=(id +j1start-1)*nx
                     xvaldata(id)=ABS(DATA(ix+index))
                  ENDIF
                  e2=e2+xvaldata(id)
               ENDDO
               DO id=1,ndep
                  xvaldata(id)= xvaldata(id)/e2
                  xvalresp(id)= xvalresp(id)/e1
               ENDDO
c     
               xmin= 1e10
               xmax=-1e10
               DO i=1,ndep
                  IF (xmin.GT.xvalresp(i)) xmin=xvalresp(i)
                  IF (xmax.LT.xvalresp(i)) xmax=xvalresp(i)
                  IF (xmin.GT.xvaldata(i)) xmin=xvaldata(i)
                  IF (xmax.LT.xvaldata(i)) xmax=xvaldata(i)
               ENDDO
               CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
               WRITE(XTXT,820) NXDIF
 820           FORMAT('Magnitude (10**',I3,')$')
               PTIT='Magnitude'
c     

C     *** WRITE PLP FILE
               CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,Xaxis,Yaxis,
     &              IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &              YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)

c***  WRITE PLM FILe
               jjplot=jjplot+1

               WRITE(13,'(a,i3,a,4f15.6,a)')
     &              ' plaxis(:,',jjplot,')=[',
     &              xleft,xright,ydown,yup,']'';'

c**   WRITE the observed DATA to the plt file
               CALL PLTWRI(ndep,0.,0.,0.,0.,xvaldata,1,yarray,1)

c**   WRITE synthetic DATA to plt file
               CALL PLTWRI(ndep,0.,0.,0.,0.,xvalresp,1,yarray,1)
c     
c-----------Phase----------------
c
c     Xright=MIN(xright,180.)
            ENDIF
c*************************************************************
c     phase plot
c****************************************************c   
            IF (iphase(2).EQ.1) THEN
               WRITE(*,*) 'plotting phase'
               PTIT='REL. PHASE'
               WRITE(XTXT,830) 
 830           FORMAT('Phase (degrees) $')
c---- calculated DATA
               xvalresp(id)=0
               DO id=1,ndep
                  index=(id +(ifrq-1)*ndep-1)*nx
c     xvalresp(id)=
c     &     57.29578*ATAN2(AIMAG(resp(1+index)),REAL(resp(1+index)))
c     &     -xvalresp(1)
                  ch =(resp(ix+index))
c     WRITE(*,*)id,resp(1+index),DATA(1+index) 
                  IF ((REAL(ch).NE.0.) .OR. 
     &                 (AIMAG(ch).NE.0.)) THEN
                     xvalresp(id)=
     &                    57.29578*ATAN2(AIMAG(ch),REAL(ch))
     &                    -xvalresp(1)           
                  ELSE
                     xvalresp(id)=0
                     WRITE(*,*)' ***Resp for  point &freq',
     &                    id,ifrq,' is zero !!!'
                     WRITE(*,*)'data point:',(resp(ix+index))
                  ENDIF
                  e1=e1+xvalresp(id)
               ENDDO
               xvalresp(1)=0
c---  scale DATA
               xvaldata(1)=0
               DO id=1,ndep
                  IF (iopt(13).EQ.1) THEN
                     ch =usvd(id,1)
                  ELSE
                     index=(id +(ifrq-1)*ndep-1)*nx
                     ch =(DATA(1+index))
                  ENDIF
c     WRITE(*,*)id,resp(1+index),DATA(1+index) 
                  IF ((REAL(ch).NE.0.) .OR. (AIMAG(ch).NE.0.)) THEN
                     xvaldata(id)=
     &                    57.29578*ATAN2(AIMAG(ch),REAL(ch))
     &                    -xvaldata(1)           
                  ELSE
                     xvaldata(id)=0
                     WRITE(*,*)' ***Data for  point &freq',
     &                    id,ifrq,' is zero !!!'
                     WRITE(*,*)'data point:',ch
                  ENDIF
                  e2=e2+xvaldata(id)
               ENDDO
               xvaldata(1)=0
               
c     *** unwrap phase
 99            XMAX=MAX(xvalresp(1),xvaldata(1))
               XMin=MIN(xvalresp(1),xvaldata(1))
               idum=0
               DO idep=2,ndep
c     WRITE(*,*)' xval,idep',xvalresp(idep),idep
                  IF ((xvaldata(idep)-xvaldata(idep-1)).GT.180.1) THEN
                     xvaldata(idep)= xvaldata(idep)-360
                     idum=1
                  ENDIF
                  IF ((xvaldata(idep)-xvaldata(idep-1)).LT.-180.1) THEN
                     xvaldata(idep)= xvaldata(idep)+360
                     idum=1
                  ENDIF
                  IF ((xvalresp(idep)-xvalresp(idep-1)).GT.180.1) THEN
                     xvalresp(idep)= xvalresp(idep)-360
                     idum=1
                  ENDIF
                  IF ((xvalresp(idep)-xvalresp(idep-1)).LT.-180.1) THEN
                     xvalresp(idep)= xvalresp(idep)+360
                     idum=1
                  ENDIF
                  XMAX=MAX(xvalresp(idep),xvaldata(idep),xmax)
                  XMin=MIN(xvalresp(idep),xvaldata(idep),xmin)
               ENDDO
               IF (idum.NE.0) GOTO 99 
C eliminate the phase difference between measured and predicted fields
C Sep 16, 2003, cfh
               DO idep = 2,ndep
                  IF ((xvaldata(idep) - xvalresp(idep)).GT.200.0) THEN
                     xvalresp(idep) = xvalresp(idep)+360.0
                  ENDIF
                  IF ((xvaldata(idep) - xvalresp(idep)) .LT.-200.0)THEN
                     xvalresp(idep) = xvalresp(idep)-360.0
                  ENDIF   
               ENDDO


c     XMAX= 180
c     XMIN=-180
               CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
               xdiv=1
               nxdif=0
               

C     *** WRITE PLP FILE
               CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,Xaxis,Yaxis,
     &              IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &              YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)

c***  WRITE PLM FILe
               jjplot=jjplot+1
               WRITE(13,'(a,i3,a,4f15.6,a)')
     &              ' plaxis(:,',jjplot,')=[',
     &              xleft,xright,ydown,yup,']'';'

c**   WRITE the observed DATA to the plt file
               CALL PLTWRI(ndep,0.,0.,0.,0.,xvaldata,1,yarray,1)
               
c**   WRITE synthetic DATA to plt file
               CALL PLTWRI(ndep,0.,0.,0.,0.,xvalresp,1,yarray,1)
c********************************************************
c     END jump loop
            ENDIF
            IF (iphase(3).EQ.1) THEN
c*********************************************************
c     
c-----------Residuals----------------
c     
c     Xright=MIN(xright,180.)
               NC=2
               PTIT='POWER'
               WRITE(XTXT,840) 
 840           FORMAT(' power$')
c---  scale DATA
               asum=0
               e1=0
               e2=0
               DO id=1,ndep
                  index=(id +(ifrq-1)*ndep-1)*nx
                  IF (iopt(13).EQ.1) THEN
                     ch =usvd(id,1)
                  ELSE
                     ch =(DATA(1+index))
                  ENDIF
                  asum=asum+CONJG(ch)*resp(ix+index)
                  e1=e1+CONJG(resp(ix+index))*resp(ix+index)
                  e2=e2+CONJG(ch)*ch
c     WRITE(*,*)'asum,e1,e2',asum,e1,e2
               ENDDO
               e1=SQRT(e1)
               e2=SQRT(e2)
c     WRITE(*,*)'asum,e1,e2',asum,e1,e2
               asum=asum/e1/e2
c---- calculated DATA
               DO id=1,ndep
                  index=(id +(ifrq-1)*ndep-1)*nx
                  IF (iopt(13).EQ.1) THEN
                     ch =usvd(id,1)
                  ELSE
                     ch =(DATA(1+index))
                  ENDIF

                  xvaldata(id)=ch*CONJG(ch)/e2**2
                  xvalresp(id)=
     1                 REAL(CONJG(asum)*CONJG(ch)*resp(ix+index)/e2/e1 )
c     WRITE(*,*)id,xvaldata(id),xvalresp(id)
               ENDDO

               XMAX=-1e10
               XMin=1e10
               DO idep=1,ndep
c     WRITE(*,*)' xval,idep',xvalresp(idep),idep
                  XMAX=MAX(xvalresp(idep),xvaldata(idep),xmax)
                  XMin=MIN(xvalresp(idep),xvaldata(idep),xmin)
               ENDDO
c     WRITE(*,*)'xmin,xmax',xmin,xmax
               CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
               xdiv=1
               nxdif=0


C     *** WRITE PLP FILE
               CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,Xaxis,Yaxis,
     &              IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &              YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)
c***  WRITE PLM FILe
               jjplot=jjplot+1
               WRITE(13,'(a,i3,a,4f15.6,a)')
     &              ' plaxis(:,',jjplot,')=[',
     &              xleft,xright,ydown,yup,']'';'


c**   WRITE the observed DATA to the plt file
               CALL PLTWRI(ndep,0.,0.,0.,0.,xvaldata,1,yarray,1)

c**   WRITE residuals to plt file
               CALL PLTWRI(ndep,0.,0.,0.,0.,xvalresp,1,yarray,1)
c**************************************************
            ENDIF
c*******************************************************
 1000    CONTINUE               ! loop over ncurv

         RETURN
         END

C**************************************************************
      SUBROUTINE PLTpdp(np1,np2,pdp)
C     *** PLOT OF PDP DISTRIBUTIONS     
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './oases/complo.f'
      REAL*8 pdp(mpar,0:(mdig-1))
      INTEGER nlab,igrid,nc,i,j,nx_plotdim
      PARAMETER(nx_plotdim=5000)
      REAL yval(nx_plotdim),xoff
      REAL xdiv,ydiv,xdum,dx,dx0,yoff
      REAL xmin,xmax,ymin,ymax
      REAL  XLEFT,XRIGHT,XAXIS,XINC ! data for the range axis
      REAL YUP,YDOWN,YAXIS,YINC
      REAL x1,x2,dx2,delpdp
      INTEGER j2,ndig,idum
      INTEGER np1,np2,ii
      INTEGER nxdif,nydif
      CHARACTER*6 OPTION(6),OPTTL(3)
      CHARACTER*80 title1
      DATA OPTTL /' PSTCK','WTLRAN','UTLRAN'/
      DATA option /6*'      '/
c---  idum =1 is used for plotting figures using physical numbers
      IF (iopt(14).EQ.1)  THEN
         idum=1
      ELSE
         idum=0
      ENDIF
      xaxis=15
c     yaxis=10
      yaxis=12
      nc=nparm
      IF (idum.EQ.1) THEN
         nc=np2-np1
         yaxis=yaxis*(nc)/nparm
         IF (yaxis.GT.1.3333*nc) yaxis=1.3333*nc  
      ENDIF
      OPTION(1)=PROGNA
      OPTION(2)=OPTTL(1)
      option(3)=',UNI,V'
      option(3)=',   ,V'
      IF (idum.EQ.1) THEN
         option(4)='TT,SEG'
         IF (iopt(19).EQ.1) THEN
            option(3)=',UNI,A'
            option(4)='DD,SEG'
         ENDIF
      ELSE
c     option(4)='TT    '
         option(4)='TT,ADD'
      ENDIF
c     IF (DEBUG) WRITE(prtfil,*) 'ENTERING PLTPDP'

c---  definition of xaxis 
      PTIT=' PDP-DISTRIBUTION'
      XMAX=1
      XMIN=0
      CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
      XDIV=1
      xinc=0.2
      XTXT=' range of parameter$'
      IF (idum.EQ.1) THEN
         xmax=-10e10
         xmin=10e10
         DO i=np1,np2-1
            xmax=MAX(xmax,fmax(i))
            xmin=MIN(xmin,fmin(i))
         ENDDO
         CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
         xright=xmax
         xleft =xmin
         XDIV=1
         xinc=xinc/2
         xtxt=phystxt(par2phy(np1))
      ENDIF
c     PTIT='    PDP'
      PTIT='     '
c     WRITE(XTXT,820) NXDIF
c     820  FORMAT('range (10**',I3,')$')
c     NXDIF=1

c     I=MAX(INR,1)
      NLAB=0
c     WRITE(LAB(1),810) FRQ
c     810  FORMAT('Freq:',F7.1,' Hz$')
      XTYP='LIN'
      YTYP='LIN'
      IGRID=0
c     
c---  scale yaxis
c     
      ymin=nc
c     ymin=nc+1
      ymax=0
c     IF (idum.EQ.1) ymin=nc
c     WRITE(prtfil,*)'plotga:nc,ymin',nc,ymin
      CALL AUTOAX(YMIN,YMAX,Ydown,YUP,YINC,YDIV,NYDIF)
      YTXT='Distribution Number$' 
      YTXT=' $'
c     YTXT='Layer $' 
      ydiv=1
      ydown=ymin
      yinc=-1
C     *** WRITE PLP FILE
c     WRITE(*,*)'idum,par2phy(np1),np1',idum,par2phy(np1),np1
      IF (idum.EQ.1) THEN
         TITLE1='                                '
         IF (par2phy(np1).EQ.9) THEN
            xLEFT=xLEFT/1000
            xright=xright/1000
            xinc=xinc/1000
         ENDIF
      ELSE
         title1=title
      ENDIF
c     WRITE(prtfil,*)'plotga:nc,ydown',nc,ydown
      CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,Xaxis,Yaxis,
c     CALL PLPWRI(OPTION,PTIT,TITLE1,NLAB,LAB,Xaxis,Yaxis,
     &     IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &     YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)

c---  find maximum digits
      ndig=0
      DO i=1,nc
         ii=i-1+np1
         ndig=MAX(ndig,ndigit(ii))   
c     WRITE(*,*)'ndig,ndigit', ndig,ndigit(ii)  
      ENDDO
      dx=1./(ndig-1)
      IF (idum.EQ.1) THEN
         dx=10e10
         DO i=1,nc
            ii=i-1+np1
            dx=MIN(dx,(fmax(ii)-fmin(ii))/(ndigit(ii)-1) )
         ENDDO
         ndig=(xmax-xmin)/dx + 1.5
      ENDIF
c     
c**   WRITE all the ppd's to file
c     
      do i=1,nc
         ii=i-1+np1
         IF (ndigit(ii).gt.nx_plotdim) stop ' DIMENSION in pltpdp'
         xoff=0
         dx0=dx
         if (idum.eq.1) then
            xoff=fmin(ii)
c     dx0=dx0*(fmax(ii)-fmin(ii))
            dx0=dx
         endif
c     write(*,*)dx
c**** for scaling of each PDP to one the following line should be used
c     fac=  1.0/real(pdp(ii,maxpdp(ii)))
         yoff=i
         delpdp=0.
         if (idum.eq.1) then
            dx2=1./(ndigit(ii)-1)*(fmax(ii)-fmin(ii)) !increase for resampling
            x2=fmin(ii)         ! next observed sampling point
            x2=0
         else
            dx2=1./(ndigit(ii)-1) !increase for resampling
            x2=xmin             ! next observed sampling point
            x2=0
         endif
         x1=xmin                ! current resampling point
         x1=0
         j2=0                   ! left sampling point
         if (x1.eq.x2) then
            x2=x2+dx2 
            j2=j2+1
            yval(1)=pdp(ii,1-1)
         else 
            yval(1)=0.
         endif
         do j=2,ndig
            x1=x1+dx
c     write(*,*)x1,x2
            if(x1.ge.x2) then
               if ( (x1.gt.((fmax(ii)-fmin(ii))*1.000002))
     1              .and.(idum.eq.1)) then
                  write(*,*)'x1,x2',x1,x2
                  write(*,*)'gt. fmax'
                  yval(j)=0
                  delpdp =0
               else
                  j2=j2+1
                  delpdp  =(pdp(ii,j2+1-1)-pdp(ii,j2-1))/dx2
                  yval(j) = pdp(ii,j2-1)+ delpdp*(x1-x2)
                  x2=x2+dx2
                  delpdp=delpdp*dx
               endif
            else
               yval(j)=yval(j-1)+delpdp
            endif
c     yval(j)=pdp(ii,j-1)*fac        ! scaled so maximum of wiggle is one
c     yval(j)=1.*pdp(ii,j-1)/real(pdp(ii,maxpdp(ii)))
         enddo
         if (idum.eq.1) then
            if (par2phy(np1).eq.9) then
               xoff=xoff/1000
               dx0=dx0/1000
            endif
         endif

c     PLTWRI(N,XOFF,DX,YOFF,DY,X,IX,Y,IY)
         CALL PLTWRI(ndig,xoff,dx0,yoff,0.,xdum,1,yval,1)
      enddo

      if (iopt(14).eq.1) then
         write(13,'(a,i2,a,i10,a)')'   npdpsamp(',ii,')=',ndig,';'
      else
         write(13,'(a,i4,a)')' npdpsamp=',ndig   ,';'
      endif
      if (((iopt(14).eq.1) .and. (ii.eq.1)).or.(iopt(14).eq.0)) then
         write(13,'(a,i4,a)')' nparm='  , nparm ,';'
      endif

      if ((iopt(14).eq.1)) then
         if (par2phy(np1).eq.9) then
            write(13,810)' f_min(',ii,')=',fmin(ii)/1000,';'
            write(13,810)' f_max(',ii,')=',fmax(ii)/1000,';'
         else
            write(13,810)' f_min(',ii,')=',fmin(ii),     ';'
            write(13,810)' f_max(',ii,')=',fmax(ii),     ';'
         endif
         write(13,'(a,i2,a,i3,a)')   'par2phy(',ii,')=',par2phy(ii),';'
      else
         write(13,'(a,i2,a,i3,a)')   'par2phy(',1,')=',par2phy(1),';'
         write(13,810)' f_min(',1,')=',fmin(1),     ';'
         write(13,810)' f_max(',1,')=',fmax(1),     ';'
         
      endif
C     *** FORMATS

 810  FORMAT(a,i2,a,e15.6,a)
 811  FORMAT('SD:',F9.1,' m$')
 812  FORMAT('RD:',F9.1,' m$')
      RETURN
      END                       ! 

      SUBROUTINE scalepdp(pdp,maxpdp)
C *** PLOT OF PDP DISTRIBUTIONS     
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER maxpdp(mpar)  ! best  model in population
      REAL*8 pdp(mpar,0:(mdig-1))
      INTEGER nc,i,j
      real fac
      nc=nparm
c--- normalize pdp to unit distribution in interval 0-1
c--- this option should not be used for equal plotting 11/12 92
      do i=1,nc
        do j=0,ndigit(i)-1
          pdp(i,j)=pdp(i,j)*ndigit(i)
        enddo
      enddo

c--- find scaling factor
      fac=0.
      do i=1,nc
c        if (i.eq.4) goto 99  ! when should that be used ???
        fac=max(fac,real(pdp(i,maxpdp(i))) )
99      continue
c        if (fac.gt.real(pdp(i,maxpdp(i))) ) then 
c           imaxpdp=i
c           fac=real(pdp(i,maxpdp(i)))  
c        endif
      enddo
      fac=1./fac
c--- normalize pdp to unit distribution in interval 0-1
      do i=1,nc
        do j=0,ndigit(i)-1
          pdp(i,j)=pdp(i,j)*fac
          if (pdp(i,j).gt.1)  pdp(i,j)=1.0 
        enddo
      enddo
      end
c*******************************************************************
      SUBROUTINE PLTmod(nc,np1,np2,xstart)
C     *** PLOT OF PDP DISTRIBUTIONS     
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './oases/complo.f'
      INTEGER nlab,igrid,nc,nc1,i,j,idum,nx_plotdim
      parameter(nx_plotdim=5000)
      REAL xval(nx_plotdim),yval(nx_plotdim),xstart(*)
      REAL xdiv,ydiv,dx,yoff
      REAL xmin,xmax,ymin,ymax
      REAL  XLEFT,XRIGHT,XAXIS,XINC ! data for the range axis
      REAL YUP,YDOWN,YAXIS,YINC
      INTEGER nxdif,nydif,istart
      integer np1,np2,np3,jj
      CHARACTER*6 OPTION(6),OPTTL(3)
      DATA OPTTL /' mpdp','WTLRAN','UTLRAN'/
      data option /6*'      '/
      xaxis=15
      yaxis=10
      yaxis=12
c---  idum =1 is used for plotting figures using physical numbers
c     idum=1
      if ((iopt(14).eq.1))  then
         idum=1
      else
         idum=0
      endif
c     idum=0
      istart=0
      nc1=abs(nc)
      if (nc.lt.0) then
         istart=1
c     nc=1
      endif
      np3=nparm
      if (idum.eq.1) then
         np3=np2-np1
         yaxis=yaxis*(np3+1)/nparm
      endif
      OPTION(1)=PROGNA
      OPTION(2)=OPTTL(1)
      option(3)=  ',UNI,V'
      option(3)=  ',   ,V'
      option(4)=  'TT,    '
      if (istart.eq.1) then
         OPTION(1)='      '
         OPTION(2)='      '
         OPTION(3)=',ADD  '
         OPTION(4)='      '
      endif
c---  definition of xaxis 
      PTIT=' BEST MODELS'
      XMAX=1
      XMIN=0     
      CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
      XDIV=1
      xinc=0.2
      XTXT=' range of parameter$'
      if (idum.eq.1) then
         xmax=0
         xmin=10e10
         do i=np1,np2-1
            xmax=max(xmax,fmax(i))
            xmin=min(xmin,fmin(i))
         enddo
         CALL AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
         xright=xmax
         xleft =xmin
         XDIV=1
         xtxt=phystxt(par2phy(np1))
      endif
      PTIT='    BEST MODELS'
c     WRITE(XTXT,820) NXDIF
c     820  FORMAT('RANGE (10**',I3,')$')
      NXDIF=1

c     I=MAX(INR,1)
      NLAB=0
c     WRITE(LAB(1),810) FRQ
c     810  FORMAT('Freq:',F7.1,' Hz$')
      XTYP='LIN'
      YTYP='LIN'
      IGRID=0
c     
c---  scale yaxis
c     
      ymin=np3+1
      ymax=0.
      if (istart.eq.1)ymin=np3
c     if (idum.eq.1) ymax=1
      CALL AUTOAX(YMIN,YMAX,Ydown,YUP,YINC,YDIV,NYDIF)
      YTXT='Distribution Number $'
      ydiv=1
      ydown=ymin
      yinc=-1

C     *** WRITE PLP FILE
      CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,Xaxis,Yaxis,
     &     IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &     YDOWN,YUP,YINC,YDIV,YTXT,YTYP,nc1)

c     write(*,*)' after plpwri'
c**   write all the ppd's to file
      IF (istart.NE.1) THEN
         DO i=1,nc
c     WRITE(*,*)' dddd',nc,dx,yoff,np1
            dx=0.
c     WRITE(*,*)' dum',nc,dx,yoff,np3
            yoff=0.
c     WRITE(*,*)np3,np1
            DO j=1,np3
               jj=j-1+np1
c     WRITE(*,*)'jj,j,i,np3,np1',jj,j,i,np3,np1
               xval(2*j-1)=1.*model(jj,i)/ndigit(jj) !
c     WRITE(*,*)'basta' 
               IF (idum.EQ.1) THEN
                  xval(2*j-1)= fmin(jj)+xval(2*j-1)*(fmax(jj)-fmin(jj))
               ENDIF        
               xval(2*j  ) = xval(2*j-1)
c     WRITE(*,*) model(j,i),ndigit(j),1.*model(j,i)/ndigit(j)
               yval(2*j-1) =j-0.5
               yval(2*j  ) =j+0.5
            ENDDO
c     PLTWRI(N,XOFF,DX,YOFF,DY,X,IX,Y,IY)
            CALL PLTWRI(2*np3,0.,dx,yoff,0.,xval,1,yval,1)
         ENDDO
      ELSE
         DO i=1,nc1
            dx=0.
            yoff=0
            DO j=1,np3
               jj=j-1+np1
               xval(2*j-1)=(XSTART(jj)) !-FMIN(jJ))/(FMAX(jJ)-FMIN(jJ))
c     pg 13/6 97          xval(2*j-1)=(XSTART(jj)-FMIN(jJ))/(FMAX(jJ)-FMIN(jJ))
               IF (idum.EQ.1) THEN
                  xval(2*j-1)=XSTART(jj)
                  IF (par2phy(jj).EQ.9) THEN
                     xval(2*j-1)=XSTART(jj)/1000
                  ENDIF
               ENDIF        
c     c           WRITE(*,*)'xval',xval(2*j-1),j
               xval(2*j  )= xval(2*j-1)
               yval(2*j-1)=j-0.5
               yval(2*j  )=j+0.5
            ENDDO
c     PLTWRI(N,XOFF,DX,YOFF,DY,X,IX,Y,IY)
            yval(1)=0.
            IF (istart.EQ.1) yval(2*np3)=np3
            CALL PLTWRI(2*np3,0.,dx,yoff,0.,xval,1,yval,1)
         ENDDO
         
      ENDIF

      RETURN
      END

