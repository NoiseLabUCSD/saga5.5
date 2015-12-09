c******************************************************
      SUBROUTINE AUTOAX(XMIN,XMAX,XLEFT,XRIGHT,XINC,XDIV,NXDIF)
C *** DETERMINE FACTOR
      xhelp=abs(XMAX-XMIN)
      ILOG=IFIX(ALOG10(xhelp))
      IF (xhelp.LT.1.0) ILOG=ILOG-1
      IFAC=IFIX((Xhelp)/10.**ILOG)+1
      NXDIF=ILOG-1
      XDIV=10.**(-NXDIF)
C *** MAKE NICE AXIS
      XFAC=IFAC
c      IF (IFAC.EQ.1) XFAC=5
      IF (IFAC.EQ.1) XFAC=5
c      IF (IFAC.EQ.2) XFAC=4
      IF (IFAC.EQ.2) XFAC=10
      IF (IFAC.GT.5) XFAC=IFAC/2.0
      XINC=IFAC*10.**ILOG/XFAC
      IF (XMIN.GE.0E0) THEN
       XLEFT=IFIX(XMIN/XINC+0.01)*XINC
      ELSE
       XLEFT=-IFIX(-XMIN/XINC+0.99)*XINC
      END IF
      IF (XMAX.GE.0E0) THEN
       XRIGHT=IFIX(XMAX/XINC+0.99)*XINC
      ELSE
       XRIGHT=-IFIX(-XMAX/XINC+0.01)*XINC
      END IF
c---  for negative values svap xmin and xmax
      if (xmax.lt.xmin) then
c        xhelp=xleft
c        xleft=xright
c        xright=xhelp
        xinc=-xinc
      endif
      RETURN
      END
c*************************************************************
      SUBROUTINE PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,XLEN,YLEN,
     &                  IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &                  YLO,YUP,YINC,YDIV,YTXT,YTYP,NC)
c   pg oct 92: changed from option(2) to allow writing options out 
      CHARACTER*6 OPTION(6)
c
      CHARACTER*16 LAB(20)
      CHARACTER*80 TITLE,PTIT,XTXT,YTXT
      CHARACTER*3 XTYP,YTYP
C
C     WRITES GENERAL PLP FILE
C
      WRITE(19,777) OPTION
      WRITE(19,778) PTIT
      WRITE(19,778) TITLE
      WRITE(19,6010) NLAB,'NUMBER OF LABELS'
      DO 10 ILAB=1,NLAB
 10   WRITE(19,779) LAB(ILAB)
      WRITE(19,6030) XLEN,'XLEN'
      WRITE(19,6030) YLEN,'YLEN'
      WRITE(19,6010) IGRID,'GRID TYPE. 0: NO GRID'
      WRITE(19,6030) XLEFT,'XLEFT'
      WRITE(19,6030) XRIGHT,'XRIGHT'
      WRITE(19,6030) XINC,'XINC'
      WRITE(19,6030) XDIV,'XDIV'
      WRITE(19,778) XTXT
      WRITE(19,780) XTYP
      WRITE(19,6030) YLO,'YDOWN'
      WRITE(19,6030) YUP,'YUP'
      WRITE(19,6030) YINC,'YINC'
      WRITE(19,6030) YDIV,'YDIV'
      WRITE(19,778) YTXT
      WRITE(19,780) YTYP
      WRITE(19,6010) NC,'NC'
C *** FORMATS
 777  FORMAT(1H ,6A6)
 778  FORMAT(1H ,A80)
 779  FORMAT(1H ,A16)
 780  FORMAT(1H ,A3)
 6010 FORMAT(1H ,I8,10X,A40)
 6020 FORMAT(1H ,F15.6,3X,A40)
 6030 FORMAT(1H ,G15.6,3X,A40)
      RETURN
      END
      SUBROUTINE PLTWRI(N,XOFF,DX,YOFF,DY,X,IX,Y,IY)
      DIMENSION X(1),Y(1)
      WRITE(19,6010) N,'N'
      WRITE(19,6030) XOFF,'XOFF'
      WRITE(19,6030) DX,'DX'
      WRITE(19,6030) YOFF,'YOFF'
      WRITE(19,6030) DY,'DY'
      IF (DX.EQ.0E0) WRITE(20,444)(X(II),II=1,N,IX)                    
      IF (DY.EQ.0E0) WRITE(20,444)(Y(II),II=1,N,IY)                    
 444  FORMAT(1H ,6G16.8)
 6010 FORMAT(1H ,I8,10X,A40)
 6030 FORMAT(1H ,G15.6,3X,A40)
      RETURN
      END
c****************************************
c condrw1 is a sligthly modifyed version of condrw
       SUBROUTINE CONDRW1(TITLE,NPX,NPY,NX,NY,XLEFT,XRIGHT   
     $,XSCALE,XINC,YUP,YDOWN,YSCALE,YINC,ZMIN              
     $,ZMAX,ZSTEP,FREQ,SD,RECUP,RECLO,X1,XL,PX,titlex,titley,conopt)        
      real XLEFT,XRIGHT,XSCALE,XINC,ydown,yup,ySCALE,yINC,
     &     dx,dy,ypointup,ypointdown,xpointleft,xpointright
       DIMENSION SECTOR(28),PX(1)      
      CHARACTER*50 FILENM,conopt
      CHARACTER*4 TITLE(20)
      character*40 TITLEX,TITLEY      
      DATA X1PL,Y1PL/1.4,1.4/,HGTPT,HGTC,LABPT,NDIV,        
     *NARC/0.1,0.14,-3,1,5/,LABC,LWGT/-1,-1/,NSM/0/         
      DATA DUMMY /0./
      
      save nplot_c
      if ((nplot_c.le.0.).or.(nplot_c.ge.100.)) then
         nplot_c=1
      else
         nplot_c= nplot_c+1
      endif
      DO  I=2,40,1
         IF((titlex(I:I).EQ.'$') .or.
     1       (titlex(I:I+1).EQ.'  ')) then
            ncx=I-1
            GO TO 8
         endif
      enddo
    8 CONTINUE
      DO  I=2,40,1
         IF((titley(I:I).EQ.'$') .or.
     1       (titley(I:I+1).EQ.'  ')) then
            ncy=I-1
c     write(*,*)phystxt2(j)
            GO TO 18
         endif
      enddo
 18   CONTINUE
 
c---  write to matlab file
      write(13,'(a,i2,a)')' nplot_c=',nplot_c,';'
      write(13,'(a,i2,a,i5,a)')' npx(',nplot_c,')=',nx,';'
      write(13,'(a,i2,a,i5,a)')' npy(',nplot_c,')=',ny,';'
      write(13,'(a,i2,a,f15.6,a,f15.6,a,i4,a)')
     1   ' yp(',nplot_c,',:)=''linspace(',reclo,',',recup,',',ny,')'';'
      write(13,'(a,i2,a,f15.6,a,f15.6,a,i4,a)')
     1   ' xp(',nplot_c,',:)=''linspace(',x1   ,',',xl   ,',',nx,')'';'
      write(13,'(a,i2,a,i5,a)')' ncx(',nplot_c,')=',ncx,';'
      write(13,'(a,i2,a,a,a)')
     1     'title_x(',nplot_c,',:)=''',titlex,''';'
      write(13,'(a,i2,a,i5,a)')' ncy(',nplot_c,')=',ncy,';'
      write(13,'(a,i2,a,a,a)')
     1     'title_y(',nplot_c,',:)=''',titley,''';'
c                
c   formats      
 401  FORMAT(1H ,F15.4,3X,'  NUMBER OF DATA POINTS ALONG THE X AXIS')         
 402  FORMAT(1H ,F15.4,3X,'  NUMBER OF DATA POINTS ALONG THE Y AXIS')         
 403  FORMAT(1H ,F15.4,3X,'  DIVX ' )
 404  FORMAT(1H ,F15.4,3X,'  DIVY ' )
 405  FORMAT(1H ,F15.4,3X,'  FLAGRC ' )
 406  FORMAT(1H ,F15.4,3X,'  RDUP ' )
 407  FORMAT(1H ,F15.4,3X,'  RDLO ' )           
 408  FORMAT(1H ,F15.4,3X,'  SOURCE DEPTH (M) ' )          
 409  FORMAT(1H ,F15.4,3X,'  NUMBER OF GRID POINTS ALONG THE X AXIS ' )        
 410  FORMAT(1H ,F15.4,3X,'  NUMBER OF GRID POINTS ALONG THE Y AXIS ' )        
 411  FORMAT(1H ,F15.4,3X,'  FREQUENCY (HZ)' )             
 412  FORMAT(1H ,F15.4,3X,'  DUMMY ' )
 413  FORMAT(1H ,F15.4,3X,'  CAY ' )
 414  FORMAT(1H ,F15.4,3X,'  NRNG ' )
 415  FORMAT(1H ,F15.4,3X,'  ZMIN ' ) 
 416  FORMAT(1H ,F15.4,3X,'  ZMAX ' ) 
 417  FORMAT(1H ,F15.4,3X,'  ZINC ' ) 
 418  FORMAT(1H ,F15.4,3X,'  X ORIGIN OF PLOT IN INCHES ' )
 419  FORMAT(1H ,F15.4,3X,'  DUMMY ' )
 420  FORMAT(1H ,F15.4,3X,'  Y ORIGIN OF PLOT IN INCHES ' )
 421  FORMAT(1H ,F15.4,3X,'  NSM   ' )
 422  FORMAT(1H ,F15.4,3X,'  HGTPT ' )
 423  FORMAT(1H ,F15.4,3X,'  HGTC ' ) 
 424  FORMAT(1H ,F15.4,3X,'  LABPT ' )
 425  FORMAT(1H ,F15.4,3X,'  NDIV ' ) 
 426  FORMAT(1H ,F15.4,3X,'  NARC ' ) 
 427  FORMAT(1H ,F15.4,3X,'  LABC ' ) 
 428  FORMAT(1H ,F15.4,3X,'  LWGT ' ) 
 800  FORMAT('CONDR,FIP,FMT,CPX,',a)             
 801  FORMAT(A50)                
 850  FORMAT(a40)                   !(20A4)                    
 851  FORMAT(20A4)                    
 900  FORMAT(1X,F15.4,3X,'  XLEFT',/,F15.4,4X,'  XRIGHT',/,F15.4,3X,           
     *'   XSCALE',/,F15.4,4X,'  XINC')
  901 FORMAT(1X,F15.4,3X,'  YUP',/,F15.4,4X,'  YDOWN',/,F15.4,3X,              
     *'   YSCALE',/,F15.4,4X,'  YINC')
  950 FORMAT(1H ,F15.4,1X,'    RMIN',/,F15.4,2X,'    RMAX')
c      write(*,*)'conopt',conopt
      WRITE(28,800)conopt             
      WRITE(28,851)TITLE              
c      CALL VCLR(SECTOR,1,28)
      do i=1,28
         sector(i)=0.
      enddo
      SECTOR(1)=NPX                   
c                
c   sector(4) is a flag which is set to zero in the range  
c   dependent version of snap for all sectors except the last                  
c   one. here is used to indicate that this is the last sector                
      SECTOR(4)=1.0
                               
      write(29,*)
      WRITE(29,444) (SECTOR(L),L=1,28)
      write(29,*)
 444  FORMAT(1H ,6G13.5)
      INQUIRE(UNIT=29,NAME=FILENM)
      WRITE(28,801) FILENM             
      DIVX=1E0
      DIVY=1E0
      CAY=5.
      NRNG=5
      FLAGRC=0.
      WRITE(28,850)TITLEX             
      R1=X1
      R2=XL
      WRITE(28,950)X1,xl             
      AX1=XLEFT
      AX2=XRIGHT
      AX3=XSCALE
      AX4=XINC
      WRITE(28,900)xleft,xright,xscale,xinc
c      WRITE(28,900)AX1,AX2,AX3,AX4    
      WRITE(28,850)TITLEY             
      WRITE(28,901)YUP,YDOWN,YSCALE,YINC                   
      WRITE(28,401)FLOAT(NPX)
      WRITE(28,402)FLOAT(NPY)
      WRITE(28,403)DIVX              
      WRITE(28,404)DIVY              
      WRITE(28,405)FLAGRC              
      WRITE(28,406)RECUP              
      WRITE(28,407)RECLO                   
      WRITE(28,408)SD     
c   number of grid points along the x axis                 
      WRITE(28,409)FLOAT(NX)          
c   number of grid points along the y axis                 
      WRITE(28,410)FLOAT(NY)          
      WRITE(28,411)FREQ      
      WRITE(28,412)DUMMY              
      WRITE(28,413)CAY              
      WRITE(28,414)FLOAT(NRNG)              
      WRITE(28,415)ZMIN               
      WRITE(28,416)ZMAX               
      WRITE(28,417)ZSTEP              
c x origin  of plot in inches         
      WRITE(28,418)X1PL               
      WRITE(28,419)DUMMY              
c y origin  of plot in inches         
      WRITE(28,420)Y1PL               
      WRITE(28,421)FLOAT(NSM)                
      WRITE(28,422)HGTPT              
      WRITE(28,423)HGTC               
      WRITE(28,424)FLOAT(LABPT)       
      WRITE(28,425)FLOAT(NDIV)        
      WRITE(28,426)FLOAT(NARC)        
      WRITE(28,427)FLOAT(LABC)        
      WRITE(28,428)FLOAT(LWGT)        
      RETURN     
      END       
       SUBROUTINE CONDRB(NP1,NP2,NPX,PX)                   
       DIMENSION PX(1)                
       DO 1000 I=NP1,NP2              
       PX(I-NP1+1)=PX(I)              
1000   CONTINUE  
       WRITE(29,444) (PX(L),L=1,NPX)
 444  FORMAT(1H ,6G13.5)
       RETURN    
       END       
C
