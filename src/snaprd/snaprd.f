      subroutine snaprd(nsect,frq,maxmod0,minmod0,R_slngth, 
     &      R_R0, R_R1, R_R2, R_BETA, R_SCATT, R_C2S, R_C2,
     &      R_C0, R_Z0, R_C1, R_Z1, 
     &      R_CC0, R_CC1, R_H0, R_H1, R_ND0, R_ND1,
     &      msec,maxdep,pressr,nptot,xranges)
c
      complex pressr(*)
      real xranges(*)
      integer nsect
      real R_slngth(msec)
  
      real R_R0(msec), R_R1(msec), R_R2(msec)
      real R_BETA(-1:3,msec), R_SCATT(2,msec), R_C2S(msec) 
      real R_C2(msec)
      DOUBLE PRECISION R_C0(maxDEP,msec), R_Z0(maxDEP,msec), 
     1     R_C1(maxDEP,msec), R_Z1(maxDEP,msec)
      DOUBLE PRECISION R_CC0(msec), R_CC1(msec),cminsnap,
     1     R_H0(msec), R_H1(msec)
      integer R_ND0(msec), R_ND1(msec)
      DOUBLE PRECISION  CC0, CC1
C
C   NDEP IS NUMBER OF DEPTHS ALLOWED IN SAMPLING 
C     THE SOUND SPEED PROFILE NDEP IS ALSO DEFINED IN  MODE
      PARAMETER ( NDEP=201) 
C
C   MSP IS NUMBER OF POINTS USED TO SAMPLE THE DEPTH
C     FUNCTIONS. IT IS ALSO DEFINED IN SUB TLDEP 
      COMMON/FLAGPULSE/FLAGPU
      PARAMETER ( MSP=301)

      PARAMETER (ICF=6,MODEN=200,NRANGE=2001,KSRD=500)
c
C   FLAG IS THE MATRIX THAT CONTAINS ALL FLAG INFORMATION
      COMMON/H/FLAG(ICF)
C
C   RD  CONTAINS THE RECEIVER DEPTH 
      COMMON/M/RD(KSRD),SD,nrd
      COMMON /I/ SECD(3)

      CHARACTER*3  MODPLT

      DIMENSION AVROLD(MODEN),AVRNEW(MODEN)
      DIMENSION XTs(moden,MSP),Xold(moden,MSP),xsou(moden,msp),SS(MSP)
      DIMENSION RNG(NRANGE)
      DIMENSION US(MODEN),URNEW(2,MODEN),UROLD(2,MODEN)

      DOUBLE PRECISION HH0,HH1,FD,T3,T2,T1,D1,
     &  F2,DELTAF,DL2
      DOUBLE PRECISION EKFST,EKOLD,EKNEW,EKAVR
      COMMON /A/ C0(NDEP),C1(NDEP),Z0(NDEP),Z1(NDEP),ND0,ND1
      COMMON /AB/ BETA(3),SCATT(2),C2S
      COMMON /B/ ALFNEW(MODEN)
      COMMON /C/ EKNEW(MODEN)
      COMMON /CFIELD/ FLDRD(3), PULRD(3)
      COMMON /D/ EKFST(MODEN), EKOLD(MODEN), EKAVR(MODEN),
     &           ALFOLD(MODEN), ALFAVR(MODEN)
      COMMON /E/ HH0,HH1,FD,T3,T2,T1,D1
      COMMON /EXPMAX/ TRESH
      COMMON /F/ NTN,LI0,N1,N1P1,N1M2,NTNM2,INC
      COMMON /FLAGS/ PLANE, EK0, SQEK0
      COMMON /FFP/ FFP,DLRAN,P0
      COMMON /G/ R0,R1,R2,C2,H0,H1,FACT0,FACT1,C11
      COMMON /LUNIT/ MSOURC, MODOLD, MODNEW, LUPRT
      COMMON /N/ MINMOD, MAXMOD, HORIG, NFREQ
      COMMON /NRNGD/ HNEW,HOLD,RSTART,REND,SLNGTH,IND1
c
      integer xplane, irflagg
      logical out_plane
      real arrayshape(10)
      common /arrayparm/out_plane,arrayshape,xplane,irflagg

      common /R_rec_sou/rng

C
      DATA HGT/0.14/
      DATA MODPLT/'   '/
  
      DATA PLANE / 0.0 /
      DATA FLAGP/-1.0/
C
      DATA fact0,fact1/1,1/
 240  FORMAT(1H ,A5,1X,4F8.0)
 241  FORMAT(1H ,A5,1X,2F8.0)
 242  FORMAT(1H ,A5)
 250  FORMAT(A80)
C
C   DISTINGUISHING BETWEEN FPS (ENTIRELY DOUBLE PRECISION) AND VAX
c     IEXP=150
      IEXP=17
      TRESH=10.0**IEXP
C     
      JUMP=0
      nfreq=1                   ! 
      DO 4400 IFREQ=1,NFREQ
         minmod=minmod0
         maxmod=maxmod0
         isect=1
         ND0=R_ND0(isect)
         r1=R_R1(isect)
         r2=R_R2(isect)
         h0= R_H0(isect)
         h1=R_H1(isect)
         nd0=R_ND0(isect)
         nd1=R_ND1(isect)
         slngth=R_SLNGTH(isect)
         c2=R_C2(isect)
         c2s=R_C2S(isect)
         do i=1,3
            beta(I)=R_beta(i,isect)
         enddo
         scatt(1)=R_SCATT(1,isect)
         scatt(2)=R_SCATT(2,isect)
         cc0=R_CC0(isect)
         cc1=R_CC1(isect)
         DO 1450 I=1,ND0
            c0(I) = R_C0(I,isect)
            z0(I) = R_z0(I,isect)/h0
 1450    CONTINUE
         if (h1.gt.0) then
            DO I=1,ND1
               c1(I) = R_C1(I,isect)
               z1(I) = (H0+R_Z1(I,isect))/H0
            enddo
         endif
         
         F2= dble(frq)
c     IF(MAXMOD.EQ.0)   GO TO 1500
         IF(MAXMOD.EQ.0)   stop 'maxmod is 0 in snaprd'
         NMES=4
c     pg 241103       deltaf is frequency step; here put to zero
         deltaf=0
         CALL MODE(1,IFREQ,F2,MODQTY,*4400,XTs,MSP,MODPLT,
     *        1.0E-3,JUMP,0,AVRNEW,ALFNEW,EKNEW,MODEN,IPARAM,NMES,
     *        DELTAF,
     $        SS,FLAGP,cc0,cc1)
C     
 1400    CONTINUE
         HORIG=HNEW
         HOLD=HNEW
         RSTART=0.
         REND=SLNGTH
         CALL FILOLD(MODQTY,MSP,XTs,EKFST,EKNEW,EKOLD,EKAVR,
     $        ALFNEW,ALFOLD,ALFAVR,AVROLD,AVRNEW,MODEN,xsou,xold)
         
         nr_start=1
         DO 4300 ISECT=1,NSECT
c     write(*,*)' sector loop ',isect
c---  copy one sector
            
            ISP1=ISECT+1
            
            IF(ISECT.NE.NSECT) THEN
               ND0=R_ND0(isp1)
               r1=R_R1(isp1)
               r2=R_R2(isp1)
               h0= R_H0(isp1)
               h1=R_H1(isp1)
               nd0=R_ND0(isp1)
               nd1=R_ND1(isp1)
               slngth=R_SLNGTH(isp1)
               c2=R_C2(isp1)
               c2s=R_C2S(isp1)
               do i=1,3
                  beta(I)=R_beta(i,isp1)
               enddo
               scatt(1)=R_SCATT(1,isp1)
               scatt(2)=R_SCATT(2,isp1)
               cc0=R_CC0(isp1)
               cc1=R_CC1(isp1)
               DO I=1,ND0
                  c0(I) = R_C0(I,isp1)
                  z0(I) = R_z0(I,isp1)/h0
               enddo
               if (h1.gt.0) then
                  DO I=1,ND1
                     c1(I) = R_C1(I,isp1)
                     z1(I) = (H0+R_Z1(I,isp1))/H0
                  enddo
               endif
               CALL MODE(ISP1,IFREQ,F2,MODQTY,*4400,
     &              XTs,MSP,MODPLT,1.0E-3,JUMP,0,AVRNEW,ALFNEW,
     &              EKNEW,
     &              MODEN,IPARAM,NMES,DELTAF,SS,FLAGP,cc0,cc1)
            ENDIF
 1500       CONTINUE
             
      IF  (out_plane) THEN
c          write(*,*)'out of plane'
             CALL TLRANPG(NPNEW,*4100,Frq,MSP,RNG,
     &           US,UROLD,URNEW,XTs,Xold,xsou,ISP1,NSECT,pressr,nptot,
     &           nr_start,xranges,isect)

      ELSE          
            CALL RANGERD(*4100,ISECT,NSECT,NPNEW,ICF,RNG,NRANGE,SECD)
c     write(*,*)'npnew',npnew
            CALL TLRAN(NPNEW,*4100,Frq,MSP,RNG,
     &           US,UROLD,URNEW,XTs,Xold,xsou,ISP1,NSECT,pressr,nptot,
     &           nr_start,xranges)
      ENDIF
 4100       CONTINUE
            
            IF(ISECT.NE.NSECT) THEN
               CALL EXCHNG(*4300,EKFST,EKNEW,EKOLD,
     $              EKAVR,ALFNEW,ALFOLD,ALFAVR,AVROLD,AVRNEW,MODEN,msp,
     $              xts,xold)
            END IF
            
 4300    CONTINUE

 4400 CONTINUE
C     
      END

c
c******************************************
c
       subroutine rotate(x1p,x2p,x1,x2,thta)

       real x1p,x2p,x1,x2,thta(2)
c       real tht,cs,sn
     
c       tht=thta*3.14159265/180
c       cs =cos(tht)
c       sn =sin(tht)
c       cs=thta(1)
c       sn=thta(2)

       x1p= x1*thta(1)-x2*thta(2)
       x2p= x1*thta(2)+x2*thta(1)

       end
c
c*******************************
c
      subroutine geom(ranout,zdepout,ndep)
      real rngref
      real     CORREC
       integer flagpu
      COMMON /FLAGPU/ FLAGPU
      COMMON /FLAGPULSE/ EXTPOL, CORREC
      integer id,ndep
      real x_coor, y_coor, x_len
      logical out_plane
      real arrayshape(10)
      integer xplane
      common /arrayparm/out_plane,arrayshape,xplane
      real xhor(1000), yhor(1000),zhor(1000),length
      real ang(2,3),tht
      common /hparmx/xhor
      common /hparmy/yhor
      common /hparmz/zhor
      real  ranout(*),zdepout(*)
      real len_inv, xlen_inv
c AARON PARABOLIC
      real xcoorarr(1000),ycoorarr(1000)
      real x, y, z, xp, yp, zp, xpp, ypp, zpp, xppp, yppp
      integer lwascomp,iWriteTrf,ierrinfile
      common /logcomp/lwascomp,iWriteTrf,ierrinfile
     
c      write(*,*)' *************** subroutine geom'
c      write(7,*)'zhor',(zhor(id),id=1,ndep)
      rngref=ranout(1)
      
      if (xplane.le.0) then
         x_len =arrayshape(1) - arrayshape(2)**2*8/3/arrayshape(1)
         length= arrayshape(1)
      elseif (xplane.eq.1) then
         length= xhor(ndep)
         x_len =xhor(ndep)*(1 - 8/3* (arrayshape(2)/xhor(ndep))**2)      
      else
         x_len =xhor(ndep)      
      endif
      xlen_inv=1./x_len
      len_inv=1./length

C  z axis is down, x axis away from source, y axis to the rigth of x axis

c  arrayshape(1) length (positive downwards)
c  arrayshape(2) maximum offset of parabola (positive away from source)
c  arrayshape(3) tilt    (counterclockwise around Y-axis ) 
c  arrayshape(4) azimuth (counterclockwise around z-axis )
c  arrayshape(5)swivel  (counterclockwise around z-axis )
        if (xplane.le.2) then
c-------- compute cosines
           do id=1,3
              tht = arrayshape(id+2)*0.0174532925
              ang(1,id) = cos(tht)
              ang(2,id) = sin(tht)
           enddo
        endif

c      write(7,*)'xplane,zhor(id),zhor(1),x_len,len_inv'
c      write(7,*)xplane,zhor(id),zhor(1),x_len,len_inv
      DO 3300   ID=1,ndep
         ZRDEP=zhor(id)
c
            if (xplane.eq.3) then
c-----   No rotations
               x_coor=xhor(id)
               y_coor=yhor(id)
            elseif (xplane.eq.2) then
c---- User defined
               z= xhor(id)
               y= 0
               x=yhor(id)
            elseif (xplane.le.1) then
c----  Parabola
               if (xplane.lt.1) Z=(zhor(id)-zhor(1))*x_len*len_inv
               if (xplane.eq.1) z= xhor(id)         *x_len*len_inv
               y=0.
               x=4*arrayshape(2)*(z*xlen_inv -(z*xlen_inv)**2)
            endif
c
c---- Rotating
c
            if (xplane.le.2) then
c              write(7,*) 'initial x,y,z',x,y,z
c rotate by swiveling around z axis, around top phone
c               call rotate(xp,yp,x,y,arrayshape(3))
               call rotate(xp,yp,x,y,ang(1,1))
               zp=z
C	       write(7,*) 'xp,yp,zp',xp,yp,zp
c tilt array by rotating around yp axis	   
c               call rotate(zpp,xpp,zp,xp,arrayshape(4))
               call rotate(zpp,xpp,zp,xp,ang(1,2))
               ypp=yp
C	       write(7,*) 'xpp,ypp,zpp',xpp,ypp,zpp
c set azimuth by rotating around zpp axis
c               call rotate(xppp,yppp,xpp,ypp,arrayshape(5))
               call rotate(xppp,yppp,xpp,ypp,ang(1,3))
               zrdep=zhor(1)+zpp
               x_coor=xppp
               y_coor=yppp
            endif
c
c---- Find range
c
            ranout(id)=sqrt( (rngref+x_coor)**2 + y_coor**2)
            zdepout(id)=zrdep
            ycoorarr(id)=y_coor
            xcoorarr(id)=x_coor
c           write(7,*)'geo',id,ranout(id),zdepout(id), x_coor,y_coor
 3300  continue

       if ((iwritetrf.eq.1) .or. (flagpu.lt.5)) then        
          write(7,*)'xplane',xplane
          write(7,*)'length,x-axis length'  ,length,  x_len    
          write(7,*)'length, top_of parabola',
     &         arrayshape(1),arrayshape(2)
          write(7,*)'Angles, tilt, rotation, swiwel',
     &         arrayshape(3), arrayshape(4),arrayshape(5)
          write(7,*)'source_range below is the total source range'
          write(7,*)'rec_no,source_range,  z_coor, x_coord,  y_coord'
          do id=1,ndep
             write(7,*)'geo',id,ranout(id),zdepout(id),
     &                          xcoorarr(id),ycoorarr(id)
          enddo
       endif
       end


