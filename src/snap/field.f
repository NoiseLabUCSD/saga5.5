      SUBROUTINE FIELD(NP,F,IFREQ,MSP,XTS,
     & TEMPOR,RNG,FLDPR,FLDPR1,US,UR,ALFA,EK,TLUS,SRD,
     & ND0,MODEN,KSRD,rdep,ndep,nsrd,nrange,sourcesp,*)
c    ksrd is the maximum number of sources
      integer flagpu
      COMMON /FLAGPU/ FLAGPU
       real     CORREC
      COMMON /FLAGPULSE/  EXTPOL, CORREC
      REAL XTS(moden,MSP)
      REAL RNG(*),rdep(*)
      REAL TLUS(ksrd,MODEN)
      REAL US(MODEN), UR(MODEN)
      REAL ALFA(MODEN)
      REAL SRD(KSRD)
      complex sourcesp(*)

      DOUBLE PRECISION H0, H1, TWOPI, PI, OMEGA
      DOUBLE PRECISION EK(MODEN), TEMPOR(MODEN)

      COMPLEX FLDPR(nrange)
      COMPLEX FLDPR1(nrange)
C
      COMMON /CFIELD/ FLDRD(3), PULRD(3)
      COMMON /FLAGS/ PLANE, FFP, EK0, SQEK0
      COMMON /GSNAP/ H0, H1, TWOPI, PI, OMEGA
      COMMON /LUNIT/ LUPLP, LUPLT, LUPRT
      COMMON /N/ MINMOD, MAXMOD, HSTART, NFREQ
      real rngref
c AARON
      logical tilt,offset,timeoffset,msource,mbeam,msourcesp,
     1     two_way,multistatic
      real dtilt,multigeo(4)
      real dr_ind(256),time_ind(256)
      common /tiltparm/tilt,offset,timeoffset,msource,mbeam,
     2    msourcesp,
     1       two_way,multistatic,dtilt,multigeo,dr_ind,time_ind
      logical out_plane
      real arrayshape(10)
      integer xplane
      common /arrayparm/out_plane,arrayshape,xplane
      real rngref2(1000)
c      write(7,*)' entering field'
C
C  DEFINITION OF CONSTANTS.
C
c       write(*,*)'field: sourcesp(1)',sourcesp(1),sourcesp(2)
      MODQTY=MAXMOD-MINMOD+1
      SDLAST=0.0
      IPLANE=PLANE
      ZRMIN=RDep(1)
      NPZ1=1
      NPZ2=ndep
      NPZT=NPZ2-NPZ1+1
      DO ID=NPZ1,NPZ2*np
          FLDPR(id)=(0.,0.)
      enddo
C
C  LOOP FOR STORING MODE AMPLITUDES AT SOURCE DEPTH.
C
      DO 2100   K=1,Nsrd
         IF(SDLAST .NE. SRD(K))   THEN
            SDLAST=SRD(K)
            CALL MODFUN(SDLAST,-1,MODQTY,XTS,MSP,US,MODEN,*4000)
         END IF
         DO   M=1,MODQTY
            TLUS(K,M)=US(M)
         enddo
 2100 CONTINUE

         ISKIP=0
C
      rngref=rng(1)
      if (tilt.or.offset) then
        do i=1,np
           rngref2(i)=rng(i)
        enddo
      endif
      ZSLAST=0.0

      DO 3300   ID=NPZ1+ISKIP,NPZ2
         ZRDEP=rdep(id)
         if (out_plane) then
            rng(1)=rng(id)
            if (flagpu.lt.2) then        
               write(7,*)'from field ',id,rng(1),zrdep
            endif
         endif
c
c---- tilted array
c
         if (tilt .and. (.not. out_plane)) then
            rng(1)=rngref+dtilt/(npzt-1)*(id-1)
            do i=2,np
               rng(i)=rngref2(i)+dtilt/(npzt-1)*(id-1)
            enddo
         endif
c AARON
         if (offset.and.tilt .and. (.not. out_plane)) then
            rng(1)=rng(1)+dr_ind(id)
            do i=2,np
               rng(i)=rng(i)+dr_ind(id)
            enddo
         endif

         if (offset.and.(.not.tilt) .and. (.not. out_plane)) then
            rng(1)=rngref+dr_ind(id)
            do i=2,np
               rng(i)=rngref+dr_ind(id)
            enddo
         endif

        CALL MODFUN(ZRDEP,-1,MODQTY,XTS,MSP,UR,MODEN,*4000)
        DO 3200   IS=1,NSRD
           DO 2500   M=1,MODQTY
              US(M)=TLUS(IS,M)
c     write(*,*)'us m,is', US(M),m,is
 2500      continue   
           IF(PLANE.GT.0.)   THEN
              IF(SRD(IS) .NE. ZSLAST)   THEN
                 ZSLAST=SRD(IS)
                 SQEK0= SQRT(OMEGA/CZS(ZSLAST,Z0,C0,ND0,HSTART))
              END IF
           END IF

      CALL PRESSR(1,NP,RNG,FLDPR1(1+(id-1)*np),ALFA,EK,US,UR,TEMPOR)
C   TRANSFER OF DATA TO OUTPUT FILE.
c      write(*,*)'hello, msourcesp',msourcesp
      if (.not.msourcesp) then
         DO  N1=1,NP
            FLDPR(N1+(id-1)*np)= FLDPR(N1+(id-1)*np)
     1           +conjg(FLDPR1(N1+(id-1)*np))
c       write(*,*)'prc,pr',FLDPR(N1+(id-1)*np),FLDPR1(N1+(id-1)*np)
         enddo
      else
         DO  N1=1,NP
            FLDPR(N1+(id-1)*np)= FLDPR(N1+(id-1)*np)
     1           +conjg(FLDPR1(N1+(id-1)*np))*conjg(sourcesp(is))
         enddo
      endif

 3200 CONTINUE
        if (timeoffset ) then   !
         DO  N1=1,NP
            FLDPR(N1+(id-1)*np)= FLDPR(N1+(id-1)*np)*
     1           exp(cmplx(0.,time_ind(id)*omega))
         enddo

        endif
 
 3300 CONTINUE
c AARON
        if (tilt.or.offset ) then   ! reset the range
          rng(1)=rngref
          do i=2,np
             rng(i)= rngref2(i)
          enddo
        endif
        if (out_plane ) then   ! reset the range
          rng(1)=rngref
        endif

   
      return
 4000 return 1
c      stop 'field: problems with the modes in field'
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
















