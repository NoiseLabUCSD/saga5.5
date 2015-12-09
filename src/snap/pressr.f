      SUBROUTINE PRESSR(NP1,NP2,RNG,PRESS,ALFA,EK,US,UR,SAVE)

      DIMENSION RNG(*)
      DIMENSION ALFA(*), EK(*), US(*), UR(*), SAVE(*)

      DOUBLE PRECISION EK,  A, SAVE
      DOUBLE PRECISION H0, H1, TWOPI, PI, OMEGA, DQ2PI

      COMPLEX ARG
      COMPLEX PRESS(*)

      COMMON /FLAGSnap/ PLANE, FFP, EK0, SQEK0
      COMMON /GSNAP/ H0, H1, TWOPI, PI, OMEGA
      COMMON /N/ MINMOD, MAXMOD, HSTART, NFREQ
      PARAMETER ( ICF=6)
      real      FLAGopt(ICF)
      common /flagarray/flagopt

c      write(*,*)' entering pressr'
     
      DQ2PI=DSQRT(TWOPI)
      IF(PLANE .LT. 1.0)   THEN
        OFFSET=TWOPI/8.0
        DO 1200   J=1,MAXMOD-MINMOD+1
           SAVE(J)=(US(J)*UR(J))/DSQRT(EK(J))
 1200   CONTINUE
      ELSE
        OFFSET=0.0
        DO 1400   J=1,MAXMOD-MINMOD+1
          SAVE(J)=(US(J)*UR(J))/EK(J)
 1400   CONTINUE
      END IF
c      write(6,*)'pressr: save(1),ek(1)',save(1),ek(1)
c       OFFSET=0.0
c       write(6,*)'J,US(J),UR(J),EK(J),maxmod',maxmod

c      write(*,*)' pressr:OFFSET,plane ',OFFSET,plane,twopi,pi
 
CCC
      DO    I=NP1,NP2
c         write(*,*)' i',i
         PRESS(I)=CMPLX(0.0,0.0)
      enddo
C
      if (flagopt(6) .eq. 1) then
c        write(*,*) 'incoherent addition'
        DO I=NP1,NP2
          DO 1600   J=1,MAXMOD-MINMOD+1
            A=SAVE(J)*EXP(-ALFA(J)*RNG(I))
            PRESS(I)=PRESS(I) + A*A
 1600     CONTINUE
          press(i)=sqrt(press(i))
        enddo
      else
        DO    I=NP1,NP2
          DO 1601   J=1,MAXMOD-MINMOD+1
            ARG=CMPLX(0.,SNGL(EK(J)*RNG(I) - OFFSET))
            PRESS(I)=PRESS(I) + CEXP(ARG)*SAVE(J)*EXP(-ALFA(J)*RNG(I))
c            write(7,*)'j,ur(j),ek(j),alfa(j),us(j),rng(i)'
c            write(7,*)j,ur(j),ek(j),alfa(j),us(j),rng(i)
 1601     CONTINUE
        enddo
      endif
c---- correct for 3d attenuation
      IF(PLANE .LT. 1.0)   THEN
        DO    I=NP1,NP2
          PRESS(I)= CMPLX(0.,-1.) * PRESS(I) * DSQRT(TWOPI/RNG(I))
        enddo
      ELSE
        DO    I=NP1,NP2
          PRESS(I)= CMPLX(0.,-SQEK0) * PRESS(I) * DQ2PI
        enddo
      END IF
C
 2000 CONTINUE
      END



