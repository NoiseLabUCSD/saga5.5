      SUBROUTINE RANGERD(*,ISECT,NSECT,NPNEW,
     &     ICF,RNG,NRANGE,SECD)
      DIMENSION RNG(NRANGE)
      DIMENSION SECD(3)
      COMMON/N/MINMOD,MAXMOD,HORIG,NFREQ
      COMMON/NRNGD/HNEW,HOLD,RSTART,REND,SLNGTH,IND1
      COMMON/LUNIT/MSOURC,MODOLD,MODNEW,LUPRT
C
      npnew=0
      DELR=SECD(3)
C     DEFINITION OF LAST RANGE VALUE (RMAX) TO BE CALCULATED
      IF((ISECT.EQ.NSECT).OR.(REND.GE.SECD(2))) then
         RMAX=SECD(2)
      else
         IF(SECD(1) .GT. REND)   RETURN 1
         RMAX=REND
      endif
c
c----- FOR DELR=0
c
      if (delr.eq.0) then
c         write(*,*)'range:', secd(1), rstart
         if (secd(1).le.rstart) then
            npnew=0
            return 1
         else
            npnew=1
            RNG(1)=SECD(1)
            GOTO 5000
         endif
      endif
C
C     DEFINITION OF FIRST RANGE VALUE (RMIN) TO BE CALCULATED
      IF((RSTART.GT.SECD(1)))   THEN
       NDELR=((RSTART-SECD(1))/DELR+1.0E-4)
       IF((RSTART-(SECD(1)+NDELR*DELR)).LE.1.0E-4*DELR)
     & NDELR=NDELR-1
       RMIN=SECD(1)+(NDELR+1)*DELR
      ELSE
       RMIN=SECD(1)
      END IF
C
C     DEFINITION OF NUMBER OF POINTS AND FILLING OF RNG ARRAY
C
C     NO POINTS FALLING IN THIS SECTOR
      IF(RMIN.GT.RMAX)  RETURN 1
C

      NPNEW=(RMAX-RMIN)/DELR+1
      IF(RMIN+(NPNEW-1)*DELR.GE.RMAX-1.0E-4*DELR)  NPNEW=NPNEW-1
      IF(NPNEW.gt.NRANGE) then
         NSIZE=NRANGE-5
         WRITE(LUPRT,300) ISECT,NSIZE
 300     FORMAT(1H,/,'  TOO MANY RANGE STEPS IN SECTOR # ',I3,'.',
     $        /,'  MAXIMUM ALLOWED NUMBER AT PRESENT IS ',I4,' .',/,
     $        '  EXECUTION TERMINATED BECAUSE OF ARRAY SIZE',
     $        ' LIMITATIONS.')
         stop
      endif
      DO 4000 IR=1,NPNEW
         RNG(IR)=(IR-1)*DELR+RMIN
 4000 CONTINUE

 5000 CONTINUE
C  ADJUSTMENT FOR THE CASE WHERE RMIN IS TOO SMALL
      IF(RNG(1).LT.1.0E-5)    RNG(1)=1.0E-5
      
      END






