      SUBROUTINE STAMP(COMMENT)

C     LOCAL VARIABLES
      CHARACTER*(*) COMMENT
      INTEGER II,IT(3)
      REAL ET(2)
      

C     SYSTEM CALL
      EXTERNAL ETIME,ITIME
      
C     COUNTERS
      DATA II /0/
      SAVE II

      CALL ETIME(ET)
      CALL ITIME(IT)
      IF (II.EQ.0) THEN
         OPEN(UNIT=12,FILE='gems.stp',STATUS='UNKNOWN')
         II=1
      ELSE
         OPEN(UNIT=12,FILE='gems.stp',STATUS='UNKNOWN',
     &        ACCESS='APPEND')
      END IF

      WRITE(12,*) COMMENT
      WRITE(12,1000)
     & 'CPU time (total)       :',ET(1),'second'
      WRITE(12,1100) 
     & 'WALL TIME (HR/MIN/SEC) :',IT(1),'/',IT(2),'/',IT(3)
      WRITE(12,*)
      CLOSE(12)
 
 1000 FORMAT(1X,A24,F15.3,2X,A6)
 1100 FORMAT(1X,A24,5X,I4,A3,I4,A3,I4)

C***** END OF STAMP
      RETURN
      END
