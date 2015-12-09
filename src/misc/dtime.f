      FUNCTION DTIME(TD)
      REAL DTIME,TD(2)
      REAL ETIME,TE(2)
      REAL TLAST(2)
      LOGICAL FIRST
      SAVE FIRST,TLAST
      EXTERNAL ETIME
      DATA FIRST /.TRUE./
C      
      IF(FIRST) THEN
         DTIME = ETIME(TLAST)
         TD(1) = TLAST(1)
         TD(2) = TLAST(2)
         FIRST = .FALSE.
      ELSE
         DTIME = ETIME(TE)
         TD(1) = TE(1) - TLAST(1)
         TD(2) = TE(2) - TLAST(2)
         DTIME = TD(1) + TD(2)
         TLAST(1) = TE(1)
         TLAST(2) = TE(2)
      ENDIF
C      
      RETURN
      END
