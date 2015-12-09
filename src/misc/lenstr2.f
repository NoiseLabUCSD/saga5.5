      FUNCTION LENSTR2(C)
      INTEGER       LENSTR2
      CHARACTER*(*) C

      LENSTR2 = LEN(C)
      IF(LENSTR2.EQ.0) RETURN

      DO 10 I=1,LENSTR2
         IF(ICHAR(C(I:I)).EQ.32) THEN
            LENSTR2 = I - 1
            RETURN
         END IF
 10         CONTINUE

      RETURN
      END
