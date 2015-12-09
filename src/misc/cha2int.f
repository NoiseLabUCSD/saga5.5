      SUBROUTINE CHA2INT(CH,N)
      CHARACTER*(*) CH
      INTEGER       N,CLEN,LENSTR2
      EXTERNAL      LENSTR2

      CLEN = LENSTR2(CH)

      N = 0
      DO 10 I=1,CLEN
         N  = N*10
         NN = ICHAR(CH(I:I)) - 48
         N  = N + NN
 10   CONTINUE
      
      RETURN
      END

