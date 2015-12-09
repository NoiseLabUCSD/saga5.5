      SUBROUTINE ARGGET(N,ARG,LARG)
      INTEGER       N,LARG
      CHARACTER*(*) ARG
      INTEGER       IARGC,LENSTR2
      EXTERNAL      IARGC,GETARG,LENSTR2
      
      IF ((N.GT.IARGC()).OR.(N.LT.1)) THEN
         WRITE(*,*) 'Invalid argument number'
         STOP
      END IF
      CALL GETARG(N,ARG)
      LARG = LENSTR2(ARG)

      RETURN
      END
