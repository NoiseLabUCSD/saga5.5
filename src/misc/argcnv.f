      SUBROUTINE ARGCNV(N,FNAME,FLEN,SUFFIX)

C     NAME    : argcnv
C     UPDATE  : November 15, 94
C     PURPOSE : Get the argument with given suffix

      INTEGER       N
      CHARACTER*(*) FNAME,SUFFIX
      INTEGER       IARGC,NMAX,FLEN,SLEN,LENSTR2
      LOGICAL       CND1,CND2,CND3
      EXTERNAL      IARGC,GETARG,LENSTR2

      NMAX = IARGC()
      CND1 = N.GT.NMAX
      CND2 = N.LT.1
      IF(CND1.OR.CND2) THEN
         WRITE(*,*) 'invalid N in argcnv'
         STOP
      END IF

      CALL GETARG(N,FNAME)
      FLEN = LENSTR2(FNAME)
      SLEN = LENSTR2(SUFFIX)

      CND1 = FLEN.GT.SLEN
      CND2 = FNAME(FLEN-SLEN+1:FLEN).EQ.SUFFIX(1:SLEN)
      CND3 = FNAME(FLEN-SLEN:FLEN-SLEN).EQ.'.'
      IF(CND1.AND.CND2.AND.CND3) RETURN
       
      FNAME(FLEN+1:FLEN+SLEN+1) = '.'//SUFFIX(1:SLEN)
      FLEN = LENSTR2(FNAME)
      
      RETURN
      END
