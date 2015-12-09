      SUBROUTINE GETCMDS(NOA,ARGS,LARG)
      LOGICAL       CMD
      CHARACTER*(*) ARGS(*)
      CHARACTER*80  DUMMY
      INTEGER       LARG(*),NOA
      INTEGER       IARGC,LENSTR2,NMAX,DLEN

C     EXTERNAL ROUTINE
      EXTERNAL      IARGC,LENSTR2,GETARG

      NMAX = IARGC()

C     READ ARGUMENTS

      II = 0
      DO 10 I=1,NMAX
         CALL GETARG(I,DUMMY)
         DLEN = LENSTR2(DUMMY)
         CMD = DUMMY(1:1).EQ.'-'
         IF(CMD) THEN
            II = II + 1
            ARGS(II) = DUMMY(2:DLEN)
            LARG(II) = DLEN-1
         END IF
 10   CONTINUE
      NOA = II

      RETURN
      END



