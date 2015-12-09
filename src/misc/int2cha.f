      SUBROUTINE INT2CHA(MM, PP, FORM, STRING, NC)

      INTEGER MM, PP, FORM
      CHARACTER*(*) STRING
      INTEGER NC
      CHARACTER*1 BSLASH
      CHARACTER*2 TIMES, UP, DOWN
      CHARACTER*20 WORK, WEXP, TEMP
      INTEGER M, P, ND, I, J, K, NBP
      LOGICAL MINUS
C
C Define backslash (escape) character and escape sequences.
C
      BSLASH = CHAR(92)
      TIMES  = BSLASH//'x'
      UP     = BSLASH//'u'
      DOWN   = BSLASH//'d'
C
C Zero is always printed as "0".
C
      IF (MM.EQ.0) THEN
          STRING = '0'
          NC = 1
          RETURN
      END IF
C
C If negative, make a note of that fact.
C
      MINUS = MM.LT.0
      M = ABS(MM)
      P = PP
C
C Convert M to a left-justified digit string in WORK. As M is a
C positive integer, it cannot use more than 10 digits (2147483647).
C
      J = 10
   10 IF (M.NE.0) THEN
          K = MOD(M,10)
          M = M/10
          WORK(J:J) = CHAR(ICHAR('0')+K)
          J = J-1
       GOTO 10
      END IF
      TEMP = WORK(J+1:)
      WORK = TEMP
      ND = 10-J
C
C Remove right-hand zeros, and increment P for each one removed.
C ND is the final number of significant digits in WORK, and P the
C power of 10 to be applied. Number of digits before decimal point
C is NBP.
C
   20 IF (WORK(ND:ND).EQ.'0') THEN
          ND = ND-1
          P = P+1
       GOTO 20
      END IF
      NBP = ND+MIN(P,0)
C
C Integral numbers of 4 or less digits are formatted as such.
C
      IF ((P.GE.0) .AND. (P+ND.LE.4)) THEN
          DO 30 I=1,P
              ND = ND+1
              WORK(ND:ND) = '0'
   30     CONTINUE
          P = 0
C
C If NBP is 4 or less, simply insert a decimal point in the right place.
C
      ELSE IF (NBP.GE.1.AND.NBP.LE.4.AND.NBP.LT.ND) THEN
          TEMP = WORK(NBP+1:ND)
          WORK(NBP+2:ND+1) = TEMP
          WORK(NBP+1:NBP+1) = '.'
          ND = ND+1
          P = 0
C
C Otherwise insert a decimal point after the first digit, and adjust P.
C
      ELSE
          P = P + ND - 1
          IF (P.EQ.-1) THEN
              TEMP = WORK
              WORK = '0'//TEMP
              ND = ND+1
              P = 0
          ELSE IF (P.EQ.-2) THEN
              TEMP = WORK
              WORK = '00'//TEMP
              ND = ND+2
              P = 0
          END IF
          IF (ND.GT.1) THEN
              TEMP = WORK(2:ND)
              WORK(3:ND+1) = TEMP
              WORK(2:2) = '.'
              ND = ND + 1
          END IF
      END IF
C
C Add exponent if necessary.
C
      IF (P.NE.0) THEN
          WORK(ND+1:ND+6) = TIMES//'10'//UP
          ND = ND+6
          IF (P.LT.0) THEN
              P = -P
              ND = ND+1
              WORK(ND:ND) = '-'
          END IF
          J = 10
   40     IF (P.NE.0) THEN
              K = MOD(P,10)
              P = P/10
              WEXP(J:J) = CHAR(ICHAR('0')+K)
              J = J-1
           GOTO 40
          END IF
          WORK(ND+1:) = WEXP(J+1:10)
          ND = ND+10-J
          IF (WORK(1:3).EQ.'1'//TIMES) THEN
              TEMP = WORK(4:)
              WORK = TEMP
              ND = ND-3
          END IF
          WORK(ND+1:ND+2) = DOWN
          ND = ND+2
      END IF
C
C Add minus sign if necessary and move result to output.
C
      IF (MINUS) THEN
         TEMP = WORK(1:ND)
         STRING = '-'//TEMP
         NC = ND+1
      ELSE
         STRING = WORK(1:ND)
         NC = ND
      END IF
C
C Check result fits.
C
      IF (NC.GT.LEN(STRING)) THEN
          STRING = '*'
          NC = 1
      END IF
      END
