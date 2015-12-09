      SUBROUTINE RFFORM(STR,LEN,A,ISTP,NUM,IERR)
C
C     READS NUM REAL NUMBERS FROM THE CHARACTER
C     STRING STR IN FREE FORMAT
C
      DIMENSION A(1)
      CHARACTER*(*) STR
      CHARACTER*80 BUF
      LOGICAL FCIF
      IERR=0
      IP=0
      FCIF=(STR(1:1).NE.' ')
      DO 10 I=1,NUM
      A((I-1)*ISTP+1)=0.
  11  IP=IP+1
      IF (IP.GE.LEN) RETURN
      IF (((STR(IP:IP).EQ.' '.OR.STR(IP:IP).EQ.',')
     1    .AND.(STR(IP+1:IP+1).NE.' '.AND.STR(IP+1:IP+1).NE.',')
     2   ).OR.FCIF) THEN
         IF (FCIF) THEN
           IP1=1
           FCIF=.FALSE.
         ELSE
         IP1=IP+1
         END IF
  12     IP=IP+1
         IF (IP.GT.LEN) RETURN
         IF ((STR(IP:IP).EQ.' '.OR.STR(IP:IP).EQ.','.OR.IP.EQ.LEN)
     1      .AND.(STR(IP-1:IP-1).NE.','.AND.
     2            STR(IP-1:IP-1).NE.' ')) THEN
           IP2=IP-1
           LBUF=IP2-IP1+1
           BUF(1:LBUF)=STR(IP1:IP2)
           CALL CON_SR(BUF,LBUF,A((I-1)*ISTP+1),IERR)
         ELSE
           GO TO 12
         END IF
         IP=IP-1
       ELSE
         GO TO 11
       END IF
 10    CONTINUE
      RETURN
      END
      SUBROUTINE CON_SR(BUF,LBUF,A,IERR)
      CHARACTER*(*) BUF
      LOGICAL DEC,SIG,CIF
      A=0
      IL=ICHAR('0')
      IH=ICHAR('9')
      DEC=.FALSE.
      SIG=.FALSE.
      SIVAL=1.
      CIF=.FALSE.
      IDEC=LBUF
      DO 10 I=1,LBUF
      IA=ICHAR(BUF(I:I))
      IF (BUF(I:I).EQ.'-') THEN
        IF (.NOT.SIG.AND.(.NOT.CIF)) THEN
          SIG=.TRUE.
          SIVAL=-1.
        ELSE
          IERR=1
          RETURN
        END IF
      ELSE IF (BUF(I:I).EQ.'+') THEN
        IF (.NOT.SIG.AND.(.NOT.CIF)) THEN
          SIG=.TRUE.
        ELSE
          IERR=1
          RETURN
        END IF
      ELSE
      IF (IA.GE.IL.AND.IA.LE.IH) THEN
        A=A*10.+(IA-IL)
      CIF=.TRUE.
      ELSE IF (BUF(I:I).EQ.'.'.AND.(.NOT.DEC)) THEN
        IDEC=I
      DEC=.TRUE.
      CIF=.TRUE.
      ELSE IF (BUF(I:I).EQ.'E') THEN
        CALL CON_SI(BUF(I+1:LBUF),LBUF-I,IE,IERR)
        IDEC=MIN0(IDEC,I-1)
        A=SIVAL*A*10.0**(IE+IDEC-I+1)
        RETURN
      ELSE
        IERR=1
        RETURN
      END IF
      END IF
 10   CONTINUE
      A=SIVAL*A*10.0**(-LBUF+IDEC)
      RETURN
      END
      SUBROUTINE CON_SI(BUF,LBUF,IA,IERR)
      CHARACTER*(*) BUF
      LOGICAL SIG,CIF
      IA=0
      IL=ICHAR('0')
      IH=ICHAR('9')
      SIG=.FALSE.
      ISIVAL=1
      CIF=.FALSE.
      DO 10 I=1,LBUF
      IC=ICHAR(BUF(I:I))
      IF (BUF(I:I).EQ.'-') THEN
        IF (.NOT.SIG.AND.(.NOT.CIF)) THEN
          SIG=.TRUE.
          ISIVAL=-1
        ELSE
          IERR=1
          RETURN
        END IF
      ELSE IF (BUF(I:I).EQ.'+') THEN
        IF (.NOT.SIG.AND.(.NOT.CIF)) THEN
          SIG=.TRUE.
        ELSE
          IERR=1
          RETURN
        END IF
      ELSE
      IF (IC.GE.IL.AND.IC.LE.IH) THEN
        IA=IA*10.+(IC-IL)
      CIF=.TRUE.
      ELSE
        IERR=1
        RETURN
      END IF
      END IF
 10   CONTINUE
      IA=ISIVAL*IA
      RETURN
      END
