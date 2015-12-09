      SUBROUTINE VELOCITY(MH0I,MSEDI,CREF,SPEED,DH0I,DSEDI,
     &                    C0,Z0,C1,Z1,NPOINT,NDEP)

C___________________________________________________________
C                                                          |
C     This routine calculates the speed of sound.          |
C__________________________________________________________|

      DOUBLE PRECISION SPEED(NPOINT), DH0I, DSEDI, DIV, Z
      DOUBLE PRECISION C0(NDEP), Z0(NDEP), C1(NDEP), Z1(NDEP)
      DOUBLE PRECISION CREF, CNORM, CMIN

      COMMON /NA/ ND0, ND1, CMIN
C
C     WATER
C
      MH0IP1=MH0I+1

      Z=0.
      N=1
      CNORM=1.0D0/CREF
      I=2
      DIV=(C0(2)-C0(1))/(Z0(2)-Z0(1))
      SPEED(1)=C0(1)*CNORM
      DO 2000 N=2,MH0I
      Z=DH0I*(N-1)
      IF(Z.GT.Z0(I)) THEN
 1000   CONTINUE
        I=I+1
        IF(Z.GT.Z0(I)) GO TO 1000
        DIV=(C0(I)-C0(I-1))/(Z0(I)-Z0(I-1))
      END IF
      SPEED(N)=(DIV*(Z-Z0(I-1))+C0(I-1))*CNORM
 2000 CONTINUE
      SPEED(MH0IP1)=C0(ND0)*CNORM
C
C     SEDIMENT
C
      IF(MSEDI.GT.0)   THEN
      I=2
      DIV=(C1(2)-C1(1))/(Z1(2)-Z1(1))
      SPEED(MH0IP1+1)=C1(1)*CNORM
      DO 4000 N=2,MSEDI
      Z=DSEDI*(N-1)+1.D0
      IF(Z.GT.Z1(I)) THEN 
 3000   CONTINUE
        I=I+1
        IF(Z.GT.Z1(I))   GO TO 3000
        DIV=(C1(I)-C1(I-1))/(Z1(I)-Z1(I-1))
      END IF
      SPEED(MH0IP1+N)=(DIV*(Z-Z1(I-1))+C1(I-1))*CNORM
 4000 CONTINUE
      SPEED(MH0IP1+MSEDI+1)=C1(ND1)*CNORM
      END IF
      RETURN
      END
