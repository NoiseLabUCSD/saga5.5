      SUBROUTINE COHRNT(DEPTH,NP1,NP2,RNG,TLOSS,
     %          US,UROLD,URNEW,MODEN,DELTA,MSP,
     %         ALFOLD,ALFNEW,ALFAVR,EKOLD,EKNEW,EKAVR)
      REAL*8 SUM,EKOLD,EKNEW,EKAVR,EIG,A,SUM1,SUM2,ARG

      DIMENSION ALFOLD(MODEN),ALFNEW(MODEN),ALFAVR(MODEN)
      DIMENSION EKOLD(MODEN),EKNEW(MODEN),EKAVR(MODEN)
      DIMENSION URNEW(2,MODEN),US(MODEN),UROLD(2,MODEN)
      DIMENSION RNG(*)
      complex   TLOSS(*)
      real*8 sum_incoh
      COMMON/NRNGD/HNEW,HOLD,RSTART,REND,SLNGTH,IND1
      COMMON/N/MINMOD,MAXMOD,HORIG,NFREQ
      real*8  twopi
      common /twopie/twopi
      logical tilt,incoh
      real dtilt
      common /tiltparm/tilt,dtilt,incoh
C
      MODQTY=MAXMOD-MINMOD+1
      FACT=0.0
      DO 2000   I=NP1,NP2
         SUM1=0.0D0
         SUM2=0.0D0
         sum_incoh=0.0D0
         RANGE=RNG(I)
         DELTA2=RNG(I)-RSTART
         IF(DEPTH .GT. 0.0)   THEN
            H1=HOLD+((HNEW-HOLD)*DELTA2)/DELTA
            DZ=H1/FLOAT(MSP-1)
            D1=(IND1-1)*DZ
            FACT=(DEPTH-D1)/DZ
         END IF
         DO 1000   J=1,MODQTY
            UR1=UROLD(1,J)+((URNEW(1,J)-UROLD(1,J))*DELTA2)/DELTA
            UR2=UROLD(2,J)+((URNEW(2,J)-UROLD(2,J))*DELTA2)/DELTA
            UR=UR1+(UR2-UR1)*FACT
            ALFA=ALFOLD(J)+((ALFNEW(J)-ALFOLD(J))*DELTA2)/DELTA
            ALFAR=(ALFAVR(J)+0.5*(ALFOLD(J)+ALFA)*DELTA2)
            EIG=EKOLD(J)  +((EKNEW(J)- EKOLD(J) )*DELTA2)/DELTA
            EIG=  (EKAVR(J) +0.5*(EKOLD(J) +EIG )*DELTA2)/RANGE
            A=(US(J)*UR)/DSQRT(EIG)*EXP(-ALFAR)
            ARG=EIG*RNG(I)-twopi/8
            SUM1=SUM1+A*DCOS(ARG)
            SUM2=SUM2+A*DSIN(ARG)
            sum_incoh=sum_incoh+A*A
c      write(30,*)'UR',J,UR1,UROLD(1,J),URNEW(1,J),DELTA2,DELTA
c      write(30,*)'j,ur,eig,alfar/rng(i),us(j),rng(i)'
c      write(30,*)j,ur,eig,alfar/rng(I),us(j),rng(i),sum1,sum2
c      if (j.eq.1) then
c      write(30,*)'alfold(j),alfnew(j),alfavr(j),delta2,delta'
c      write(30,*) alfold(j),alfnew(j),alfavr(j),delta2,delta
c      endif
 1000    CONTINUE
         if (incoh) then
c     write(*,*)'*************incoherent averaging'
            TLOSS(I)=sqrt(sum_incoh)* sqrt(twopi/rng(i))
         else
            TLOSS(I)=cmplx(sum2,sum1)* sqrt(twopi/rng(i))
         endif
c      write(30,*)'Tloss',tloss(I),I
c      write(*,*)'Tloss bef',tloss(I),I
2000  CONTINUE
C  
c      write(*,*)'Tloss(1)',tloss(1)
      RETURN
      END
