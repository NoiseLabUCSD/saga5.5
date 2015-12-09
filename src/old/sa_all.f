      SUBROUTINE SA()
C     
C     THIS CODE USES SIMULATED ANNEALING AND THE PE METHOD TO
C     INVERT OCEAN SEDIMENT PARAMETERS.
C     Adapted from Mike Collins
      USE global
      INCLUDE 'comopt.h'
      
      
      INTEGER i,j,jj,iter,jruns
      REAL temp
      REAL e,de,ebest,arg
      REAL olde,oldf
      INTEGER ibest
C     
      REAL xran
      REAL FBEST(mpar)
      INTEGER iforwtotal


      DO 1000 jruns=1,npop
         iforwtotal=0
c---  generate the random starting point
         DO 4 J=1,Nparm
            F(J)=FMIN(J)+(FMAX(J)-FMIN(J))*RAN2(1)
            FBEST(J)=F(J)
 4       CONTINUE
C     
c     
         CALL setmodelreal(f)
         CALL forw2
         CALL cost(e)
         iforwtotal=iforwtotal+1
c-----For some aplications the sqrt energy should be used.
c     It is not completely logic !
c     
c     e=SQRT(e) 
C     
         IBEST=0
         EBEST=E
         ITER=0
         WRITE(*,*)ITER,ITER,E
         WRITE(*,*)(F(J),J=1,Nparm)
C     
C     THE SIMULATED ANNEALING LOOP.
C     
         TEMP=TEMP0
         WRITE(*,*)' niter',niter

         WRITE(*,*)' temp',temp0
         WRITE(prtfil,*)' temp',temp
         WRITE(prtfil,*)' Mter',mter
         WRITE(prtfil,*)' Niter',Niter
         DO 11 ITER=1,NiTER
C     
C     LINEAR SEARCH (QUENCHING) FOR ITER > MTER.
C     
            IF(ITER.EQ.MTER+1)THEN
               E=EBEST
               DO 8 J=1,Nparm
                  ILIN(J)=1
                  F(J)=FBEST(J)
 8             CONTINUE
            END IF
C     
C     THE PARAMETERS ARE PERTURBED ONE AT A TIME.
C     
            DO 10 J=1,Nparm
               OLDE=E
               OLDF=F(J)
C     -1< xran <1
               XRAN=-1.0+2.0*RAN2(1)
               F(J)=F(J)+DF(J)*XRAN**3
               IF(F(J).LT.FMIN(J))F(J)=2.0*FMIN(J)-F(J)
               IF(F(J).GT.FMAX(J))F(J)=2.0*FMAX(J)-F(J)
C     
               CALL setmodelreal(f)
               CALL forw2
               CALL cost(e)
               iforwtotal=iforwtotal+1
               E=SQRT(E)
               DE=E-OLDE
               IF (de.GT.0) THEN
c---  should an ophill move be accepted ?
                  IF (ilin(j).EQ.0) THEN
                     ARG=AMIN1(60.0,DE/TEMP)
                     XRAN=RAN2(1)*EXP(ARG)
                  ELSE
C     LINEARIZED INVERSION FOR THE J-TH PARAMETER IF ILIN(J)=1.
                     xran=1.
                  ENDIF
C     
C     REVERT TO OLD VALUES IF THE PERTURBATION IS NOT ACCEPTED.
C     
                  IF(XRAN.GT.1.0)THEN
                     E=OLDE
                     F(J)=OLDF
                  ENDIF
               END IF
C     
               IF(E.LT.EBEST)THEN
                  WRITE(60,'(1x,I7,I7,G12.5,100f12.2)')
     &                 jruns,iforwtotal,e,(f(jj),jj=1,nparm)
                  WRITE(*,*)'   **** NEW BEST ENERGY ****'
c     ABEST=1.0/AMP
                  EBEST=E
                  IBEST=ITER
                  DO 9 JJ=1,Nparm
                     FBEST(JJ)=F(JJ)
 9                CONTINUE
                  WRITE(prtfil,*)ITER,E
                  WRITE(prtfil,*)iter,(f(jj),jj=1,nparm)
                  WRITE(*,*)ITER,E
                  WRITE(*,*)(f(jj),jj=1,nparm)
               END IF
C     
 10         CONTINUE            ! Finished metropolis
C     
            TEMP=TEMP0/FLOAT(ITER)
C     
c     WRITE(prtfil,*)ITER,E
c     WRITE(prtfil,*)(F(J),J=1,Nparm)
C     
 11      CONTINUE               ! Finnished with the iterations.
C     
C     BEST PARAMETER VALUES.
C     
c     OPEN(UNIT=11,STATUS='UNKNOWN',FILE='best.dat')
         WRITE(*,*)' Best energy & parameters'
         WRITE(*,*)EBEST,IBEST
c     WRITE(11,*)EBEST,IBEST
         DO 12 J=1,Nparm
            WRITE(*,*)J,FBEST(J)
            WRITE(prtfil,*)J,FBEST(J)
 12      CONTINUE
 200     CONTINUE
 1000 CONTINUE
      WRITE(*,*)'number of forward modelling calls',iforwtotal
      WRITE(*,*)' End of sa'


      END





