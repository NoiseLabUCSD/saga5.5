      SUBROUTINE SA()
C     
C     THIS CODE USES SIMULATED ANNEALING
C     Adapted from Mike Collins, 1992
c     Peter Gerstoft
      USE global
      INCLUDE 'comopt.h'
      
      
      INTEGER i,j,jj,iter,jruns
      REAL temp,t1
      REAL e,de,ebest,arg
      REAL olde
      INTEGER ibest
C     
      REAL xran,oldf(mpar)
      REAL FBEST(mpar),ftemp(mpar)
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
         WRITE(98,*)e,f(1),f(2),f(1),f(2),fbest(1),fbest(2)
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
               DO jj=1,nparm
                  OLDF(jj)=F(jJ)
               ENDDO
c     WRITE(*,*)' isubopt', isubopt(4)
               IF (isubopt(4).EQ.0) THEN
C     -1< xran <1
                  XRAN=-1.0+2.0*RAN2(1)
                  F(J)=F(J)+DF(J)*XRAN**3
                  IF(F(J).LT.FMIN(J)) F(J)=2.0*FMIN(J)-F(J)
                  IF(F(J).GT.FMAX(J)) F(J)=2.0*FMAX(J)-F(J)
               ELSE
C     -1< xran <1
                  DO jj=1,Nparm
                     XRAN=-1.0+2.0*RAN2(1)
                     F(JJ)=F(JJ)+DF(JJ)*XRAN**3
                     IF(F(JJ).LT.FMIN(JJ)) F(JJ)=2.0*FMIN(JJ)-F(JJ)
                     IF(F(JJ).GT.FMAX(JJ)) F(JJ)=2.0*FMAX(JJ)-F(JJ)
                  ENDDO
               ENDIF
C     
               
               CALL setmodelreal(f)
               CALL forw2
               CALL cost(e)
               iforwtotal=iforwtotal+1
c     E=SQRT(E)
               DE=E-OLDE
               xran=0
               DO i=1,nparm
                  ftemp(i)=f(i)
               ENDDO
               IF (de.GT.0) THEN
c---  should an ophill move be accepted ?
                  IF (ilin(j).EQ.0) THEN
                     ARG=AMIN1(60.0,DE/TEMP)
                     XRAN=RAN2(1)*EXP(ARG)
                  ELSE
C     LINEARIZED INVERSION FOR THE J-TH PARAMETER IF ILIN(J)=1.
                     xran=1.
                     WRITE(*,*)' linear parameter,using quencing ',
     &                    ilin(j)
                  ENDIF
C     
C     REVERT TO OLD VALUES IF THE PERTURBATION IS NOT ACCEPTED.
C     
                  IF(XRAN.GE.1.0)THEN
                     E=OLDE
                     DO jj=1,nparm
                        F(jj)=OLDF(jj)
                     ENDDO
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
                  WRITE(prtfil,*)'Iteration, energy',ITER,E
                  WRITE(prtfil,*)iter,(f(jj),jj=1,nparm)
                  WRITE(*,*)ITER,E
                  WRITE(*,*)(f(jj),jj=1,nparm)
               END IF
C     
               WRITE(98,*)e,ftemp(1),ftemp(2),f(1),f(2),fbest(1),
     &              fbest(2)
 10         CONTINUE            ! Finished metropolis
C     
            IF ((jruns.EQ.1).AND.(iter.EQ.1)) THEN
               CALL rdtime(t1)
               WRITE(prtfil,320) t1*(niter*npop*nparm)/nparm
               WRITE(*,320) t1*(niter*npop*nparm)/nparm
 320           FORMAT(/1H ,' >>> Estimated serial Inversion time: ',
     1              F12.1,' secs.')
            ENDIF
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
