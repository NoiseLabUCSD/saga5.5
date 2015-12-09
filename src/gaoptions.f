      SUBROUTINE gaoptions
c     decodes the option string in optinv
      USE global
      INCLUDE 'comopt.h'
      INTEGER i,j,itemp
c---  default values
      iopt (7) = 1              ! normal stress
      iopt(14) = 1
      iopt (2) = 99
      iopt(13) = 99
      iopt (5) = 99
      isubopt(36) = 0
      DO i=1,mopt
         isubopt(i)=0
         itrans(i)=0
      ENDDO
c     
c---  decode the option string
c     
      j=0
      DO 50 i=1,mopt-2
         j=j+1
         IF (j.GT.mopt) GOTO 50

         IF( optinv(j).EQ.'w') THEN
            iopt(1)=1
            WRITE(prtfil,'('' inversion in wavenumber-depth domain'')')
         ELSE IF (optinv(j).EQ.'r') THEN
            iopt(1)=2
            WRITE(prtfil,'('' inversion in range-depth domain'')') 
         ELSE IF (optinv(j).EQ.'d') THEN
            IF ((iopt(2).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies
     &              two data formats'
            ENDIF
            iopt(2)=1
            WRITE(prtfil,'('' inversion reading real data'')') 
         ELSE IF (optinv(j).EQ.'D') THEN
            IF ((iopt(2).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              data formats'
            ENDIF
            iopt(2)=2
            WRITE(prtfil,'('' inversion reading real data'')') 
         ELSE IF (optinv(j).EQ.'e') THEN
            IF ((iopt(2).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              data formats'
            ENDIF
            iopt(2)=3
            WRITE(prtfil,*)' inversion reading real data in VA-format' 
         ELSE IF (optinv(j).EQ.'T') THEN
            IF ((iopt(2).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              data formats'
            ENDIF
            iopt(2)=4
            WRITE(prtfil,*)' inversion reading real data in HA-format' 
c     
c     ----- transforming DATA
c     
         ELSE IF (optinv(j).EQ.'s') THEN
            iopt(3)=1
            itrans(1)=1
            WRITE(prtfil,'('' inversion multiplied with sqrt(r)'')') 
         ELSE IF (optinv(j).EQ.'R') THEN
            iopt(3)=3
            itrans(2)=3
            WRITE(prtfil,'('' inversion divided by sqrt(r)'')') 
         ELSE IF (optinv(j).EQ.'l') THEN
            iopt(3)=2
            itrans(3)=2
            WRITE(prtfil,'('' inversion on a dB scale'')') 
            IF ((ICHAR(optinv(j+1)).GE.48).and
     &           .(ICHAR(optinv(j+1)).LE.57)) THEN
               
               READ(optinv(j+1),'(i3)')itemp
c     WRITE(*,*)'itemp',itemp,ICHAR(optinv(j+1))
               IF (itemp.EQ.2 ) THEN 
                  WRITE(prtfil,*)
     &                 ' ...both observed and calculated',
     &                 ' data are scaled' 
                  iopt(3)=12
                  itrans(3)=12
                  j=j+1
               ENDIF
            ENDIF
         ELSE IF (optinv(j).EQ.'M') THEN
            IF ((ICHAR(optinv(j+1)).GE.48) .and.
     &           (ICHAR(optinv(j+1)).LE.57)) THEN
               
               READ(optinv(j+1),'(i3)') itemp
c     WRITE(*,*)'itemp',itemp,ICHAR(optinv(j+1))
               IF (itemp.EQ.1 ) THEN 
                  WRITE(prtfil,*)
     &                 ' ...observed data scaled with user',
     &                 ' supplied weight' 
                  iopt(3)=40
                  itrans(4)=40
                  j=j+1
               ELSEIF (itemp.EQ.2 ) THEN
                  WRITE(prtfil,*)
     &                 ' ...calculated data scaled with user',
     &                 ' supplied weight' 
                  iopt(3)=4
                  itrans(4)=4
                  j=j+1
               ELSE
                  WRITE(*,*)'suboption',itemp,'not valid for opt M'
                  STOP
               ENDIF
            ELSE
               iopt(3)=44
               itrans(4)=44
               WRITE(prtfil,*)
     &              ' ...both calculated and observed',
     &              ' data scaled with user supplied weight' 
               WRITE(*,*)
     &              ' ...both calculated and observed',
     &              ' data scaled with user supplied weight' 
            ENDIF
c     WRITE(*,*)iopt(3)-10*INT(iopt(3)/10),iopt(3), INT(iopt(3)/10)
         ELSE IF (optinv(j).EQ.'G') THEN
            iopt(25)=1
            iopt(3)=5
            itrans(5)=5
            WRITE(prtfil,'('' Based on magnitude of observations'')') 
c     
c---- trans transformations
c     
         ELSE IF (optinv(j).EQ.'a') THEN
            iopt(4)=1
            WRITE(prtfil,'('' inversion using simulated annealing'')') 
            IF (ICHAR(optinv(j+1)).EQ.49) THEN ! 49 =1
               WRITE(prtfil,*)
     &              ' ...in each step all parameters are changed' 
               isubopt(4)=1
               j=j+1
            ENDIF
         ELSE IF (optinv(j).EQ.'v') THEN
            iopt(4)=3
            isubopt(4)=1
            WRITE(prtfil,
     $           '('' inversion using very fast'',
     &           '' simulated annealing'')') 
            IF (ICHAR(optinv(j+1)).EQ.49) THEN ! 49 =1
               WRITE(prtfil,*)
     &              ' ...in each step ONE parameter is changed' 
               isubopt(4)=0
               j=j+1
            ENDIF
         ELSE IF (optinv(j).EQ.'S') THEN
            iopt(4)=4
            isubopt(4)=1
            WRITE(prtfil,
     $           '('' MCMC sampling of posteriori probability'')') 
            IF (ICHAR(optinv(j+1)).EQ.49) THEN ! 49 =1
               WRITE(prtfil,*)
     &              ' ...in each step ONE parameter is changed' 
               isubopt(4)=0
            ELSEIF (ICHAR(optinv(j+1)).EQ.50) THEN ! 50 = 1
               WRITE(prtfil,*)
     &              ' ...Enumerative integration' 
               isubopt(4)=2
               j=j+1
            ENDIF
            IF (ICHAR(optinv(j+1)).EQ.49 .OR. 
     &           ICHAR(optinv(j+1)).EQ.48) THEN ! 49 =1, 48=0
               j=j+1
               IF (ICHAR(optinv(j)).EQ.49) THEN
                  WRITE(prtfil,*)
     &                 ' ...Adaptive search interval' 
                  isubopt(35)=1
                  j=j+1
               ENDIF
            ENDIF
c     .cfh.  isubopt(36)=1 is used for optimizing nu
         ELSE IF (optinv(j).EQ.'*') THEN     
            isubopt(36) = 1
            WRITE(prtfil,'('' Optimizing for nu (for orca & snap)'')') 
         ELSE IF (optinv(j).EQ.'g') THEN
            iopt(4)=2
            WRITE(prtfil,'('' inversion using Gauss Newton'')') 
         ELSE IF (optinv(j).EQ.'N') THEN
            iopt(5)=0
            WRITE(prtfil,'('' Cost based on |a-b|^2'')') 
         ELSE IF (optinv(j).EQ.'n') THEN
            IF ((iopt(5).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies
     &              two obj functions'
            ENDIF
            iopt(5)=1
            WRITE(prtfil,'('' cost function (abs(a)-abs(b))'')') 
         ELSE IF (optinv(j).EQ.'X') THEN
            IF ((iopt(5).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              obj functions'
            ENDIF
            iopt(5)=3
            WRITE(prtfil,*)' cost function (|a|-|b|)^2, ',
     &           'no offset correc.' 
         ELSE IF (optinv(j).EQ.'f') THEN
            IF ((iopt(5).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              obj functions'
            ENDIF
            iopt(5)=4
            WRITE(prtfil,'('' Incoherent addition over'',
     &           '' frequencies'')') 
            IF (ICHAR(optinv(j+1)).EQ.49) THEN ! 49 =1
               WRITE(prtfil,*)
     &              ' ...with each frequency weighted identicallly' 
               isubopt(5)=1
               j=j+1
            ENDIF
         ELSE IF (optinv(j).EQ.'F') THEN
            IF ((iopt(5).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              obj functions'
            ENDIF
            iopt(5)=5   
            isubopt(5)=1
            WRITE(prtfil,'(''  Frequency domain matched filter'')') 
            IF (ICHAR(optinv(j+1)).EQ.49) THEN ! 49 =1
               WRITE(prtfil,*)
     &    ' ...with each frequency summed and weigthed by power' 
               isubopt(5)=1
               j=j+1
            ENDIF
            IF (ICHAR(optinv(j+1)).EQ.50) THEN ! 49 =2
               WRITE(prtfil,*)
     &    ' ...with each frequency summed and weighted identicallly' 
               isubopt(5)=2
               j=j+1
            ENDIF
           IF (ICHAR(optinv(j+1)).EQ.51) THEN ! 49 =3
               WRITE(prtfil,*)
     &    ' ...with each frequency multiplied and weighted by power' 
               isubopt(5)=3
               j=j+1
            ENDIF
            IF (ICHAR(optinv(j+1)).EQ.52) THEN ! 49 =4
               WRITE(prtfil,*)
     &    ' ...with each frequency multiplied, weighted identicallly' 
               isubopt(5)=4
               j=j+1
            ENDIF
          ELSE IF (optinv(j).EQ.'j') THEN
            IF ((iopt(5).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              obj functions'
            ENDIF
            iopt(5)=6
            WRITE(prtfil,'('' Phase and magnitude depth processor'')') 
         ELSE IF (optinv(j).EQ.'k') THEN
            IF ((iopt(5).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              obj functions'
            ENDIF
            iopt(5)=7
            WRITE(prtfil,'(''Coherent in range and incoherent'',
     &           '' in freq'')') 
            IF (ICHAR(optinv(j+1)).EQ.49) THEN ! 49 =1
               WRITE(prtfil,*)
     &              ' ...with each frequency weithed identicallly' 
               isubopt(5)=1
               j=j+1
            ENDIF
         ELSE IF (optinv(j).EQ.'Y') THEN
            IF ((iopt(5).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              obj functions'
            ENDIF
            iopt(5)=8
            WRITE(prtfil,'(''Search for lag'')') 
            IF (ICHAR(optinv(j+1)).EQ.49) THEN ! 49 =1
               WRITE(prtfil,*)
     &              ' ...with each frequency weighted identicallly' 
               isubopt(5)=1
               j=j+1
            ENDIF
         ELSE IF (optinv(j).EQ.'J') THEN
            IF ((iopt(5).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              obj functions'
            ENDIF
            iopt(5)=10
            WRITE(prtfil,'(''Range uncertainty incorporated'')') 
c     IF (ICHAR(optinv(j+1)).EQ.49) THEN   ! 49 =1
c     WRITE(prtfil,*)
c     &             ' ...with each frequency weighted identicallly' 
c     isubopt(5)=1
c     j=j+1
c     ENDIF

         ELSEIF( optinv(j).EQ.'Q') THEN
            iopt(6)=1
            WRITE(prtfil,'('' inversion debugging'')') 
         ELSE IF (optinv(j).EQ.'C') THEN
            iopt(8) = 1
            WRITE(prtfil,'(''contour of objective function'')') 
            IF ((ICHAR(optinv(j+1)).EQ.49)) THEN
               iopt(8)=2
               j=j+1
               WRITE(prtfil,'(''   ...No scaling of maximum '')') 
            ELSEIF ((ICHAR(optinv(j+1)).EQ.50)) THEN
               j=j+1
               iopt(8)=3
               WRITE(prtfil,'(''   ...Plotting as log(1-phi) '')') 
            ELSEIF ((ICHAR(optinv(j+1)).EQ.51)) THEN
               iopt(8)=4
               j=j+1
               WRITE(prtfil,'(''   ...Plotting as log(1-phi) '')') 
            ENDIF
         ELSE IF (optinv(j).EQ.'p') THEN
            iopt(10)=1
            WRITE(prtfil,*)' plot of best model and the initial model' 
            itemp=0
            IF ((ICHAR(optinv(j+1)).GE.48).and
     &           .(ICHAR(optinv(j+1)).LE.57)) THEN       
               READ(optinv(j+1),'(i3)')itemp
c     WRITE(*,*)'itemp',itemp,ICHAR(optinv(j+1))
               IF (itemp.EQ.2 ) THEN 
                  WRITE(prtfil,*)
     &                 ' ...both phase and magnitude is plotted' 
               ELSEIF  (itemp.EQ.3) THEN
                  WRITE(prtfil,*)' ...both phase and magnitude',
     &                 ' and bartlett power is plotted' 
               ELSE
                  WRITE(*,*) ' Error in specifing option ''p'''
                  STOP
               ENDIF
               j=j+1
               iopt(10)=itemp
            ENDIF
         ELSE IF (optinv(j).EQ.'z') THEN
            iopt(11)=1
            WRITE(prtfil,*)' Noise is added to the data in',
     &           ' the *.obs file' 
         ELSE IF (optinv(j).EQ.'c') THEN
            IF ((iopt(2).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              data formats'
            ENDIF
            IF ((iopt(5).NE.99).OR.(iopt(13).NE.99)) THEN
               STOP 'The options in input file specifies two
     &              obj functions'
            ENDIF
            iopt(13)=1
            iopt(2)=1
            WRITE(prtfil,'('' reading in covariance matrix '')') 
         ELSE IF (optinv(j).EQ.'I') THEN
            iopt(14)=0
            WRITE(prtfil,'('' plotting ppd condensed on a 0-1'',
     &           '' scale '')') 
         ELSE IF (optinv(j).EQ.'i') THEN
            iopt(14)=1
            WRITE(prtfil,'('' plotting ppd on a full scale '')') 
         ELSE IF (optinv(j).EQ.'H') THEN
            iopt(15)=1
            WRITE(prtfil,'('' Using the hybrid optimization scheme'')')
         ELSE IF (optinv(j).EQ.'m') THEN
            iopt(16)=1
            READ(optinv(j+1),'(i3)')ippd1
            READ(optinv(j+2),'(i3)')ippd2
            j=j+2
            WRITE(prtfil,'(a,i3,'' and'',i3)')
     &           'computing 2D marginal PPD between parameters',
     &           ippd1,ippd2 
         ELSE IF (optinv(j).EQ.'E') THEN
            iopt(17)=1
            WRITE(prtfil,'('' Using EOF.... '')') 
         ELSE IF (optinv(j).EQ.'A') THEN
            iopt(19)=1
            WRITE(prtfil,'('' Adding the input values to PPD-plot '')')
         ELSE IF (optinv(j).EQ.'b') THEN
            iopt(20)=1
            WRITE(prtfil,'('' The covariance is divided by Ndep '')') 
         ELSE IF (optinv(j).EQ.'B') THEN
            iopt(20)=2
            WRITE(prtfil,
     &           '('' The covariance is divided'',
     &           '' by lagest eigenvalue '')') 
         ELSE IF (optinv(j).EQ.'W') THEN
            iopt(21)=1
            WRITE(prtfil,'('' Observed data written '',
     &           ''out to unit 30 '')') 
         ELSE IF (optinv(j).EQ.'U') THEN
            iopt(22)=1
            WRITE(prtfil,'('' uncertainty for local methods '')') 
         ELSE IF (optinv(j).EQ.'q') THEN
            iopt(23)=1
            WRITE(prtfil,'('' A priori uncertainty included '')') 
         ELSE IF (optinv(j).EQ.'x') THEN
            iopt(24)=1
            WRITE(prtfil,'(a,a)')' Starting model INCLUDED,',
     &           ' as an initial',
     &           ' member of the temp0 populations' 
         ELSE IF (optinv(j).EQ.'L') THEN
            iopt(9)=1
            WRITE(prtfil,'('' Regularization....'')') 
         ELSE IF (optinv(j).EQ.'h') THEN
            iopt(27)=1
c     iopt(3)=10
            WRITE(prtfil,*)' Obseverd and calculated data norm.',
     &           ' to unit'
         ELSE IF (optinv(j).EQ.'O') THEN
            iopt(18)=1
            WRITE(prtfil,*)' Bartlett powers are summed ',
     &           'logarithmically over the frequncies'
         ELSE IF (optinv(j).EQ.'t') THEN
            iopt(26)=1
            WRITE(prtfil,*)' Bartlett power is written to *.env file',
     &           '          only for options c,f,F'
         ELSE IF (optinv(j).EQ.'o') THEN
            iopt(18)=2
            WRITE(prtfil,*)' The source power is constant',
     &           ' over the frequncies'
         ELSE IF (optinv(j).EQ.'u') THEN
            iopt(28)=2
            WRITE(prtfil,*)' The probability distributions ',
     &           'are unscaled'
         ELSE IF (optinv(j).EQ.'P') THEN
            iopt(28)=1
            nummodes=0
 111        IF ((ICHAR(optinv(j+1)).GE.48).and
     &           .(ICHAR(optinv(j+1)).LE.57)) THEN
               
               READ(optinv(j+1),'(i3)',err=112) itemp
               WRITE(*,*)'itemp',itemp,ICHAR(optinv(j+1))
               j=j+1
               nummodes=nummodes*10 + itemp
               GOTO 111
            ENDIF
 112        IF (nummodes.EQ.0) THEN   
               nummodes=5
            ENDIF
            WRITE(*,*)' number of modes for objective fct', nummodes
c     iopt(3)=10
            WRITE(prtfil,*)' New post processing'
         ELSE IF (optinv(j).EQ.'Z') THEN
            iopt(29)=1
            WRITE(prtfil,*)' Single point crossover'
         ELSE IF (optinv(j).EQ.'V') THEN
            iopt(31)=1
            WRITE(prtfil,*)' No checking of data from cost'
         ELSE IF (optinv(j).EQ.'K') THEN
            iopt(32)=1
            WRITE(prtfil,*)' Lineplot of sensitivity',
     &           ' for each parameter'
            IF ((ICHAR(optinv(j+1)).EQ.49)) THEN
               iopt(32)=2
               j=j+1
               WRITE(prtfil,'(''   ...No scaling of maximum '')') 
            ELSEIF ((ICHAR(optinv(j+1)).EQ.50)) THEN
               iopt(32)=3
               j=j+1
               WRITE(prtfil,'(''   ...Plotting as log(1-phi) '')') 
            ENDIF
         ELSE IF (optinv(j).EQ.'y') THEN
            iopt(33)=1
            WRITE(prtfil,*)' Gray coding is used in',
     &           ' the GA optimization'
         ELSE IF (optinv(j).EQ.'?') THEN
            iopt(34)=1
            WRITE(prtfil,*)' random seeds for each run'
         ELSE IF (optinv(j).EQ.'+') THEN
            iopt(35) = 1
            WRITE(prtfil,*)' uncertainty covariance matrix is included'
            WRITE(prtfil,'('' reading in unc. cov. matrix '')') 
         ELSE IF (optinv(j).EQ.'!') THEN
            GOTO 60
         ELSE IF (optinv(j).NE.' ') THEN
            WRITE(prtfil,399) optinv(j)
 399        FORMAT(' >>>> unknown inversion OPTION: ',A1,' <<<<')
         END IF
 50   CONTINUE
 60   CONTINUE
      IF ((iopt(5).NE.99)) THEN
         iopt(13)=0
      ELSEIF ((iopt(13).NE.99)) THEN
         iopt(5)=0
      ENDIF
      IF ((iopt(5).EQ.99).AND. (iopt(13).EQ.99)) THEN
         WRITE(*,*)' The objective function must be specified in'
         WRITE(*,*)'  option line. Specify one of the objective '
         WRITE(*,*)'  functions in the manual. '
         STOP
      ENDIF
      IF ((iopt(2).EQ.99).AND. (iopt(13).EQ.0)) THEN
         WRITE(*,*) 'The data format must be specified in option line'
         WRITE(*,*) ' Specify either of options: c, d, D,  e or T'
         STOP
      ENDIF
      IF ((iopt(13).EQ.1) .AND. (itrans(3).GE.1)) THEN
         WRITE(*,*) ' option c (covariance matrix) not compatible with'
         WRITE(*,*) ' the transformation of data'
         STOP     
      ENDIF
c     
c---- check for line and coutour plot for option f
c     
      IF  (iopt(5).EQ.4) THEN
c         IF ((iopt(32).NE.0) .OR.(iopt(8).NE.0)) THEN
c     THEN we should USE option f1
c            IF (isubopt(5).NE.1) THEN
c               WRITE(*,*) ' *************'
c               WRITE(*,*) ' Line or contur plotting only with '
c               WRITE(*,*) ' objective function f1'
c               STOP
c           ENDIF
c            IF ((iopt(8).eq.1).or.(iopt(8).eq.2)) THEN
c               WRITE(*,*) '*************'
c               WRITE(*,*) 'Contour  plotting only with option C2 or C3',
c     1              'for objective function f1'
c               STOP
c            ELSEIF ((iopt(32).eq.1).or.(iopt(32).eq.2)) THEN
c               WRITE(*,*) '*************'
c               WRITE(*,*) 'Line  plotting only with option K2 or K3',
c     1              ' for objective function f1'
c               STOP
c            ENDIF
c         ENDIF
      ENDIF
c     
c---- check for line and coutour plot for option c
c     
      IF  (iopt(13).EQ.1) THEN
         IF ((iopt(32).NE.0) .OR.(iopt(8).NE.0)) THEN
c     THEN we should USE option f1
            IF (iopt(20).NE.1 .and.(iopt(8).eq.3)) THEN
               WRITE(*,*) '*************'
               WRITE(*,*) 'Line or contur plotting C2 only with ',
     1              'objective function cb ',
     2               'or alternatively with contour plotting C3 alone'
               STOP
            ENDIF
            IF ((iopt(8).eq.1).or.(iopt(8).eq.2)) THEN
               WRITE(*,*) '*************'
               WRITE(*,*) 'Contour  plotting only with option C2 or C3',
     1              ' for objective function c'
               STOP
            ELSEIF ((iopt(32).eq.1).or.(iopt(32).eq.2)) THEN
               WRITE(*,*) '*************'
               WRITE(*,*) 'Line  plotting only with option K2 or K3',
     1              ' for objective function c'
               STOP
            ENDIF
         ENDIF
      ENDIF
      END





