      SUBROUTINE MODFUN(DD,INDEX,MODQTY,XTs,MSP,U,MODEN,*)
      COMMON/N/MINMOD,MAXMOD,HSTART,NFREQ
      DIMENSION XTs(moden,MSP),U(MODEN)

  100 FORMAT(1X,' ***  ERROR DETECTED IN SUB MODFUN :  ',/,
     & '      THIS SUB MAY NOT BE CALLED FOR A DEPTH WHICH EXCEEDS THE',
     & '      WATER DEPTH.',/,
     & '      WATER DEPTH    : ',F10.3,' m',/,
     & '      MODE AMP DEPTH : ',F10.3,' m')          
      DEPTH=DD

c      write(*,*) ' ent sub modfun, depth, hstart : ', depth, hstart
      IF(DEPTH .GT. HSTART)   THEN
        IF( (DEPTH - HSTART) .GT. 1.0E-3)   THEN
          WRITE(6,100) HSTART, DEPTH
          RETURN 1
        ELSE 
          IND= MSP
        END IF
      END IF

      IND= INDEX
      IF(IND .LT. 0)   THEN
        DZ= HSTART/FLOAT(MSP-1)
        IND= DEPTH/DZ+1
        FACT= (DEPTH-(IND-1)*DZ)
c        write(*,*)'modfun,ind,msp',ind,msp
        DO 2000   MODE= 1, MODQTY
          if (ind.ne.msp) then
        U(MODE)= XTS(MODE,IND)+ ((XTS(MODE,IND+1)-XTS(MODE,IND))/DZ)*
     &           FACT
          else
        U(MODE)= XTS(MODE,IND)
          endif
 2000   CONTINUE
      ELSE
        DO 1000   MODE= 1, MODQTY
        U(MODE)= XTS(MODE,IND)
 1000   CONTINUE
      END IF

      RETURN
      END
