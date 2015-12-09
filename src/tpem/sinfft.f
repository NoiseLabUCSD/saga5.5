       SUBROUTINE SINFFT ( N, X )
C
C***********************************************************************
C*                                                                     *
C*                                                                     *
C* PURPOSE:   SINFFT replaces the real array X()                       *
C*            by its finite discrete sine transform                    *
C*                                                                     *
C* METHOD :                                                            *
C*                                                                     *
C*     The algorithm is based on a mixed radix (8-4-2) real vector     *
C*     fast Fourier synthesis routine published by Bergland:           *
C*                                                                     *
C*     ( G.D. Bergland, 'A Radix-eight Fast Fourier Transform          *
C*      Subroutine for Real-valued Series,' IEEE Transactions on       *
C*      Audio and Electro-acoustics', vol. AU-17, pp. 138-144, 1969 )  *
C*                                                                     *
C*     and sine and cosine transform algorithms for real series        *
C*     published by Cooley, Lewis, and Welch:                          *
C*                                                                     *
C*     (J.W. COOLEY, P.A.W. LEWIS AND P.D. WELSH, 'The Fast Fourier    *
C*     Transform Algorithm: Programming Considerations in the          *
C*     Calculation of Sine, Cosine and Laplace Transforms',            *
C*     J. SOUND VIB., vol. 12, pp. 315-337, 1970 ).                    *
C*                                                                     *
C*                                                                     *
C* ARGUMENTS:                                                          *
C*                       -- INPUT --                                   *
C*                                                                     *
C*      N...... transform size ( = 2**N )                              *
C*                                                                     *
C*      X().... data array dimensioned 2**N in calling program         *
C*                                                                     *
C*                      -- OUTPUT --                                   *
C*                                                                     *
C*      X().... sine transform                                         *
C*                                                                     *
C*  TABLES:  array     required size                                   *
C*                                                                     *
C*              B         2**N                                         *
C*              JINDX     2**(N-1)                                     *
C*              COSTBL    2**(N-4)                                     *
C*              SINTBL    2**(N-4)                                     *
C*                                                                     *
C*  Sub-programs called:  -                                            *
C*                                                                     *
C*                R8SYN..... (radix 8 synthesis)                       *
C*                                                                     *
C***********************************************************************
C
      include 'fftsiz.inc'
      
       INTEGER*4   N
C
       DIMENSION   X(0:*)
       INTEGER*4   NMAX2, NMAX16, NP, NPD2, NPD4
C
       PARAMETER   ( NMAX2  = MAXPTS/2 )
       PARAMETER   ( NMAX16 = MAXPTS/16 )
       DIMENSION   B(MAXPTS), JINDX (NMAX2)
       DIMENSION   SINES (MAXPTS)
       DIMENSION   COSTBL (NMAX16), SINTBL (NMAX16)
C
       SAVE B, COSTBL, JINDX, SINES, SINTBL
       SAVE NSAVE, N4, N8, NP, NPD2, NPD4, NPD16, NPM1
C
C
       DOUBLE PRECISION ARG, DT, PI
       DATA NSAVE  / 0 /
       DATA PI     / 3.1415 92653 58979 32D0 /
C
C----------------------------------------------------------------------+
C
                  
       IF ( N .NE. NSAVE ) THEN 
C                                 compute constants and construct tables
          NSAVE = N
          N8 = NSAVE / 3
          N4 = NSAVE - 3 * N8 - 1
          NP = 2**N
          NPD2  = NP / 2
          NPD4  = NP / 4
          NPD16 = NP / 16
          NPM1  = NP - 1
C                                      build reciprical sine table
          DT    = PI / FLOAT ( NP )
          DO 10 J = 1, NPM1
             ARG = DT * J
             SINES ( J ) = (0.5D0 / SIN ( ARG ))
   10     CONTINUE
C                                 construct bit reversed subscript table
          J1 = 0
          DO 30 J = 1, NPD2 - 1
             J2 = NPD2
   20        CONTINUE
             IF ( IAND ( J1, J2 ) .NE. 0 ) THEN
                J1 = IABS ( J1 - J2 )
                J2 = J2 / 2
                GO TO 20
             ENDIF
             J1 = J1 + J2
             JINDX ( J ) = J1
   30     CONTINUE 
C
C                          form the trig tables for the radix-8 passes; 
C                          tables are stored in bit reversed order.
          J1 = 0
          DO 50 J = 1, NPD16 - 1
             J2 = NPD16
   40        CONTINUE
             IF ( IAND ( J1, J2 ) .NE. 0 ) THEN
                J1 = IABS ( J1 - J2 )
                J2 = J2 / 2
                GO TO 40
             ENDIF
             J1  = J1 + J2
             ARG = DT * FLOAT (J1)
             COSTBL ( J ) =  COS ( ARG )
             SINTBL ( J ) = -SIN ( ARG )
   50     CONTINUE
C
       ENDIF
C
C                     ***  form the input Fourier coefficients ***
C
C                                                       sine transform
          B ( 1 ) = -2. * X ( 1 )
          B ( 2 ) =  2. * X ( NPM1 )
          J1 = 0
          DO 110 J = 3, NPM1, 2
             J1 = J1 + 1
             J2 = JINDX ( J1 )
             B ( J )     = X ( J2 - 1 ) - X ( J2 + 1 )
             B ( J + 1 ) = X ( NP-J2 )
  110     CONTINUE
C
C                   **************************************  
C                   *                                    *
C                   *    Begin Fast Fourier Synthesis    *
C                   *                                    *
C                   **************************************
C
       IF ( N8 .NE. 0 ) THEN
C                                                radix-8 iterations
          INTT = 1
          NT = NPD16
          DO 130 J = 1, N8
             J1 =  1 + INTT
             J2 = J1 + INTT
             J3 = J2 + INTT
             J4 = J3 + INTT
             J5 = J4 + INTT
             J6 = J5 + INTT
             J7 = J6 + INTT
C***
             CALL R8SYN (INTT, NT, COSTBL, SINTBL, B(1), B(J1), B(J2), 
     *                       B(J3), B(J4), B(J5), B(J6), B(J7)  )
C***
             NT = NT / 8
             INTT = 8 * INTT
  130     CONTINUE
       ENDIF
C                                                radix-4 iteration
       IF ( N4 .GT. 0 ) THEN
          J1 = NPD4
          J2 = 2*NPD4
          J3 = 3*NPD4
          DO 140 J = 1, NPD4   
             T0 = B(J) + B(J + J1) 
             T1 = B(J) - B(J + J1) 
             T2 = 2. * B(J + J2)
             T3 = 2. * B(J + J3)
             B(J)      = T0 + T2    
             B(J + J2) = T0 - T2    
             B(J + J1) = T1 + T3    
             B(J + J3) = T1 - T3    
  140     CONTINUE             
C
       ELSE IF ( N4 .EQ. 0 ) THEN
C                                                radix-2 iteration    
          K = NPD2
          DO 150 J = 1, NPD2
             K    = K + 1
             T    = B(J) + B (K)       
             B(K) = B(J) - B (K)
             B(J) = T 
  150     CONTINUE 
       ENDIF
C
C                   ************************  
C                   *                      *
C                   *    Form Transform    *
C                   *                      *
C                   ************************
C
C
C                                                sine transform
          J1 = NP          
          DO 160 J = 1, NPM1
             X(J) = .25*(( B(J+1) + B(J1)) * SINES(J) - B(J+1) + B(J1))
               J1 = J1 - 1
  160     CONTINUE

       RETURN
       END
C
C
       SUBROUTINE R8SYN ( INTT, NT, COSTBL, SINTBL, B0, B1, B2, B3, 
     *                    B4, B5, B6, B7 )
C
C***********************************************************************
C
C  PURPOSE:    Radix-8 synthesis subroutine used by mixed radix driver.
C
C
C
C***********************************************************************
C
C
       DIMENSION   COSTBL(*), SINTBL(*)
       DIMENSION B0(*), B1(*), B2(*), B3(*), B4(*), B5(*), B6(*), B7(*)
C
C
C            ///     Local variables     ///
C
C
       DOUBLE PRECISION C1, C2, C3, C4, C5, C6, C7
       DOUBLE PRECISION S1, S2, S3, S4, S5, S6, S7
       DOUBLE PRECISION CPI4, CPI8, R2, SPI8
C
       SAVE  CPI4, CPI8,  R2, SPI8 
C
       DATA  R2    / 1.41421 35623 7310D+0 /, 
     *       CPI4  / 0.70710 67811 8655D+0 /,
     *       CPI8  / 0.92387 95325 1129D+0 /, 
     *       SPI8  / 0.38268 34323 6509D+0 /
C
C
C----------------------------------------------------------------------+
C
      JT = 0
      JL = 2
      JR = 2
      JI = 3
      INT8 = 8 * INTT
C
      DO 60 K = 1, INTT
        T0 = B0(K) + B1(K)
        T1 = B0(K) - B1(K)
        T2 = B2(K) + B2(K)
        T3 = B3(K) + B3(K)
        T4 = B4(K) + B6(K)
        T5 = B4(K) - B6(K)
        T6 = B7(K) - B5(K)
        T7 = B7(K) + B5(K)
        T8 = R2 * (T7 - T5)
        T5 = R2 * (T7 + T5)
        TT0 = T0 + T2
        T2  = T0 - T2
        TT1 = T1 + T3
        T3  = T1 - T3
        T4  = T4 + T4
        T6  = T6 + T6
C
        B0(K) = TT0 + T4
        B4(K) = TT0 - T4
        B1(K) = TT1 + T5
        B5(K) = TT1 - T5
        B2(K) = T2 + T6
        B6(K) = T2 - T6
        B3(K) = T3 + T8
        B7(K) = T3 - T8
   60 CONTINUE
C
       IF ( NT .EQ. 0 )                RETURN
C
       K0 = INT8 + 1
       KLAST = INT8 + INTT
C
       DO 70 K = K0, KLAST
          T1 = B0(K) + B6(K)
          T3 = B0(K) - B6(K)
          T2 = B7(K) - B1(K)
          T4 = B7(K) + B1(K)
          T5 = B2(K) + B4(K)
          T7 = B2(K) - B4(K)
          T6 = B5(K) - B3(K)
          T8 = B5(K) + B3(K)
C
          B0(K) = (T1 + T5) + (T1 + T5)
          B4(K) = (T2 + T6) + (T2 + T6)
          T5    = T1 - T5
          T6    = T2 - T6
          B2(K) = R2 * (T6 + T5)
          B6(K) = R2 * (T6 - T5)
          T1    = T3 * CPI8 + T4 * SPI8
          T2    = T4 * CPI8 - T3 * SPI8
          T3    = T8 * CPI8 - T7 * SPI8
          T4    =  - T7 * CPI8 - T8 * SPI8
          B1(K) = (T1 + T3) + (T1 + T3)
          B5(K) = (T2 + T4) + (T2 + T4)
          T3    = T1 - T3
          T4    = T2 - T4
          B3(K) = R2 * (T4 + T3)
          B7(K) = R2 * (T4 - T3)
   70  CONTINUE
C
      DO 90 JT = 1, NT-1
       C1 = COSTBL(JT)
       S1 = SINTBL(JT)
       C2 = C1 * C1 - S1 * S1
       S2 = C1 * S1 + C1 * S1
       C3 = C1 * C2 - S1 * S2
       S3 = C2 * S1 + S2 * C1
       C4 = C2 * C2 - S2 * S2
       S4 = C2 * S2 + C2 * S2
       C5 = C2 * C3 - S2 * S3
       S5 = C3 * S2 + S3 * C2
       C6 = C3 * C3 - S3 * S3
       S6 = C3 * S3 + C3 * S3
       C7 = C3 * C4 - S3 * S4
       S7 = C4 * S3 + S4 * C3
C
       K  = JI * INT8
       J0 = JR * INT8 + 1
       JLAST = J0 + INTT - 1
C
       DO 80 J = J0, JLAST
C
          K   = K + 1
          TR0 = B0(J) + B6(K)
          TR1 = B0(J) - B6(K)
          TI0 = B7(K) - B1(J)
          TI1 = B7(K) + B1(J)
          TR2 = B4(K) + B2(J)
          TI3 = B4(K) - B2(J)
          TI2 = B5(K) - B3(J)
          TR3 = B5(K) + B3(J)
          TR4 = B4(J) + B2(K)
          T0  = B4(J) - B2(K)
          TI4 = B3(K) - B5(J)
          T1  = B3(K) + B5(J)
          TR5 = CPI4 * (T1 + T0)
          TI5 = CPI4 * (T1 - T0)
          TR6 = B6(J) + B0(K)
          T0  = B6(J) - B0(K)
          TI6 = B1(K) - B7(J)
          T1  = B1(K) + B7(J)
          TR7  =  - CPI4 * (T0 - T1)
          TI7  =  - CPI4 * (T0 + T1)
          T0   = TR0 + TR2
          TR2  = TR0 - TR2
          T1   = TI0 + TI2
          TI2  = TI0 - TI2
          T2   = TR1 + TR3
          TR3  = TR1 - TR3
          T3   = TI1 + TI3
          TI3  = TI1 - TI3
          T5   = TI4 + TI6
          TTR6 = TI4 - TI6
          TI6  = TR6 - TR4
          T4   = TR4 + TR6
          T7   = TI5 + TI7
          TTR7 = TI5 - TI7
          TI7  = TR7 - TR5
          T6   = TR5 + TR7
C
          B0(J) = T0 + T4
          B0(K) = T1 + T5
          B4(J) = C4 * (T0 - T4)    - S4 * (T1 - T5)
          B4(K) = C4 * (T1 - T5)    + S4 * (T0 - T4)
C
          B1(J) = C1 * (T2 + T6)    - S1 * (T3 + T7)
          B1(K) = C1 * (T3 + T7)    + S1 * (T2 + T6)
          B5(J) = C5 * (T2 - T6)    - S5 * (T3 - T7)
          B5(K) = C5 * (T3 - T7)    + S5 * (T2 - T6)
C
          B2(J) = C2 * (TR2 + TTR6) - S2 * (TI2 + TI6)
          B2(K) = C2 * (TI2 + TI6)  + S2 * (TR2 + TTR6)
          B6(J) = C6 * (TR2 - TTR6) - S6 * (TI2 - TI6)
          B6(K) = C6 * (TI2 - TI6)  + S6 * (TR2 - TTR6)
C
          B3(J) = C3 * (TR3 + TTR7) - S3 * (TI3 + TI7)
          B3(K) = C3 * (TI3 + TI7)  + S3 * (TR3 + TTR7)
          B7(J) = C7 * (TR3 - TTR7) - S7 * (TI3 - TI7)
          B7(K) = C7 * (TI3 - TI7)  + S7 * (TR3 - TTR7)
C
   80  CONTINUE
C
       JR = JR + 2
       JI = JI - 2
       IF ( JI .GT. JL) GOTO 90
       JI = JR + JR - 1
       JL = JR
   90 CONTINUE
C
      RETURN
      END
