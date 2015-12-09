C*--------------------------------------------------------------------*C
C* Title          : RAN2                                              *C
C* Purpose        : Random number generator                           *C
C*                                                                    *C
C* Computer       : ALLIANT                                           *C
C* Compiler       : fortran                                           *C
C*                                                                    *C
C* Usage          : CALL RAN2(ISEED)                                  *C
C* Formal parm.   : I  INTEGER       ISEED                            *C
C* Description    : Generates uniform random deviates in [0,1)        *C
C*                  Properties: no sequential correlations, but       *C
C*                  only 714025 different values.                     *C
C*                  See Numerical Recipies, chapter 7.1.              *C
C* See also       :                                                   *C
C* External calls : none                                              *C
C* Cautions       : ...............................................   *C
C* Latest revision 901122 by Theo Zimmermann                          *C
C*--------------------------------------------------------------------*C
C* HISTORY:                                                           *C
C* =======                                                            *C
C*                                                                    *C
C* Recording of modifications, changes, updates, extentions and       *C
C* releases of the module stating the most recent first.              *C
C*                                                                    *C
C* Date.  Who.                    What.                               *C
C*        .                       .                                   *C
C* 901122 Theo Zimmermann         Adapted from Numerical Recipies.    *C
C*                                Changes: SAVE statements            *C
C*                                and declarations added;             *C
C*                                ISEED introduced and IDUM SAVEd     *C
C*                                such that use of RAN2 is exactly    *C
C*                                like RAN1, namely ISEED is only     *C
C*                                input variable.                     *C
C*                                                                    *C
C*                                                                    *C
C* EXTERNAL MODULES:                                                  *C
C* ================                                                   *C
C*                                                                    *C
C* CONSTANTS:                                                         *C
C* =========                                                          *C
C*                                                                    *C
C* XXXXXX         : ...........................................       *C
C*                                                                    *C
C*                                                                    *C
C* GLOBAL VARIABLES:                                                  *C
C* ================                                                   *C
C*                                                                    *C
C* XXXXXX         : ...........................................       *C
C*                                                                    *C
C*                                                                    *C
C* LOCAL VARIABLES:                                                   *C
C* ===============                                                    *C
C*                                                                    *C
C* XXXXXX         : ...........................................       *C
C*                                                                    *C
C*       1         2         3         4         5         6         7*C
C*34567890123456789012345678901234567890123456789012345678901234567890*C
      REAL FUNCTION RAN2(ISEED)
      IMPLICIT       NONE
C*                           Formal parameters:                       *C
      INTEGER      ISEED
C*                           External function declarations:          *C
C*                           Local variables:                         *C
      INTEGER      IA, IC, M, IFF, J, IR(97), IY, IDUM
      REAL         RM
C*                           Parameters and constants:                *C
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112E-6)
C*                           Global variables (common areas):         *C
C*                           Initialized local variables:             *C
      DATA IFF /0/
C*                           Global parameters and constants:         *C
C*                           Declaration of static local variables:   *C
      integer ii
C
C
C
      SAVE           IFF, IR, IDUM, ii, iy
C*                                                                    *C
C*                           the following 2 statements               *C
C*                           have been changed                        *C
       write(*,*)'iseed',iseed
      IF(ISEED.LT.0.OR.IFF.EQ.0)THEN
        IDUM = ISEED
        IFF=1
        IDUM=MOD(IC-IDUM,M)
        DO 11 J=1,97
          IDUM=MOD(IA*IDUM+IC,M) 
          IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
      ENDIF
      J=1+(97*IY)/M
c      write(*,*)'J,IY,M,IFf,idum,ii'
c      write(*,*)J,IY,M,IFf,idum,ii
      IF(J.GT.97.OR.J.LT.1)PAUSE
      IY=IR(J)
      RAN2=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      ii=ii+1
c      write(11,*)'randomnumber',ii,ran2 
      RETURN
      END





















