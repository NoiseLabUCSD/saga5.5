c     
c***************************************************************
c     tHis is for all normal machines!!!!      
      SUBROUTINE  crossover(n1,n2,jx,nbit)    
c---  crossover of 2 set of genes
      INTEGER n1,n2,jx,nbit
      INTEGER len,ntemp         ! local variables
      len=nbit-jx
      ntemp=n1
      CALL MVBITS(n2   ,jx,len,n1,jx)
      CALL MVBITS(ntemp,jx,len,n2,jx)
      END
c     
c     this version is for SUN OSF
c     
c     SUBROUTINE  crossover(n1,n2,jx,nbit)    
c---  crossover of 2 set of genes
c     INTEGER n1,n2,jx,nbit
c     INTEGER len,ntemp ! local variables
c     INTEGER imask,n1l,n1u,n2l,n2u
c     len=nbit-jx
c     ntemp=n1
c     imask=2**jx-1
c     WRITE(*,*)'n1,n2,jx',n1,n2,jx
c     n1l=and(imask,n1)
c     n2l=and(imask,n2)
c     n1u=n1-n1l
c     n2u=n2-n2l
c     n1=n1l+n2u
c     n2=n2l+n1u
ccc   WRITE(*,*)'n1,n2,n1l,n1u,n2l,n2u',n1,n2,n1l,n1u,n2l,n2u
ccc   CALL MVBITS(n2   ,jx,len,n1,jx)
ccc   CALL MVBITS(ntemp,jx,len,n2,jx)
c     END
c     
c************************************************
c     
      INTEGER FUNCTION mutate(n1,jm)
      INTEGER n1,jm
      IF (BTEST(n1,jm)) THEN
         mutate=IBCLR(n1,jm)
      ELSE
         mutate=IBSET(n1,jm)
      ENDIF
      END

c**************************************
c     gray incoding the parents    
      SUBROUTINE grayincode(qpar)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'graycode.h'
      INTEGER iq,jn,qpar       
c     
      DO iq=1,qpar
         DO jn=1,nparm
c     WRITE(prtfil,*)'bin ', parents(jn,iq)
            parents(jn,iq)=bin2g(parents(jn,iq))    
c     WRITE(prtfil,*)'gray', parents(jn,iq)
         ENDDO                  !jn
      ENDDO                     !iq     
      RETURN
      END      
c     gray decoding the parents    
      SUBROUTINE graydecode(qpar)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'graycode.h'
      INTEGER iq,jn,qpar       
c     
      DO iq=1,qpar
         DO jn=1,nparm
c     WRITE(prtfil,*)'gray', parents(jn,iq)
            parents(jn,iq)=gray2b(parents(jn,iq))    
c     WRITE(prtfil,*)'back', parents(jn,iq)
         ENDDO                  !jn
      ENDDO                     !iq     
      RETURN
      END      
c**************************************
c---  test for crossover(qpar)
      SUBROUTINE cross(qpar)
      USE global
      INCLUDE 'comopt.h'
      INTEGER iq,jn,jx,qpar       
      REAL ran2
      IF (iopt(29).EQ.0) THEN
c     multiple point crossover
         DO iq=1,qpar,2
            DO jn=1,nparm
               IF (ran2(1).LT.px) THEN
                  jx=ifix(ran2(1)*(nbit(jn)-1))+1
c     WRITE(prtfil,*)'crossover,iq,jx',iq,jx
                  CALL crossover(parents(jn,iq),parents(jn,iq+1),jx,
     &                 nbit(jn))    
               ENDIF                
            ENDDO               !jn
         ENDDO                  !iq           
      ELSE
c     single crossover
         DO iq=1,qpar,2
            jn=ifix(ran2(1)*(nparm-1))+1
            jx=ifix(ran2(1)*(nbit(jn)-1))+1
c     WRITE(*,*)'crossover,iq,jx,jn',iq,jx,jn
            CALL crossover(parents(jn,iq),parents(jn,iq+1),jx,
     &           nbit(jn))    
         ENDDO                  !iq           
      ENDIF            
      END
c*******************************************
      SUBROUTINE mutation(qpar)
      USE global
      INCLUDE 'comopt.h'
      INTEGER iq,jn,ib,qpar
      INTEGER mutate
      REAL ran2
c---  test for mutations
      DO iq=1,qpar
         DO jn=1,nparm
            DO ib=0,nbit(jn)-1
c--   ib is the mutation bit, the first bit is zero
               IF (ran2(1).LT.pm) THEN
c     WRITE(prtfil,*)' MUTATION'
                  parents(jn,iq)=mutate(parents(jn,iq),ib)  
               ENDIF 
            ENDDO               ! ib
         ENDDO                  !jn
      ENDDO                     !iq  
      END
c**************************************************************
c     
      SUBROUTINE sortfitPost(fit,fittmp)
c     sorts according to fit the model-array. Both fit and model are returned 
c     sorted. On RETURN all models does ONLY appear once.
c     
      USE global
      INCLUDE 'comopt.h'
      REAL fit(mq_post),fittmp(mq_post)
      INTEGER j,i,l
      INTEGER parentsPost(mpar,mq_post) ! integer discretization of the models
      INTEGER idtot,iact,ipoint(mq_post)

      idtot=1
c     
c     order fitness in increasing order according to POINTER ipoint
c     
      ipoint(1)=1
      DO 10 I=2,q
         DO 5 j=1,idtot
            IF (fit(ipoint(J)).GT.fit(i)) THEN
               iact=j
               idtot=idtot+1
               DO 3 L=idtot,iact+1,-1
 3                ipoint(l)=ipoint(l-1)
                  ipoint(iact)=i
                  GOTO 10
               ELSE IF (j.EQ.idtot) THEN
                  idtot=idtot+1
                  ipoint(idtot)=i
                  GOTO 10
               END IF
 5          CONTINUE
 10      CONTINUE
c---  copy fit and model.
         DO i=1,q
            fittmp(i)=fit((i))
            DO j=1,nparm
               parentsPost(j,i)=model(j,i)
            ENDDO
         ENDDO 
c---  copy into fit and model using the constructed POINTER array.
         DO i=1,q
            fit(i)=fittmp(ipoint(i))
            DO j=1,nparm
               model(j,i)=parentsPost(j,ipoint(i))
            ENDDO
         ENDDO      
         END
c**************************************************************
c     
      SUBROUTINE sortfit(fit,fittmp)
c     sorts according to fit the model-array. Both fit and model are returned 
c     sorted.
c     
      USE global
      INCLUDE 'comopt.h'
      REAL fit(*),fittmp(*)
      INTEGER j,i,l
      INTEGER idtot,iact,ipoint(mq)

      idtot=1
c     
c     order fitness in increasing order according to POINTER ipoint
c     
      ipoint(1)=1
      DO 10 I=2,q
         DO 5 j=1,idtot
            IF (fit(ipoint(J)).GT.fit(i)) THEN
               iact=j
               idtot=idtot+1
               DO 3 L=idtot,iact+1,-1
 3                ipoint(l)=ipoint(l-1)
                  ipoint(iact)=i
                  GOTO 10
               ELSE IF (j.EQ.idtot) THEN
                  idtot=idtot+1
                  ipoint(idtot)=i
                  GOTO 10
               END IF
 5          CONTINUE
 10      CONTINUE
c---  copy fit and model.
         DO i=1,q
            fittmp(i)=fit((i))
            DO j=1,nparm
               parents(j,i)=model(j,i)
            ENDDO
         ENDDO 
c---  copy into fit and model using the constructed POINTER array.
         DO i=1,q
            fit(i)=fittmp(ipoint(i))
            DO j=1,nparm
               model(j,i)=parents(j,ipoint(i))
            ENDDO
         ENDDO      
         END

c*****************************************************************
      SUBROUTINE insertmodel(qpar,irep)
      USE global
      INCLUDE 'comopt.h'
      INTEGER qpar
      INTEGER iq,iiq,j,irep
      LOGICAL lcomp(mq)

c---- reset
      DO iq=1,q
         lcomp(iq)=.TRUE.
      ENDDO
c     
c---  test IF each PARAMETER in a model vector is within the bounds
      DO iq=1,qpar
         DO j=1,nparm
            IF (parents(j,iq).GT.ndigit(j)-1) THEN
c     WRITE(*,*)'parents,before',parents(j,iq),j
c     parents(j,iq)=IBCLR(parents(j,iq),nbit(j)-1)
c     WRITE(*,*)'parents,after ',parents(j,iq),j
               lcomp(iq)=.FALSE.
            ENDIF
         ENDDO
      ENDDO
c     
c     
c---  test does the parent exist in the last  part of the model vector ?
      DO iq=2,qpar
         IF (lcomp(iq)) THEN
            DO iiq=1,iq-1
               DO j=1,nparm
                  IF (parents(j,iiq).NE. parents(j,iq)) GOTO 11 ! identical ?
               ENDDO
               lcomp(iq)=.FALSE.
               GOTO 21
 11            CONTINUE
            ENDDO
 21         CONTINUE
         ENDIF
      ENDDO
c     
c     test does the parent exist in the  model vector ?
      DO iq=1,qpar
         IF (lcomp(iq)) THEN
            DO iiq=1,q
               DO j=1,nparm
                  IF (model(j,iiq).NE. parents(j,iq)) GOTO 10 ! identical ?
               ENDDO
               lcomp(iq)=.FALSE.
               GOTO 20
 10            CONTINUE
            ENDDO
 20         CONTINUE
         ENDIF
      ENDDO
c     
c     test does the child exist in any previous  model vector ?
      DO iq=1,qpar
         IF (lcomp(iq)) THEN
            DO iiq=iallmodel,1,-1
               DO j=1,nparm
                  IF (allmodel(j,iiq).NE. parents(j,iq)) GOTO 12 ! identical ?
               ENDDO
               lcomp(iq)=.FALSE.
c     WRITE(*,*)'rejecting child',iq,' was in model',iiq
               GOTO 22
 12            CONTINUE
            ENDDO
 22         CONTINUE
         ENDIF
      ENDDO

c     ok test now put a NEW model
      irep=0
      DO iq=1,qpar
         IF (lcomp(iq)) THEN
c     iallmodel=iallmodel+1
c     allfit(iallmodel)=fit(iq)
            DO j=1,nparm
               model(j,q-irep)= parents(j,iq) ! move
c     allmodel(j,iallmodel)= parents(j,iq)  ! move
            ENDDO
            irep=irep+1
         ENDIF
      ENDDO
c     WRITE(7,*)'Models replaced:',irep
c     WRITE(*,*)'iallmodel:',iallmodel
c     WRITE(90,*)'iallmodel:',iallmodel


c     WRITE(90,*)(allmodel(1,iq), iq=1,iallmodel)
      END

c*****************************************************************
      SUBROUTINE elimiatemodels(qpar,fit)
      USE global
      INCLUDE 'comopt.h'
      INTEGER qpar
      INTEGER iq,iiq,j,irep
      LOGICAL lcomp(mq_post)
      REAL fit(mq_post)

c---- reset
      DO iq=1,qpar
         lcomp(iq)=.TRUE.
      ENDDO
c     
c     
c---  test does the parent exist in the last  part of the model vector ?
      DO iq=2,qpar
         IF (lcomp(iq)) THEN
            DO iiq=1,iq-1
               DO j=1,nparm
                  IF (model(j,iiq).NE. model(j,iq)) GOTO 11 ! identical ?
               ENDDO
               lcomp(iq)=.FALSE.
               GOTO 21
 11            CONTINUE
            ENDDO
 21         CONTINUE
         ENDIF
      ENDDO


c     ok test now put a NEW model
      irep=0
      DO iq=1,qpar
         IF (lcomp(iq)) THEN
            irep=irep+1
            fit(irep)=fit(iq)
            DO j=1,nparm
               model(j,irep)= model(j,iq) ! move
            ENDDO
         ENDIF
      ENDDO
      q=irep
      WRITE(*,*)' After elimination q=',q
      END

c*******************************************************
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
C*                  ONLY 714025 different values.                     *C
C*                  See Numerical Recipies, chapter 7.1.              *C
C* See also       :                                                   *C
C* EXTERNAL calls : NONE                                              *C
C* Cautions       : ...............................................   *C
C* Latest revision 901122 by Theo Zimmermann                          *C
C*--------------------------------------------------------------------*C
C* HISTORY:                                                           *C
C* =======                                                            *C
C     *                                                                    *C
C     * Recording of modifications, changes, updates, extentions and       *C
C     * releases of the MODULE stating the most recent first.              *C
C     *                                                                    *C
C     * Date.  Who.                    What.                               *C
C     *        .                       .                                   *C
C     * 901122 Theo Zimmermann         Adapted from Numerical Recipies.    *C
C     *                                Changes: SAVE statements            *C
C     *                                and declarations added;             *C
C     *                                ISEED introduced and IDUM SAVEd     *C
C     *                                such that USE of RAN2 is exactly    *C
C     *                                like RAN1, namely ISEED is ONLY     *C
C     *                                input variable.                     *C
C     *                                                                    *C
C     *                                                                    *C
C     * EXTERNAL MODULES:                                                  *C
C     * ================                                                   *C
C     *                                                                    *C
C     * CONSTANTS:                                                         *C
C     * =========                                                          *C
C     *                                                                    *C
C     * XXXXXX         : ...........................................       *C
C     *                                                                    *C
C     *                                                                    *C
C     * GLOBAL VARIABLES:                                                  *C
C     * ================                                                   *C
C     *                                                                    *C
C     * XXXXXX         : ...........................................       *C
C     *                                                                    *C
C     *                                                                    *C
C     * LOCAL VARIABLES:                                                   *C
C     * ===============                                                    *C
C     *                                                                    *C
C     * XXXXXX         : ...........................................       *C
C     *                                                                    *C
C     *       1         2         3         4         5         6         7*C
C     *34567890123456789012345678901234567890123456789012345678901234567890*C

      REAL FUNCTION RAN2(ISEED)
      IMPLICIT       NONE
C     *                           Formal parameters:                       *C
      INTEGER      ISEED
C     *                           EXTERNAL FUNCTION declarations:          *C
C     *                           Local variables:                         *C
      INTEGER      IA, IC, M, IFF, J, IR(97), IY, IDUM
      REAL         RM
C     *                           Parameters and constants:                *C
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112E-6)
C     *                           Global variables (COMMON areas):         *C
C     *                           Initialized local variables:             *C
C     *                           Global parameters and constants:         *C
C     *                           Declaration of static local variables:   *C
      INTEGER ii
      DATA IFF /0/
C     
C     
C     
      SAVE           IFF, IR, IDUM, ii, iy
C     *                                                                    *C
C     *                           the following 2 statements               *C
C     *                           have been changed                        *C
c     WRITE(*,*)'iseed',iseed
      IF(ISEED.LT.0.OR.IFF.EQ.0)THEN
         ii=0
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
c     WRITE(*,*)'J,IY,M,IFf,idum,ii'
c     WRITE(*,*)J,IY,M,IFf,idum,ii
      IF(J.GT.97.OR.J.LT.1)PAUSE
      IY=IR(J)
      RAN2=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      ii=ii+1
c     WRITE(*,*)'randomnumber',ii,ran2 
      RETURN
      END
