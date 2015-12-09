c     vfsa.f sa using very fast simulated annealing
c     This SUBROUTINE is adapted by Sen and inspired by Ingber
c     Peter Gerstoft 96

      SUBROUTINE vfsa
      USE global
      INCLUDE 'comopt.h'
c     PARAMETER (mpar=2)
c     COMMON /amod/nparam,xmin(mpar),xmax(mpar),dx(mpar),t0(mpar),
c     .               nx(mpar)

      
      INTEGER i,jmov,jtemp,jj,jruns,maxtemps,nmov
      REAL temp,tmp(mpar),t1
      REAL e,ebest,arg
      REAL eold,fold(mpar),fbest(mpar),ftemp(mpar)
c     functions
      REAL sample,walk,ran2
      REAL pde,decay
      INTEGER ibest,iforwtotal
      iforwtotal=0

c     simple plotting of sampling & histogram of RESULT
c     FUNCTION is func in SUBROUTINE error.f


      WRITE(*,*)' --- Using very fast simulated annealing --- '  
c     get the tmperature 

c     WRITE(*,*)'enter starting temperature for x & y '
c     READ(*,*)t0(1),t0(2)

c     WRITE(*,*)'enter starting temp for use in acceptance rule'
c     READ(*,*)temp0

c     get the # of temperatures & the number per temperature

c     WRITE(*,*)'# of temperatures, # of random walks/temperature'
c     READ(*,*)maxtemps,nmov

c     WRITE(*,*)'no of runs ????'
c     READ(*,*)nruns

c     define the FUNCTION range in x

c     DO i=1,nparam
c     nx(i)=(xmax(i)-xmin(i))/dx(i)+1
c     END DO

c-----------------------------------------------------------------

c     VFSA

c     plot the FUNCTION, the temperature for the probability is templt
c     nruns=1
      DO jruns=1,npop
         iforwtotal=0
c     initial random guess at xmod, ymod
c     & get emod, the initial error

         DO i=1,nparm
            f(i)=
     &           sample((fmin(i)+fmax(i))/2,df(i),ndigit(i),
     &           fmin(i),fmax(i))
         ENDDO
         CALL setmodelreal(f)
         CALL forw2
         iforwtotal=iforwtotal+1
         CALL cost(e)
         EBEST=E
         WRITE(*,*)' starting parameters',(f(i),i=1,nparm)
         WRITE(*,*)'  and energy',e
         WRITE(97,*)e,f(1),f(2),f(1),f(2),fbest(1),fbest(2)
         
         decay = 1.0
         nmov=nparm
         maxtemps=niter
         WRITE(prtfil,*)' decay    ',decay
         WRITE(prtfil,*)' temp     ', temp0
         WRITE(prtfil,*)' nmov     ', nmov
         WRITE(prtfil,*)' maxtemps ', maxtemps
c     WRITE(prtfil,*)' Niter    ',Niter
         IF (temp0.LE.0) THEN
            WRITE(*,*)'temperature less than zero',temp0
            STOP
         ENDIF

         DO jtemp=1,maxtemps

c     initialize tempratures
            temp   = temp0*EXP(-decay*(jtemp-1)**(1./nparm))
            DO i=1,nparm 
               tmp(i) = temp0*EXP(-decay*(jtemp-1)**(1./nparm))
            ENDDO

            DO jmov=1,nmov      ! metropolis

c     try a discrete random walk for xmod & ymod
c     & get etrial, the NEW error
               
               eold=e
               DO i=1,nparm
                  fold(i)=f(i)
               ENDDO
               IF (isubopt(4).EQ.0) THEN
                  f(jmov)=walk(fold(jmov),df(jmov)
     1                 ,fmin(jmov),fmax(jmov),tmp(jmov))
               ELSE
                  DO i=1,nparm
                     f(i)=walk(fold(i),df(i),fmin(i),fmax(i),tmp(i))
                  ENDDO
               ENDIF
c     WRITE(*,*)'parameters',(f(i),i=1,nparm)

               CALL setmodelreal(f)
               CALL forw2
               CALL cost(e)
               iforwtotal=iforwtotal+1
c     WRITE(*,*)'new energy',e
               ftemp(1)=f(1)
               ftemp(2)=f(2)
               
c     IF e < eold, accept the NEW model
c     but,IF e > eold compute Boltzann 
               IF(e.GT.eold) THEN

c     ? accept even though etrial > emod
                  arg=(e-eold)/temp
                  IF(arg .GT. 1.e6) THEN
                     pde = 0.001
                  ELSE
                     pde = EXP(-arg)
                  ENDIF
                  
                  IF(pde.LT.ran2(1)) THEN
c     move rejected;    revert to old solution
c     WRITE(*,*)' move rejected'
                     E=EOLD 
                     DO i=1,nparm
                        F(i)=fOLD(i)
                     ENDDO
                  ENDIF
               END IF
c     1000		CONTINUE
               IF(E.LT.EBEST)THEN
                  WRITE(60,'(1x,I7,I7,G12.5,100f12.2)')
     &                 jruns,iforwtotal,e,(f(jj),jj=1,nparm)
                  WRITE(*,*)'   **** NEW BEST ENERGY ****'
                  EBEST=E
                  IBEST=jtemp
                  DO 9 JJ=1,Nparm
                     FBEST(JJ)=F(JJ)
 9                CONTINUE
                  WRITE(prtfil,*)jtemp,Ebest
                  WRITE(prtfil,*)jtemp,(f(jj),jj=1,nparm)
                  WRITE(*,*)jtemp,temp,iforwtotal,ebest
                  WRITE(*,*)(f(jj),jj=1,nparm)
               END IF
               WRITE(97,*)e,ftemp(1),ftemp(2),f(1),f(2),
     &              fbest(1),fbest(2)

               IF ((jtemp.EQ.1).AND.(jruns.EQ.1).AND.
     &              (jmov.EQ.nmov)) THEN
                  CALL rdtime(t1)
                  WRITE(prtfil,320) t1*(maxtemps*nmov*npop)/nmov
                  WRITE(*,320) t1*(maxtemps*nmov*npop)/nmov
 320              FORMAT(/1H ,' >>> Estimated serial Inversion time: ',
     1                 F12.1,' secs.')
               ENDIF

            ENDDO               ! metropolis step
         ENDDO                  ! temperature loop


         WRITE(*,*)' Best energy & parameters'
         WRITE(*,*)EBEST,IBEST
         DO 12 Jj=1,Nparm
            WRITE(*,*)Jj,FBEST(Jj)
            WRITE(prtfil,*)Jj,FBEST(jJ)
 12      CONTINUE

         WRITE(*,*)'number of forward modelling calls',iforwtotal
         WRITE(*,*)' End of vfsa'

      ENDDO                     ! 

      END	

      REAL FUNCTION sample(xmod,delx,nx,xmin,xmax)
      REAL xmod,delx,xmin,xmax
      INTEGER nx,ix
      REAL ran2

c     get an integerized random sample from xmod
c     this could be simplified or changed
c     note try agin IF at xmin or xmax
      
 1    ix = 2.*(ran2(1)-.5)*nx
      sample = xmod + ix*delx
      IF(sample.GT.xmax) GOTO 1
      IF(sample.LT.xmin) GOTO 1
      RETURN
      END


      REAL FUNCTION walk(xmod,dx,xmin,xmax,tmp)
      REAL xmod,dx,xmin,xmax,tmp
c---  local variables
      REAL ran2
      REAL arand,ayy,dif,xmod1,yy,pwr
      INTEGER ntry

c     generate a NEW state following Ingber ----

      ntry = 1
 123  CONTINUE

      arand = ran2(1)
      ayy = 0.0
      dif = arand - 0.5
      IF (dif.LT.0.0) ayy = -1.0
      IF (dif.GE.0.0) ayy = +1.0

      pwr = ABS(2*arand-1.)
      yy = ayy*tmp*( (1+1/tmp)**pwr - 1.)
      xmod1 = xmod + yy*(xmax-xmin)

      IF (xmod1.LT.xmin.OR.xmod1.GT.xmax) THEN
         ntry = ntry + 1
         IF( ntry.LT.100) go to 123
      END IF

      walk = xmod1
      RETURN
      END




	

	
	

