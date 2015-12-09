      subroutine forwardmodel(iopt,mopt)
      integer   mopt,i,iopt(mopt)
      DO i=1,mopt
         iopt(i)=0.
      ENDDO
      iopt(1)=2
      end

c****************************************
      subroutine forwinitgrad
C     
      USE global
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnp.f'
      INCLUDE './oases/comnrd.f'
      INCLUDE './oases/comfip.f'
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE 'comgrad.h'
      Complex wavenograd(np,nrd,mpar)
      common  /wavenogr/wavenograd
c     prognm
      real theta1,theta2,dtheta,xrange_ref
      COMPLEX SLOW
      DIMENSION X(NP2,3),FF(2,NP3)              
      DIMENSION FFS(2,NP),XS(NP2),AKM(100),AKM1(100),AKM2(100)     
      LOGICAL NFLAG,ICONTU,lincoh,nofft,ausamp
      EQUIVALENCE (NREC,ISPACE),(LF,NUMFR)
      EQUIVALENCE (FF(1,1),CFF(1,1)),(FFS(1,1),CFFS(1))
      EQUIVALENCE (X(1,1),CFF(1,1)),(XS(1),CFFS(1))  
      integer ibody      
      common /oasesvar/ IPOW,DELTA,THETA,FOCDEP,LTYP,NFLAG,ICNTIN,
     .     ICONTU,lincoh,nofft,ausamp,
     .     ibody
      real              freqm,cminin,cmaxin,cref,ranref,rm,rmi,rr
      common /oasesgard/freqm,cminin,cmaxin,cref,ranref,rm,rmi,rr
C     c
c     The initialization has been done in the forward-model
c     
c**************************
      entry forw2grad( )

c     
c---- quick fix to obtain gradients
c     
      if (ngrad.eq.0) ngrad=nparm

      write(*,*)'ngrad=',ngrad
      write(*,*)'nparm=',nparm
c     nofft=.false.
      if (iopt(1).eq.2.and.  (icut2-icut1+1).lt.mx) then
         mindex=mx
      else
         mindex=(icut2-icut1+1)
      endif
c     if (mindex*ndep.gt. mobs) then
c     write(*,*)' vmsoastgrad: Mobs is too small; it is contains' 
c     write(*,*)' the wavenumber differentiations for each depth'
c     write(*,*)'  Mobs=',mobs,'  mindex*ndep=',mindex*ndep
c     endif
      if (flagpu.le.0) then     
         call recderiv          ! check if derivatives in source or receiver layers  
c     IF (iopt(3).eq.3) THEN
c     xhelp=1./sqrt(xranges(1))
c     do i=1,nx
c     weight(i)=xhelp/sqrt(xranges(i))
c     enddo
c     else
c     do i=1,nx
c     weight(i)=1.0
c     enddo
c     ENDIF
      endif

      flagpu=flagpu+1
      
      lwascomp=1
c     
c---- Use of EOF
c     
      if (iopt(17).eq.1) then
         call eofval_nocheck
      endif

      CALL INENVI



      if (flagpu.lt.0) WRITE(prtfil,918)
 918  FORMAT(//1H ,'RECEIVER DATA:',//1H ,'  NO. ','        DEPTH  ',
     1     'LAYER','           Z')
C     
C     EQUIDISTANT RECEIVER DEPTHS
C     
      DO 921 JJ=1,IR
         CALL RECEIV(V,NUML,RDC(JJ),LAY(JJ),Z(JJ))                 
         if (flagpu.lt.0) WRITE(prtfil,907) JJ,RDC(JJ),LAY(JJ),Z(JJ)
 921  CONTINUE
c---  for a new source reset the depths
      if (flagpu.lt.0)WRITE(prtfil,908)
 908  FORMAT(//1H ,'SOURCE DATA:',//1H ,'  N0. ','        DEPTH  ',
     1     'LAYER','           ZU             ZL')
      DO 906 I=1,LS
         CALL SOURCE(V,NUML,SDC(I),LAYS(I),ZUS(I),ZLS(I))
 906  CONTINUE
      do iparm=1,nparm
         if (par2phy(iparm).eq.8) then
            par2lay(iparm)=lays(1)
         endif
      enddo
      if (flagpu.lt.0) then
         WRITE(prtfil,907) 1,SDC(1),LAYS(1),ZUS(1),ZLS(1)
         IF (LS.GT.1) THEN
            WRITE(prtfil,907) LS,SDC(LS),LAYS(LS),ZUS(LS),ZLS(LS)
         END IF
 907     FORMAT(1H ,I6,G15.6,I6,2G15.6)
      endif

c---  for analytic gradient we uses pointers between input and forward model.
c     w/o gradient it is a dummy routine in oastsub
      call map2forw(1)


      if (tilt) then
         if (nplots.gt.1) stop 'tilt is only allowed for one rage'
         xrange_ref=xranges(1)
      endif


      IF (AUSAMP) THEN
         AUSAMP=.TRUE.
         cref=0e0
         do 981 l=1,numl
            if (v(L,2).lt.2E4) cref=max(cref,v(L,2))
 981     continue
         RANREF=0
c     >>> max and min ranges
         rm =0e0
         rmi=1e20
         do 988 ii=1,nplots
c     rr=abs(1E3*(r0+(ii-1)*rspace))
            rr=abs(xranges(ii))
            if (tilt) then
               rr=rr+dtilt
            endif
            rm=max(rm,rr)
            if (rr.gt.0e0) rmi=min(rmi,rr)
 988     continue
         if (rmi.gt.1e19) rmi=0e0
c     RANREF=RANREF+RM
         ranref=min(ranref+2*rm,6*rm)
         OFFDBIN=0E0
         if (flagpu.lt.0) then
            WRITE(6,*)
            WRITE(6,*) '>>> AUTOMATIC SAMPLING '
            write(6,*) '    REFERENCE SPEED:',CREF
            WRITE(6,*) '    REFERENCE RANGE:',RANREF
         endif
         FREQ=FREQM
         ibody=0
         CALL AUTSAM(CMININ,CMAXIN,rmi,RANREF,CMIN,CMAX,
     &        NWVNO,ICW1,ICW2,ibody)
c     WRITE(6,*) '    MAX NO. OF WAVENUMBERS:',NWVNO
         ICUT1=1
         ICUT2=NWVNO
      ELSE
         CMIN=CMININ
         CMAX=CMAXIN
      END IF                    ! Atomatic sampling loop

      CALL PINIT1

c---------------------Frequency loop
      do 15 jj=1,nfrq,nsensor
         freq=frq(jj)

         DSQ=2E0*PI*FREQ             
         CSQ=DSQ*DSQ                 
         RDSQ=DSQ

         IF (AUSAMP) THEN
C     ***  AUTOMATIC SAMPLING
            
            IF (IBODY.GT.0) THEN
C     ***   SLOWNESS INTEGRATION
               CALL AUTSAM(CMININ,CMAXIN,rmi,RANREF,CMIN,CMAX,
     &              NWVNO,ICW1,ICW2,ibody)
            ELSE
C     ***   WAVENUMBER INTEGRATION
               CALL AUTSAM(CMININ*FREQ/FREQM,CMAXIN,rmi,RANREF,
     &              CMIN,CMAX,NWVNO,ICW1,ICW2,ibody)
            END IF
            if (flagpu.lt.2) then
               WRITE(6,*)
               WRITE(6,*) '>>> FREQUENCY:',FREQ
               WRITE(6,*) '    CMIN,CMAX:  ',CMIN,CMAX
               WRITE(6,*) '    WAVENUMBERS:',NWVNO,ICW1,ICW2
            endif
            ICUT1=1
            ICUT2=NWVNO
            WK0=RDSQ/CMAX           
            WKMAX=RDSQ/CMIN
            IF (NWVNO.GT.1) THEN
               DLWVNO = ( WKMAX - WK0 ) / ( FLOAT(NWVNO-1) )   
            ELSE
               DLWVNO=1.0
            END IF
            IF (ICNTIN.GT.0) THEN
               OFFDB = 60.0*V(LAYS((LS-1)/2+1),2)*
     &              (1E0/CMIN-1E0/CMAX)/NWVNO
               if (flagpu.lt.2) then
                  WRITE(6,*) 'DEFAULT CONTOUR OFFSET APPLIED,',OFFDB,
     &                 ' dB/wavelength'
               endif
            END IF
            IF (NFLAG) THEN
               WRITE(6,*) 'NEGATIVE SPECTRUM BY SYMMETRY'
            END IF
         ELSE
            IF (.NOT.NFLAG) THEN
               IF (CMAX.GT.0) THEN
                  WK0=RDSQ/CMAX           
               ELSE
                  WK0=RDSQ/CMIN-(NWVNO-1)*DLWVNO
               END IF
            END IF
            IF (IBODY.GT.0) THEN
               WK0=RDSQ/CMAX           
               WKMAX=RDSQ/CMIN
               IF (NWVNO.GT.1) THEN
                  DLWVNO = ( WKMAX - WK0 ) / ( FLOAT(NWVNO-1) )     
               ELSE
                  DLWVNO=1.0
               END IF
            END IF
         END IF

         DLRAN=2E0*PI/(NWVNO*DLWVNO) 
         R1=dlran
         RANMAX=NWVNO*DLRAN
         
         FNIFAC =  SQRT( 2.0/PI ) * 0.5               
         IF (ICDR.EQ.1) THEN
            FNI5=DLWVNO*SQRT(FREQ/V(LAYS((LS-1)/2+1),2))
         ELSE
            FNI5=DLWVNO*FNIFAC
         END IF 
c     
c***  and now set the discrete sampling points
c     
         if (iopt(1).eq.2) then

            do i=1,nx
               xhelp=xranges(i)/dlran
               ixpoints(i)=xhelp          
               xweight(i)=xhelp-ixpoints(i)
c     write(*,*)' from oast',i,xranges(i),xweight(i),ixpoints(i)
               if (xranges(i).gt.0.66*ranmax) then
                  write(*,*)' frequency',freq
                  WRITE(*,360) DLRAN,RANMAX*1E-3
 360              FORMAT(1H ,' ',/1H ,'RANGE STEP:   ',F12.3,' m',
     &                 /1H ,'MAX FFT RANGE:',F12.3,' km')
                  write(*,*)'** ERROR ** range is greater than 0.66',
     &                 ' times the max fft-range, for range',
     &                 xranges(i)       
                  stop
               endif
               if (ixpoints(i).eq.0) then
                  write(*,*)'range to small relative to sampling'
                  write(*,*)'for range',xranges(i)       
                  stop
               endif
            enddo

         endif

         CALL PINIT2                 

         CALL PINIT2                 
         CALL PHASES(LS,FREQ,V,DELTA,THETA,LTYP,FOCDEP)
         CALL CALINTgrad(mindex)
C******************
c     
         DO 120 NREC=1,IR
c     index2=(nrec-1)*mindex
c     if (iopt(1).eq.1) then
c     index=(nrec-1)*(icut2-icut1+1)
c     else
c     index=(nrec-1)*(nx)
c     endif
            do 120 isensor=1,nsensor ! number of vector measurements 
               if (iopt(1).eq.1) then
                  index=(nrec-1+(jj-1)*ir)*(icut2-icut1+1)
               else
                  index=(nrec-1+(jj-1)*ir)*(nx) 
     &                 + (isensor-1)*ir
               endif
               ir0=ir*(isensor-1)
               if (iopt(1).eq.1) then
c     For iopt(1)=1 (wavenumber inversion) 
                  do ii=icut1,icut2
                     resp(ii+index)=wavenoint(ii,nrec)
                  enddo
               elseif (iopt(1).eq.2) then
c---  multiply with sqrt (wavenumber) (moved from calint)
                  do ii=icut1,icut2
                     wavenoint(ii,nrec+ir0)=
     &                    wavenoint(ii,nrec+ir0)*facsqrt(ii)
                  enddo
                  if (nofft) then
c     if((flagpu.lt.0).and.(nrec.eq.1))
c     &              write(*,*)' Using trapezoidal integration'
                     if (tilt) then
                        xranges(1)=xrange_ref+dtilt/(ir-1)*(nrec-1)
                        if (flagpu.lt.0) 
     &                       write(prtfil,*) ' response for range',
     &                       xranges(1)
                     endif
                     call intgrn_pg(xranges,nrec+ir0,nflag)
c     write(*,*)' exited intgrn',nrec
                  else
                     CALL PHINT(NFLAG)
                     CALL TLOSS2
                  endif
c     
c     
c     
                  if (nofft) then
                     do i=1,nx
                        resp(i+index)=cff(i  ,1)
                     enddo
                     if (iopt(6).eq.1)
     .                    write(*,*)'resp(ran=1,dep,frq):'
     .                    ,nrec,jj,resp(1+index)
                  else

                     do i=1,nx
                        theta1=
     &                       ATAN2(AIMAG(cff(ixpoints(i)  ,1)),
     &                       REAL(cff(ixpoints(i)  ,1)))
                        theta2=
     &                       ATAN2(AIMAG(cff(ixpoints(i)+1,1)),
     &                       REAL(cff(ixpoints(i)+1,1)))
                        if (abs(theta1-theta2).gt.pi) then
                           if ((theta1.lt.0).and.(theta2.gt.0))then
                              theta1=theta1+2.*pi
                           elseif ((theta2.lt.0).and.(theta1.gt.0))then
                              theta2=theta2+2.*pi
                           endif
                        endif
                        theta=theta1*(1.-xweight(i))+
     &                       theta2*(xweight(i))
                        
                        resp(i+index)=(abs(cff(ixpoints(i) ,1))*
     &                       (1.-xweight(i))
     &                       +abs(cff(ixpoints(i)+1,1))*xweight(i))
     &                       *exp(cmplx(0.,theta))
                     enddo
                  endif
c     
c---- for the gradients
c     
                  do iparm=1,ngrad
                     do ii=icut1,icut2
                        wavenoint(ii,nrec)=
     &                       wavenograd(ii,nrec+ir0,iparm)*facsqrt(ii)
                     enddo
                     if (nofft) then
c     if((flagpu.lt.0).and.(nrec.eq.1))
c     &                 write(*,*)' Using trapezoidal integration'
                        if (tilt) then
                           xranges(1)=xrange_ref+dtilt/(ir-1)*(nrec-1)
                           if (flagpu.lt.0) 
     &                          write(prtfil,*)' response for range',
     &                          xranges(1)
                        endif
                        call intgrn_pg(xranges,nrec,nflag)
c     write(*,*)' exited intgrn',nrec
                     else
                        CALL PHINT(NFLAG)
                        CALL TLOSS2
                     endif
c     
c---- Transfer the gradients to grad
c     
                     if (nofft) then
c     write(*,*)'writing resp'
                        do i=1,nx
                           grad(i+index,iparm)=cff(i  ,1)
c     write(*,*)'grad', grad(i+index,iparm),i,index
                        enddo
                     else
                        do i=1,nx
                           grad(i+index,iparm)=
     &                          (cff(ixpoints(i),1))*(1.-xweight(i))
     &                          +(cff(ixpoints(i)+1,1))*xweight(i)
c     write(*,*)'grad',i+index,iparm,grad(i+index,iparm)
                        enddo
                     endif
                  enddo         ! gradinets
               endif
               do i=1,nx
                  write(80,'(i4,6e13.4)')jj,real(resp(i+index)), 
     1                 imag(resp(i+index)),
     1                 (real(grad(i+index,iparm)),
     1                 imag(grad(i+index,iparm)),
     2                 iparm=1,nparm)
               enddo            ! gradinets
 120        continue
            if (tilt) then
               xranges(1)=xrange_ref
            endif
 15      continue               ! frequency
     l           wascomp=1
c     if (1.eq.1) stop "stopped artificially"
         END  
c****************************************************************
      SUBROUTINE CALINTgrad(mindex)
      USE global
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnp.f'
      INCLUDE './oases/comnrd.f'
      INCLUDE 'comopt.h'
      INCLUDE 'comgrad.h'
      Complex wavenograd(np,nrd,mpar)
      common  /wavenogr/wavenograd
      COMPLEX FACSQ
C     *** OPEN SCRATCH FILE FOR KERNELS
      IF (IOUT(1).GT.0) THEN
         INDXCF=IR*3*(0)
      else IF (IOUT(2).GT.0) THEN
         INDXCF=IR*3*(1)
      else IF (IOUT(3).GT.0) THEN
         INDXCF=IR*3*(2)
      else
         write(*,*)'calint:unknown parameter'
      endif 

C     *** WAVENUMBER RAMP
      CALL VRAMP(WK0,DLWVNO,FAC,1,NWVNO)
      NGVALS=2*IR*3
C     *** WAVENUMBER LOOP

      facsq=1  
      DO 20 II=ICUT1,ICUT2         
         WVNO=CMPLX(FAC(II),OFFIMA)
         IF (ICDR.EQ.0) FACSQ=CSQRT(WVNO)
         facsqrt(ii)=facsq
         

         CALL INITS
         CALL BUILD  
         
         CALL SOLVE    

         CALL KERNEL(CFILE,IR)
c     write(*,*)'wavno resp',cfile(1)
c     *** tapering
         IF (II.LT.ICW1) THEN
            TFAC=(0.5*(1E0+COS((II-ICW1)*PI/(ICUT1-ICW1-1))))**4
            CALL VSMUL(CFILE(1),1,TFAC,CFILE(1),1,NGVALS)
         ELSE IF (II.GT.ICW2) THEN
            TFAC=(0.5*(1E0+COS((II-ICW2)*PI/(ICUT2-ICW2+1))))**4
            CALL VSMUL(CFILE(1),1,TFAC,CFILE(1),1,NGVALS)
         else
            tfac=1e0
         END IF
c     if (mod(ii,100).eq.0) write(*,*)ii,tap,icw1,icw2
         
         i0=0
         disp_to_vel= dsq*1000*1500*1e-6*1000 ! dsq*dens*Cref*um* (1000 ?)
c     write(*,*)'iout',iout(1) ,iout(2) ,iout(3), disp_to_vel,dsq
         IF (IOUT(1).GT.0) THEN
            INDXCF=IR*3*(0)
            do i=1,ir
               wavenoint(ii,i+i0)=cfile(indxcf+i)
            enddo
            i0=i0+ir
         endif 
         IF (IOUT(2).GT.0) THEN
            INDXCF=IR*(1)
            do i=1,ir
               wavenoint(ii,i+i0)=cfile(indxcf+i)* disp_to_vel
c     write(*,*)'cfile(indxcf+i)',cfile(indxcf+i)
            enddo
            i0=i0+ir
         endif 
         IF (IOUT(3).GT.0) THEN
            INDXCF=IR*(2)
            do i=1,ir
               wavenoint(ii,i+i0)=cfile(indxcf+i)* disp_to_vel
c     write(*,*)'cfile(indxcf+i)',cfile(indxcf+i)
            enddo
         endif 
c     
c---- for computations of gradient
c     
c--------move source potentials
         do jlay=1,numl
            ssr(jlay,1)=ss(jlay,1)
c     ssr(jlay,2)=ss(jlay,2)
            ssr(jlay,3)=ss(jlay,3)
c     ssr(jlay,4)=ss(jlay,4)
c     write(*,*)'ssr', ssr(jlay,1), ssr(jlay,3)
         enddo
         

         call gradoas

         do iparm=1,ngrad
c     write(*,*)' Iparm=',iparm
            if (par2phy_forw(iparm).eq.9) then
c     write(*,*)'range'
               do i=1,ir*nout
                  wavenograd(ii,i,iparm) = 
     &                 wavenoint(ii,i)*wvno*(0.,-1.)
               enddo

            else
               call kernelgrad(cfile,ir,iparm)
c     write(*,*)'after kernelgrad',cfile(1)
               if (irecderiv.eq.1)
     &              call kernelgradrec(cfile,ir,iparm)
c     write(*,*)'after kernelgradrec',cfile(1)

               IF ((II.LT.ICW1) .or.  (II.GT.ICW2)) THEN
                  CALL VSMUL(CFILE(1),1,TFAC,CFILE(1),1,NGVALS)
               END IF
               i0=0
               disp_to_vel= dsq*1000*1500*1e-6*1000 ! dsq*dens*Cref*um* (1000 ?)
c     write(*,*)'iout',iout(1) ,iout(2) ,iout(3), disp_to_vel,dsq
               IF (IOUT(1).GT.0) THEN
                  INDXCF=IR*3*(0)
                  do i=1,ir
                     wavenograd(ii,i+i0,iparm) = cfile(indxcf+i)
                  enddo
                  i0=i0+ir
               endif 
               IF (IOUT(2).GT.0) THEN
                  INDXCF=IR*(1)
                  do i=1,ir
                     wavenograd(ii,i+i0,iparm) = cfile(indxcf+i)
     .                    *disp_to_vel
c     write(*,*)'cfile(indxcf+i)',cfile(indxcf+i)
                  enddo
                  i0=i0+ir
               endif 
               IF (IOUT(3).GT.0) THEN
                  INDXCF=IR*(2)
                  do i=1,ir
                     wavenograd(ii,i+i0,iparm) = cfile(indxcf+i)
     .                    *disp_to_vel
c     write(*,*)'cfile(indxcf+i)',cfile(indxcf+i)
                  enddo
               endif 
            endif
         enddo


 20   CONTINUE     
C     
      end

c*******************************************************************

      SUBROUTINE BSOLVE                  
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnp.f'
c---  
      CALL CVIMOV(ALO,INDA,1,WORK1,2,NNA)
      CALL CVFILL(CNUL,WORK2,2,NNB)
      CALL CVMOVI(WORK1,2,INDB,1,WORK2,NNA)
      DO 5 IS=1,LSTOT
         CALL CVIMOV(Rgrad(1,1,IS),INDR,1,RHS(1+(IS-1)*NEQ),2,NEQ)
 5    CONTINUE
c     IF (DEBUG) THEN
c     DO 10 IS=1,LSTOT
c     WRITE(*,*) 'right hand side for source no',IS
c     DO 10 I=1,NEQ
c     WRITE(*,*)I,RHS(I+(IS-1)*NEQ)
c     10      CONTINUE
c     ENDIF
      CALL CBGEMR(WORK2,RHS,NEQ,NEQ,LSTOT,IBW,EPS)
      DO 30 IS=1,LSTOT
         CALL CVMOVI(RHS(1+(IS-1)*NEQ),2,INDS,1,SSgrad(1,1,IS),NEQ)
 30   CONTINUE
      RETURN        
      END           
c**************************************************************
      SUBROUTINE KERNELgrad(CKERN,NRCV,iparm)        
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnrd.f'
      COMPLEX CKERN(NRCV,3)
      COMPLEX ERALFA,ERBETA,ERALFM,ERBETM 
      COMPLEX CC,CC1,CC2,CC3,CC4,CWUFAC
      COMPLEX ZETA(NRD),AIRY(NRD),BIRY(NRD),
     &     AIRYD(NRD),BIRYD(NRD),ZTAM(NRD)

      CWUFAC=AI*DSQ*PCORR
C     *** RECEIVERS IN ISOVELOCITY FLUID LAYERS
c     vd$  permutation(NUMTR)
c     vd$  nodepchk
C     DEC$ INIT_DEP_FWD
      DO 10 J=1,NUMTR(1)
         INT=NRPNT(J,1)
         LL=LAY(INT)      
         ZZ=Z(INT)        
         IF (LL.NE.1) THEN
            ERALFM=CEXP(-ZZ*ALFA(LL))
         ELSE
            ERALFM=0E0
         END IF                 
         IF (LL.NE.NUML) THEN
            ERALFA=CEXP((ZZ-THICK(LL))*ALFA(LL))
         ELSE
            ERALFA=0E0
         END IF
         CC1=ssgrad(LL,1,iparm)*ERALFM
         CC3=ssgrad(LL,3,iparm)*ERALFA
         CKERN(INT,1)=-CON1(LL)*(CC1+CC3)
         CKERN(INT,2)=ALFA(LL)*(-CC1+CC3)
         CKERN(INT,3)=-WVNO*(CC1+CC3)
 10   CONTINUE
C     
C     AIRY SOLUTION IMPLEMENTED 840907
C     
C     *** RECEIVERS IN non-ISOVELOCITY FLUID LAYERS
c     vd$  permutation(NRPNT)
c     vd$  cncall
c     vd$  nodepchk
C     DEC$ INIT_DEP_FWD
c     vd$  select(concur)
      DO 20 J=1,NUMTR(2)
         INT=NRPNT(J,2)
         LL=LAY(INT)      
         ZZ=Z(INT)        
         ZETA(INT)=CCO(LL)*S2-ZZ*ACO(LL)-BCO(LL)
         CALL SCAIRY(ZETA(INT),AIRY(INT),BIRY(INT),
     &        AIRYD(INT),BIRYD(INT),ZTAM(INT))
         CC1=CEXP(AISC(LL)-ZTAM(INT))
         CC2=CEXP(ZTAM(INT)-BISC(LL))
         IF ((REAL(ACO(LL))).LT.0) THEN
            CC3=ssgrad(LL,1,iparm)*CC1*AIRY(INT)
     &           +ssgrad(LL,3,iparm)*CC2*BIRY(INT)
            CC4=ssgrad(LL,1,iparm)*CC1*AIRYD(INT)
     &           +ssgrad(LL,3,iparm)*CC2*BIRYD(INT)
         ELSE
            CC3=ssgrad(LL,1,iparm)*CC2*BIRY(INT)
     &           +ssgrad(LL,3,iparm)*CC1*AIRY(INT)
            CC4=ssgrad(LL,1,iparm)*CC2*BIRYD(INT)
     &           +ssgrad(LL,3,iparm)*CC1*AIRYD(INT)
         END IF
         CKERN(INT,1)=-CON1(LL)*CC3
         CKERN(INT,2)=-ACO(LL)*CC4
         CKERN(INT,3)=-WVNO*CC3          
 20   CONTINUE
C     *** RECEIVERS IN SOLID LAYERS
c     vd$  permutation(NRPNT)
c     vd$  nodepchk
C     DEC$ INIT_DEP_FWD
      DO 30 J=1,NUMTR(3)
         INT=NRPNT(J,3)
         LL=LAY(INT)      
         ZZ=Z(INT)        
         IF (LL.NE.1) THEN
            ERALFM=CEXP(-ZZ*ALFA(LL))                 
            ERBETM=CEXP(-ZZ*BETA(LL))            
         ELSE
            ERALFM=0E0
            ERBETM=0E0
         END IF
         IF (LL.NE.NUML) THEN
            ERALFA=CEXP((ZZ-THICK(LL))*ALFA(LL))
            ERBETA=CEXP((ZZ-THICK(LL))*BETA(LL))
         ELSE
            ERALFA=0E0
            ERBETA=0E0
         END IF
         CC1=ssgrad(LL,1,iparm)*ERALFM
         CC2=ssgrad(LL,2,iparm)*ERBETM
         CC3=ssgrad(LL,3,iparm)*ERALFA
         CC4=ssgrad(LL,4,iparm)*ERBETA
         CKERN(INT,1)=CON2(LL)*(CC1+CC3)+CON4(LL)*(CC4-CC2)
         CKERN(INT,2)=ALFA(LL)*(CC3-CC1)+WVNO*(CC2+CC4)
         CKERN(INT,3)=-WVNO*(CC1+CC3)+BETA(LL)*(CC2-CC4)
 30   CONTINUE

C     *** RECEIVERS IN TRANSVERSILY ISOTROPIC LAYERS
c     vd$  permutation(NRPNT)
c     vd$  nodepchk
C     DEC$ INIT_DEP_FWD
      DO 40 J=1,NUMTR(4)
         INT=NRPNT(J,4)
         LL=LAY(INT)      
         ZZ=Z(INT)        
         IF (LL.NE.1) THEN
            ERALFM=CEXP(-ZZ*ALFA(LL))                 
            ERBETM=CEXP(-ZZ*BETA(LL))            
         ELSE
            ERALFM=0E0
            ERBETM=0E0
         END IF
         IF (LL.NE.NUML) THEN
            ERALFA=CEXP((ZZ-THICK(LL))*ALFA(LL))
            ERBETA=CEXP((ZZ-THICK(LL))*BETA(LL))
         ELSE
            ERALFA=0E0
            ERBETA=0E0
         END IF
         CC1=ssgrad(LL,1,iparm)*ERALFM
         CC2=ssgrad(LL,2,iparm)*ERBETM
         CC3=ssgrad(LL,3,iparm)*ERALFA
         CC4=ssgrad(LL,4,iparm)*ERBETA
         IF (IOUT(1).GT.0) CKERN(INT,1)=-CON1(LL)*
     &        (ANSTD(13,LL)*CC1+ANSTD(14,LL)*CC2+
     &        ANSTD(17,LL)*CC3+ANSTD(18,LL)*CC4)
         IF (IOUT(2).GT.0) CKERN(INT,2)=-DSQ*
     &        (ANSTD(7,LL)*CC1+ANSTD(8,LL)*CC2+
     &        ANSTD(11,LL)*CC3+ANSTD(12,LL)*CC4)
         IF (IOUT(3).GT.0) CKERN(INT,3)=-DSQ*
     &        (ANSTD(5,LL)*CC1+ANSTD(6,LL)*CC2+
     &        ANSTD(9,LL)*CC3+ANSTD(10,LL)*CC4)
 40   continue
C     *** CONVERT DISPLACEMENTS TO VELOCITIES
      DO 500 J=1,IR
         CKERN(J,2)=CKERN(J,2)*CWUFAC
         CKERN(J,3)=-DSQ*PCORR*CKERN(J,3)
 500  CONTINUE
      RETURN           
      END              



c****************************************************************

      SUBROUTINE gradoas    
      USE global
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnp.f'
      INCLUDE './oases/comnrd.f'
      INCLUDE 'comopt.h'
      INCLUDE 'comgrad.h'
      complex  a11,a21,a31,a41,a12,a22,a32,a42,cc1
      complex  amat(4,2),cc, chelp, cdisp, cpres

c     write(*,*)' entering gradoases'
      do iparm=1,ngrad
         ilay=par2lay_forw(iparm)
         zdep= (v(ilay+1,1)-v(ilay,1))
         if (zdep.gt.0) then
            cc=exp(-alfa(ilay)* zdep)
         else
            cc=0
         endif
         if (par2phy_forw(iparm).eq.2 .or. 
     1        par2phy_forw(iparm).eq.4)then
c---  P-speed or P-att
            if (par2phy_forw(iparm).eq.2 ) then
c--   P-speed 
               chelp=ak2(ilay,1)/(alfa(ilay)*v(ilay,2))
            else
c--   P-att
               chelp=(ak(ilay,1)*DSQ)/(alfa(ilay)*v(ilay,2))
     1              *cmplx(0.,1.0/54.57505)
            endif
            cdisp = (1-alfa(ilay)*zdep)*cc*chelp
            cpres = -con1(ilay)*zdep*cc*chelp
            a11= chelp
            a12=-cdisp
            a32= chelp
            a31=-cdisp
            a22= cpres          ! rho* omega**2/rho_norm
            a41=-cpres          ! rho* omega**2/rho_norm
            a21=0
            a42=0
         elseif (par2phy_forw(iparm).eq.6) then
c     write(*,*)' Densities are opt'
c---  Density
            a11=0
            a12=0
            a31=0
            a32=0
            a21= con1(ilay)/v(ilay,6)
            a22= cc*con1(ilay)/v(ilay,6)
            a41=-cc*con1(ilay)/v(ilay,6)
            a42=-1.0*con1(ilay)/v(ilay,6)
         elseif (par2phy_forw(iparm).eq.7) then
c---  Thickness
            cdisp =alfa(ilay)*(-alfa(ilay)*cc)
            cpres =con1(ilay)*(-alfa(ilay)*cc)
            a12=-cdisp
            a31=-cdisp
            a22= cpres          ! rho* omega**2/rho_norm
            a41=-cpres          ! rho* omega**2/rho_norm
            a11=0
            a21=0
            a32=0
            a42=0
         else
            a11=0
            a12=0
            a21=0
            a22=0
            a31=0
            a32=0
            a41=0
            a42=0
         endif
c     
c---- The system matrix times the local derivatives
         if (ilay.ne.1) then        
            rgrad(ilay-1,1,iparm)= 
     1           -(a11*ss(ilay,1)+ a12*ss(ilay,3))*disnrm
            rgrad(ilay-1,3,iparm)= 
     1           -(a21*ss(ilay,1) + a22*ss(ilay,3))
         endif
         if (ilay.ne.numl) then
            rgrad(ilay  ,1,iparm)= 
     1           -(a31*ss(ilay,1)+ a32*ss(ilay,3))*disnrm
            rgrad(ilay  ,3,iparm)= 
     1           -(a41*ss(ilay,1) + a42*ss(ilay,3))         
         endif
      enddo
c     
c-----derivatives in source-layers
c     
      if (irecderiv.eq.1) then
         do iloop=1,Nparmlay(lays(1))
            iparm=lay2parm(lays(1),iloop)
c     write(*,*)'test:.iloop',iparm
C     *** SOURCE TERMS
            I=NSPNT(1,1)
            ilay=LAYS(I)
            LN=ilay
            CC1=1E0/ALFA(LN)
            if (par2phy_forw(iparm).eq.2 .or. 
     1           par2phy_forw(iparm).eq.4)then
c---  P-speed or P-att
               if (par2phy_forw(iparm).eq.2) then
c--   P-speed 
                  chelp=ak2(ilay,1)/(alfa(ilay)*v(ilay,2))
               else
c--   P-att
                  chelp=(ak(ilay,1)*DSQ)/(alfa(ilay)*v(ilay,2))
     1                 *cmplx(0.,1.0/54.57505)
               endif
               IF (LN.GT.1) THEN
                  CC=CPHFAC(I)*CEXP(-ZUS(I)*ALFA(LN))*chelp
                  rgrad(ilay-1,1,iparm)=(-zus(i))*CC *disnrm
     1                 + rgrad(ilay-1,1,iparm)
                  rgrad(ilay-1,3,iparm)=
     1                 -CON1(LN)*CC*(CC1*(-zus(i))-CC1**2)
     1                 + rgrad(ilay-1,3,iparm) 
               END IF
               IF (LN.LT.NUML) THEN
                  CC=CPHFAC(I)*CEXP(-ZLS(I)*ALFA(LN))*chelp
                  rgrad(ilay  ,1,iparm)=(-zls(i))*CC *disnrm
     1                 +rgrad(ilay  ,1,iparm)
                  rgrad(ilay  ,3,iparm)=
     1                 CON1(LN)*CC*(CC1*(-zls(i))-CC1**2)
     1                 +rgrad(ilay  ,3,iparm)
               END IF
            elseif (par2phy_forw(iparm).eq.6) then
c---  Density
               IF (LN.GT.1) THEN
                  CC=CPHFAC(I)*CEXP(-ZUS(I)*ALFA(LN))
                  rgrad(ilay-1,3,iparm)=-CON1(LN)*CC*CC1/v(ilay,6) !*(-1)
     1                 + rgrad(ilay-1,3,iparm) 
               END IF
               IF (LN.LT.NUML) THEN
                  CC=CPHFAC(I)*CEXP(-ZLS(I)*ALFA(LN))
                  rgrad(ilay  ,3,iparm)=CON1(LN)*CC*CC1/v(ilay,6) !*(-1)
     1                 +rgrad(ilay  ,3,iparm)
               END IF
            elseif (par2phy_forw(iparm).eq.7) then
c---  Thickness
c     note that thickness is just changed by moving the lower interface.
               IF (LN.LT.NUML) THEN
                  CC=CPHFAC(I)*CEXP(-ZLS(I)*ALFA(LN))
                  rgrad(ilay  ,1,iparm)=(-alfa(ilay))*CC *disnrm
     1                 +rgrad(ilay  ,1,iparm)
                  rgrad(ilay  ,3,iparm)=CON1(LN)*CC*cc1*(-alfa(ilay))
     1                 +rgrad(ilay  ,3,iparm)
               END IF
            elseif (par2phy_forw(iparm).eq.8) then
c---  Source depth
               IF (LN.GT.1) THEN
c     this correspond to the one for ilay, but with a change of sign; further 
c     the stress has also changed sign (because of the structure of the
c     original rhs.
                  CC=CPHFAC(I)*CEXP(-ZUS(I)*ALFA(LN))
                  rgrad(ilay-1,1,iparm)= (-alfa(ilay))*CC *disnrm
                  rgrad(ilay-1,3,iparm)=- CON1(LN)*CC*cc1*(-alfa(ilay))
c     write(*,*)' press at upper interface',rgrad(ilay-1,3,iparm)
               END IF
               IF (LN.LT.NUML) THEN
c     this correspond to the thickness-derivative* (-1)
                  CC=CPHFAC(I)*CEXP(-ZLS(I)*ALFA(LN))
                  rgrad(ilay  ,1,iparm)=- (-alfa(ilay))*CC *disnrm
                  rgrad(ilay  ,3,iparm)=
     1                 - CON1(LN)*CC*cc1*(-alfa(ilay))
c     write(*,*)' press at lower interface',rgrad(ilay,3,iparm)
               END IF
            endif
         enddo                  ! loop over iparm - in source layer
      endif                     ! derivatives in source layer
      
      lstot=ngrad
c     DEBUG=.TRUE. 
      call bsolve
      end
c*************************************************

      SUBROUTINE KERNELgradrec(CKERN,NRCV,iparm)        
      USE global
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnp.f'
      INCLUDE './oases/comnrd.f'
      INCLUDE 'comopt.h'
      INCLUDE 'comgrad.h'
      COMPLEX CKERN(NRCV,3)
      COMPLEX ERALFA,ERBETA,ERALFM,ERBETM 
      complex  amat(3,2),cc, chelp, cdisp, cpres,cc1

      CWUFAC=AI*DSQ*PCORR
      DO 10 J=1,NUMTR(1)
         INT=NRPNT(J,1)
         LL=LAY(INT)      
         if (par2lay_forw(iparm).eq.ll) then
c     write(*,*) 'derivatives in rec layer' 
c--------Start of receiver layer
            ilay=ll
            do i=1,2
               do k=1,3
                  amat(k,i)=(0.,0.)
               enddo
            enddo
            ZZ=Z(INT)        
            IF (LL.NE.1) THEN
               ERALFM=CEXP(-ZZ*ALFA(LL)) ! down (1)
            ELSE
               ERALFM=0E0
            END IF                 
            IF (LL.NE.NUML) THEN
               ERALFA=CEXP((ZZ-THICK(LL))*ALFA(LL)) !up (3)
               zdep=-(ZZ-THICK(LL))
            ELSE
               ERALFA=0E0
            END IF
            if (par2phy_forw(iparm).eq.2 .or. 
     &           par2phy_forw(iparm).eq.4)then
c---  P-speed or P-att
               if (par2phy_forw(iparm).eq.2 ) then
c--   P-speed 
                  chelp=ak2(ilay,1)/(alfa(ilay)*v(ilay,2))
               else
c--   P-att
                  chelp=(ak(ilay,1)*DSQ)/(alfa(ilay)*v(ilay,2))
     1                 *cmplx(0.,1.0/54.57505)
               endif
               cdisp = (1-alfa(ilay)*zdep)*chelp
               cpres = -con1(ilay)*zdep*chelp
               amat(1,2)= (1-alfa(ilay)*zdep)*chelp*eralfa
               amat(1,1)=-(1-alfa(ilay)*zz )*chelp*eralfm
               amat(3,1)=-(-con1(ilay) *zz  )*chelp*eralfm !rho*omega**2/rho_norm
               amat(3,2)=-(-con1(ilay) *zdep)*chelp*eralfa !rho*omega**2/rho_norm
            elseif (par2phy_forw(iparm).eq.6) then
c     write(*,*)' Densities are opt'
c---  Density
               amat(3,1)=-con1(ilay)/v(ilay,6)*eralfm 
               amat(3,2)=-con1(ilay)/v(ilay,6)*eralfa 
            elseif (par2phy_forw(iparm).eq.7) then
c---  Thickness
               cdisp =alfa(ilay)*(-alfa(ilay))
               cpres =con1(ilay)*(-alfa(ilay))
               amat(1,2)= cdisp*eralfa             
               amat(3,2)=-cpres*eralfa ! rho* omega**2/rho_norm
            endif
c     
c---- The system matrix times the local derivatives
c     
            CKERN(INT,2)=   CKERN(INT,2)
     1           +(amat(1,1)*ssr(ilay,1)+ amat(1,2)*ssr(ilay,3))*CWUFAC
            CKERN(INT,3)=   CKERN(INT,3)
     1           +(amat(2,1)*ssr(ilay,1)+ amat(2,2)*ssr(ilay,3))*
     1           (-DSQ*PCORR)
            CKERN(INT,1)=   CKERN(INT,1)
     1           +(amat(3,1)*ssr(ilay,1)+ amat(3,2)*ssr(ilay,3))         
         endif
c     Finished receiver layer.
c     
c     Now source layer.
c     
         if ((par2lay_forw(iparm).eq.lays(1)).and.
     &        (ll.eq.lays(1))) then
c     write(*,*) 'derivatives in source layer' 
            I=NSPNT(1,1)
            ilay=LAYS(I)
            LN=ilay
            CC1=1E0/ALFA(LN)
            zz=abs(zus(i)-z(int))
            if (par2phy_forw(iparm).eq.2 .or. 
     &           par2phy_forw(iparm).eq.4)then
c---  P-speed or P-att
               if (par2phy_forw(iparm).eq.2 ) then
c--   P-speed 
                  chelp=ak2(ilay,1)/(alfa(ilay)*v(ilay,2))
               else
c--   P-att
                  chelp=(ak(ilay,1)*DSQ)/(alfa(ilay)*v(ilay,2))
     1                 *cmplx(0.,1.0/54.57505)
               endif
               CC=CPHFAC(I)*CEXP(-Zz*ALFA(LN))*chelp
               CKERN(INT,2)=   CKERN(INT,2)
     1              +sign(1.0,z(int)-zus(i))*(-zz)*CC*CWUFAC
               CKERN(INT,1)=   CKERN(INT,1)
     1              -CON1(LN)*CC*(CC1*(-zz)-CC1**2)
               write(*,*)' not tested....'
            elseif (par2phy_forw(iparm).eq.6) then
c---  Density
               CC=CPHFAC(I)*CEXP(-Zz*ALFA(LN))
               CKERN(INT,1)=   CKERN(INT,1)
     1              -CON1(LN)*CC*CC1/v(ilay,6) !*(-1)
c---  Thickness--- the derivative is zero !
c     
c-----source depth
            elseif (par2phy_forw(iparm).eq.8) then
c     write(*,*) 'derivatives in source layer'
c     chelp=(-alfa(ln))
               chelp=(-alfa(ln))*SIGN(1.0,ZUS(I)-Z(int))
               CC=CPHFAC(I)*CEXP(-ZZ*ALFA(Ln))*chelp
c     write(*,*)'before',ckern(INT,1)
               CKERN(INT,1)=CKERN(INT,1)+CC*(-CON1(Ln))*CC1
               CKERN(INT,2)=CKERN(INT,2)
     1              -SIGN(1.0,Z(int)-ZUS(I))*CC*CWUFAC
               CKERN(INT,3)=CKERN(INT,3)+CC*(-WVNO)*CC1*(-DSQ*PCORR)
c     write(*,*)'after direct',ckern(INT,1)
            endif
         ENDIF                  ! source layer     

 10   CONTINUE
      end
      SUBROUTINE recderiv()        
      USE global
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnp.f'
      INCLUDE './oases/comnrd.f'
      INCLUDE 'comopt.h'
      INCLUDE 'comgrad.h'


      irecderiv=0               ! no derivatives in layers with rec or source

      do iparm=1,ngrad
         DO 10 J=1,NUMTR(1)
            INT=NRPNT(J,1)
            LL=LAY(INT)      
            if (par2phy_forw(iparm).eq.ll) then
               irecderiv=1   
            endif
 10      CONTINUE
         if (par2phy_forw(iparm).eq.8) then
            irecderiv=1   
         endif
         if (par2lay_forw(iparm).eq.lays(1)) then
            irecderiv=1   
         endif
      enddo
      write(*,*)' irecderiv=', irecderiv
      end
c*************************************
      subroutine setmodelrealjac(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comgrad.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnrd.f'     
      integer i,j,k
      real x(*)
      real  deltathick
      do k=1,nthet
        i=forw2opt(k)
        if (par2phy(i).eq.7) then
c***    this is the thickness of each layer
c          write(*,*)'fval', fval(model(i,iq),i)
          deltathick= 
     &        x(k)-(v(par2lay(i)+1,1)-v(par2lay(i),1))
c          write(*,*)'deltathick', deltathick 
          do j=par2lay(i)+1,nlay
            v(j,1)=v(j,1)+deltathick
          enddo
        elseif (par2phy(i).eq.8) then
c***    source depth...
          sdc(1)= x(k)
        elseif (par2phy(i).eq.9) then
c***    receiver range...
          xranges(1)= x(k)
        elseif (par2phy(i).eq.11) then
c***    this is the EOF
           aeof(par2lay(i)) = x(k)
        else
          v(par2lay(i),par2phy(i))=x(k)
        endif
      enddo
      end

