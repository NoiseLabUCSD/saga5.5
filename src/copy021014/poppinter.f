C*************
      SUBROUTINE setmodel(iq)
      USE global
      INTEGER iq
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'compopp.h'
      REAL fval                 ! function for computation of real values
      INTEGER i
      DO i=1,nparm
         IF (par2phy(i).EQ.1) THEN
c***  this is the water depth
            h=fval(model(i,iq),i)
c     WRITE(*,*)i,iq,model(i,iq),fval(model(i,iq),i),hw
         ELSEIF (par2phy(i).EQ.2) THEN
c***  this is the water velocity
c     WRITE(*,*)'fval', fval(model(i,iq),i)
            cp(par2lay(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.3) THEN
c***  this is the sediment velocity
            cbot(par2lay(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.4) THEN
c***  this is the attenuation
            alpbot(par2lay(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.5) THEN
c***  this is the thickness
            hbot(par2lay(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.6) THEN
c***  this is the sediment density
            rhobot(par2lay(i)) = fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.8) THEN
c***  source depth
            zs(1)=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.11) THEN
c***  this is the EOF
            aeof(par2lay(i)) = fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.15) THEN
c***  this is the receiver depth
            zr(1) = fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.16) THEN
c***  this is the velocity points in water
            zp(par2lay(i)) = fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.17) THEN
c***  this is the lambert scattering coefficient
            dmu(par2lay(i)) = fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.18) THEN
c***  this is the lambert power
            lampow = fval(model(i,iq),i)
         ELSE
            WRITE(*,*) ' setmodel: option not defined for parameter',i
         ENDIF
c     WRITE(prtfil,*)' setmodel:, i,fval',i,fval(model(i,iq),i)
      ENDDO
c     WRITE(*,*)'parameters'
c     WRITE(*,'(8f10.3)')(v(par2lay(i),par2phy(i)),i=1,nparm)      
c     DO i=1,nlay
c     WRITE(*,'(6f12.3)')(v(i,j),j=1,6)
c     ENDDO
      END
c*************
      SUBROUTINE setmodelbest(modelb)
      USE global
      INTEGER i,modelb(*)
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'compopp.h'
      REAL fval                 ! function for computation of real values
      DO i=1,nparm
         IF (par2phy(i).EQ.1) THEN
c***  this is the water depth
            h=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.2) THEN
c***  this is the water velocity
            cp(par2lay(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.3) THEN
c***  this is the sediment velocity
            cbot(par2lay(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.4) THEN
c***  this is the attenuation
            alpbot(par2lay(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.5) THEN
c***  this is the thickness
            hbot(par2lay(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.6) THEN
c***  this is the sediment density
            rhobot(par2lay(i)) = fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.8) THEN
c***  source depth
            zs(1)=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.11) THEN
c***  this is the EOF
            aeof(par2lay(i)) = fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.15) THEN
c***  this is the receiver depth
            zr(1) = fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.16) THEN
c***  this is the velocity points in water
            zp(par2lay(i)) = fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.17) THEN
c***  this is the lambert scattering coefficient
            dmu(par2lay(i)) = fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.18) THEN
c***  this is the lambert power
            lampow = fval(modelb(i),i)
         ELSE
            WRITE(*,*) ' setmodel: option not defined for parameter',i
         ENDIF
c     WRITE(prtfil,*)' setmodel:, i,fval',i,fval(model(i,iq),i)
      ENDDO
      END

      SUBROUTINE setmodelreal(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'compopp.h'
      INTEGER i
      REAL x(*)
      DO i=1,nparm
        IF (par2phy(i).EQ.1) THEN
c***    this is the water depth
           h=x(i)
c           WRITE(*,*)i,iq,model(i,iq),x(i),hw
        ELSEIF (par2phy(i).EQ.2) THEN
c***    this is the water velocity
c          WRITE(*,*)'fval', x(i)
           cp(par2lay(i))=x(i)
        ELSEIF (par2phy(i).EQ.3) THEN
c***    this is the sediment velocity
           cbot(par2lay(i))=x(i)
        ELSEIF (par2phy(i).EQ.4) THEN
c***    this is the attenuation
           alpbot(par2lay(i))=x(i)
        ELSEIF (par2phy(i).EQ.5) THEN
c***    this is the thickness
           hbot(par2lay(i))=x(i)
        ELSEIF (par2phy(i).EQ.6) THEN
c***    this is the sediment density
           rhobot(par2lay(i)) = x(i)
        ELSEIF (par2phy(i).EQ.8) THEN
c***    source depth
           zs(1)=x(i)
        ELSEIF (par2phy(i).EQ.11) THEN
c***    this is the EOF
           aeof(par2lay(i)) = x(i)
        ELSEIF (par2phy(i).EQ.15) THEN
c***    this is the receiver depth
           zr(1) = x(i)
        ELSEIF (par2phy(i).EQ.16) THEN
c***    this is the velocity points in water
           zp(par2lay(i)) = x(i)
        ELSEIF (par2phy(i).EQ.17) THEN
c***    this is the lambert scattering coefficient
           dmu(par2lay(i)) = x(i)
        ELSEIF (par2phy(i).EQ.18) THEN
c***    this is the lambert power
           lampow = x(i)
        ELSE
          WRITE(*,*) ' setmodel: option not defined for parameter',i
        ENDIF
c        WRITE(prtfil,*)' setmodel:, i,fval',i,x(i)
      ENDDO
      END
c********************
      SUBROUTINE getmodelreal(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'compopp.h'
      INTEGER i
      REAL x(*)
      DO i=1,nparm
         IF (par2phy(i).EQ.1) THEN
c***  this is the water depth
            x(i)=h
         ELSEIF (par2phy(i).EQ.2) THEN
c***  this is the water velocity
            x(i)=cp(par2lay(i))
         ELSEIF (par2phy(i).EQ.3) THEN
c***  this is the sediment velocity
            x(i)=cbot(par2lay(i))
         ELSEIF (par2phy(i).EQ.4) THEN
c***  this is the attenuation
            x(i)=alpbot(par2lay(i))
         ELSEIF (par2phy(i).EQ.5) THEN
c***  this is the thickness
            x(i)=hbot(par2lay(i))
         ELSEIF (par2phy(i).EQ.6) THEN
c***  this is the sediment density
            x(i)=rhobot(par2lay(i))
         ELSEIF (par2phy(i).EQ.8) THEN
c***  source depth
            x(i)=zs(1)
         ELSEIF (par2phy(i).EQ.11) THEN
c***  this is the EOF
            x(i)=aeof(par2lay(i)) 
         ELSEIF (par2phy(i).EQ.15) THEN
c***  this is the receiver depth
            x(i)=zr(1) 
         ELSEIF (par2phy(i).EQ.16) THEN
c***  this is the velocity points in water
            x(i)=zp(par2lay(i))
         ELSEIF (par2phy(i).EQ.17) THEN
c***  this is the lambert scattering coefficient
            x(i)=dmu(par2lay(i)) 
         ELSEIF (par2phy(i).EQ.18) THEN
c***  this is the lambert power
            x(i)=lampow
         ELSE
            WRITE(*,*) ' setmodel: option not defined for parameter',i
         ENDIF
c     WRITE(prtfil,*)' setmodel:, i,fval',i,fval(model(i,iq),i)
      ENDDO
      END
c*********************
      SUBROUTINE setmodelx(ixlay,ixpar,theta)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'compopp.h'     
      INTEGER ixlay,ixpar
      REAL theta
      IF (ixpar.EQ.1) THEN
c***  this is the water depth
         h=theta
c     WRITE(*,*)i,iq,model(i,iq),theta,hw
      ELSEIF (ixpar.EQ.2) THEN
c***  this is the water velocity
c     WRITE(*,*)'fval', theta
         cp(ixlay)=theta
      ELSEIF (ixpar.EQ.3) THEN
c***  this is the sediment velocity
         cbot(ixlay)=theta
      ELSEIF (ixpar.EQ.4) THEN
c***  this is the attenuation
         alpbot(ixlay)=theta
      ELSEIF (ixpar.EQ.5) THEN
c***  this is the thickness
         hbot(ixlay)=theta
      ELSEIF (ixpar.EQ.6) THEN
c***  this is the sediment density
         rhobot(ixlay) = theta
      ELSEIF (ixpar.EQ.8) THEN
c***  source depth
         zs(1)=theta
      ELSEIF (ixpar.EQ.11) THEN
c***  this is the EOF
         aeof(ixlay) = theta
      ELSEIF (ixpar.EQ.15) THEN
c***  this is the receiver depth
         zr(1) = theta
      ELSEIF (ixpar.EQ.16) THEN
c***  this is the velocity points in water
         zp(ixlay) = theta
      ELSEIF (ixpar.EQ.17) THEN
c***  this is the lambert scattering coefficient
         dmu(ixlay) = theta
      ELSEIF (ixpar.EQ.18) THEN
c***  this is the lambert power
         lampow = theta
      ELSE
         WRITE(*,*) ' setmodel: option not defined for parameter'
      ENDIF
      END

c*************
      SUBROUTINE setmodeleof(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'compopp.h'
      INTEGER i
      REAL x(*)
      DO i=1,neofvar
         IF (par2phy_eof(i).EQ.1) THEN
c***  this is the water depth
            h=x(i)
c     WRITE(*,*)i,iq,model(i,iq),x(i),hw
         ELSEIF (par2phy_eof(i).EQ.2) THEN
c***  this is the water velocity
c     WRITE(*,*)'fval', x(i)
            cp(par2lay_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.3) THEN
c***  this is the sediment velocity
            cbot(par2lay_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.4) THEN
c***  this is the attenuation
            alpbot(par2lay_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.5) THEN
c***  this is the thickness
            hbot(par2lay_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.6) THEN
c***  this is the sediment density
            rhobot(par2lay_eof(i)) = x(i)
         ELSEIF (par2phy_eof(i).EQ.8) THEN
c***  source depth
            zs(1)=x(i)
         ELSEIF (par2phy_eof(i).EQ.11) THEN
c***  this is the EOF
            aeof(par2lay_eof(i)) = x(i)
         ELSEIF (par2phy_eof(i).EQ.15) THEN
c***  this is the receiver depth
            zr(1) = x(i)
         ELSEIF (par2phy_eof(i).EQ.16) THEN
c***  this is the velocity points in water
            zp(par2lay_eof(i)) = x(i)
         ELSEIF (par2phy_eof(i).EQ.17) THEN
c***  this is the lambert scattering coefficient
            dmu(par2lay_eof(i)) = x(i)
         ELSEIF (par2phy_eof(i).EQ.18) THEN
c***  this is the lambert power
            lampow = x(i)
         ELSE
            WRITE(*,*) ' setmodel: option not defined for parameter',i
         ENDIF
c     WRITE(prtfil,*)' setmodel:, i,fval',i,x(i)
      ENDDO
      END












