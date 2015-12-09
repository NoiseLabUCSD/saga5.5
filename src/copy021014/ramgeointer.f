C*************
      SUBROUTINE setmodel(iq)
      USE global
      INTEGER iq
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './ramgeo/ram.h'
      INCLUDE 'comramgeo.h'
      REAL fval                 ! function for computation of real values
      INTEGER i
      DO i=1,nparm
         IF (par2phy(i).EQ.1) THEN
c***  this is the ocean sound speed
            geopar(par2lay(i),1,par3(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.2) THEN
c***  this is the sediment sound speed
            geopar(par2lay(i),2,par3(i))=fval(model(i,iq),i)
            WRITE(*,*)'geopar', geopar(par2lay(i),2,par3(i))
            WRITE(*,*)par2lay(i),par3(i)
         ELSEIF (par2phy(i).EQ.3) THEN
c***  this is the bottom density
            geopar(par2lay(i),3,par3(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.4) THEN
c***  this is bottom attenuation
            geopar(par2lay(i),4,par3(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.5) THEN
c***  this is the ocean sound speed depth
            geodep(par2lay(i),1,par3(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.6) THEN
c***  this is the sediment sound speed depth
            geodep(par2lay(i),2,par3(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.7) THEN
c***  this is the bottom density depth
            geodep(par2lay(i),3,par3(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.10) THEN
            geodep(par2lay(i),4,par3(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.8) THEN
c***  this is bottom attenuation depth
            zs=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.9) THEN
c***  this is max range
            rmin=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.12) THEN
c***  this is bathymetry
            zb(par2lay(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.12) THEN
c***  this is bathymetry
            zb(par2lay(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.13) THEN
c***  this is bathymetry range
            rb(par2lay(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.11) THEN
c***  this is the EOF
            aeof(par2lay(i)) = fval(model(i,iq),i)
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
      INCLUDE './ramgeo/ram.h'
      INCLUDE 'comramgeo.h'
      REAL fval                 ! function for computation of real values
      DO i=1,nparm
         IF (par2phy(i).EQ.1) THEN
c***  this is the ocean sound speed
            geopar(par2lay(i),1,par3(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.2) THEN
c***  this is the sediment sound speed
            geopar(par2lay(i),2,par3(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.3) THEN
c***  this is the bottom density
            geopar(par2lay(i),3,par3(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.4) THEN
c***  this is bottom attenuation
            geopar(par2lay(i),4,par3(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.5) THEN
c***  this is the ocean sound speed depth
            geodep(par2lay(i),1,par3(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.6) THEN
c***  this is the sediment sound speed depth
            geodep(par2lay(i),2,par3(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.7) THEN
c***  this is the bottom density depth
            geodep(par2lay(i),3,par3(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.10) THEN
            geodep(par2lay(i),4,par3(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.8) THEN
c***  this is bottom attenuation depth
            zs=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.9) THEN
c***  this is max range
            rmin=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.12) THEN
c***  this is bathymetry
            zb(par2lay(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.13) THEN
c***  this is bathymetry range
            rb(par2lay(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.11) THEN
c***  this is the EOF
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
      INCLUDE './ramgeo/ram.h'
      INCLUDE 'comramgeo.h'
      INTEGER i
      REAL x(*)
      DO i=1,nparm
         IF (par2phy(i).EQ.1) THEN
c***  this is the ocean sound speed
            geopar(par2lay(i),1,par3(i))=x(i)
         ELSEIF (par2phy(i).EQ.2) THEN
c***  this is the sediment sound speed
            geopar(par2lay(i),2,par3(i))=x(i)
         ELSEIF (par2phy(i).EQ.3) THEN
c***  this is the bottom density
            geopar(par2lay(i),3,par3(i))=x(i)
         ELSEIF (par2phy(i).EQ.4) THEN
c***  this is bottom attenuation
            geopar(par2lay(i),4,par3(i))=x(i)
         ELSEIF (par2phy(i).EQ.5) THEN
c***  this is the ocean sound speed depth
            geodep(par2lay(i),1,par3(i))=x(i)
         ELSEIF (par2phy(i).EQ.6) THEN
c***  this is the sediment sound speed depth
            geodep(par2lay(i),2,par3(i))=x(i)
         ELSEIF (par2phy(i).EQ.7) THEN
c***  this is the bottom density depth
            geodep(par2lay(i),3,par3(i))=x(i)
         ELSEIF (par2phy(i).EQ.10) THEN
            geodep(par2lay(i),4,par3(i))=x(i)
         ELSEIF (par2phy(i).EQ.8) THEN
c***  this is bottom attenuation depth
            zs=x(i)
         ELSEIF (par2phy(i).EQ.9) THEN
c***  this is max range
            rmin=x(i)
         ELSEIF (par2phy(i).EQ.12) THEN
c***  this is bathymetry
            zb(par2lay(i))=x(i)
         ELSEIF (par2phy(i).EQ.13) THEN
c***  this is bathymetry range
            rb(par2lay(i))=x(i)
         ELSEIF (par2phy(i).EQ.11) THEN
c***  this is the EOF
            aeof(par2lay(i)) = x(i)
         ELSE
            WRITE(*,*)'setmodelreal: opt. not defined for parameter',i
         ENDIF
      ENDDO
      END
c********************
      SUBROUTINE getmodelreal(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './ramgeo/ram.h'
      INCLUDE 'comramgeo.h'
      INTEGER i
      REAL x(*)
      DO i=1,nparm

         IF (par2phy(i).EQ.1) THEN
c***  this is the ocean sound speed
            x(i)=geopar(par2lay(i),1,par3(i))
         ELSEIF (par2phy(i).EQ.2) THEN
c***  this is the sediment sound speed
            x(i)=geopar(par2lay(i),2,par3(i))
         ELSEIF (par2phy(i).EQ.3) THEN
c***  this is the bottom density
            x(i)=geopar(par2lay(i),3,par3(i))
         ELSEIF (par2phy(i).EQ.4) THEN
c***  this is bottom attenuation
            x(i)=geopar(par2lay(i),4,par3(i))
         ELSEIF (par2phy(i).EQ.5) THEN
c***  this is the ocean sound speed depth
            x(i)= geodep(par2lay(i),1,par3(i))
         ELSEIF (par2phy(i).EQ.6) THEN
c***  this is the sediment sound speed depth
            x(i)=geodep(par2lay(i),2,par3(i))
         ELSEIF (par2phy(i).EQ.7) THEN
c***  this is the bottom density depth
            x(i)=geodep(par2lay(i),3,par3(i))
         ELSEIF (par2phy(i).EQ.10) THEN
c***  this is the bottom attenuation depth
            x(i)=geodep(par2lay(i),4,par3(i))
         ELSEIF (par2phy(i).EQ.8) THEN
c***  this is depth
            x(i)=zs
         ELSEIF (par2phy(i).EQ.9) THEN
c***  this is max range
            x(i)=rmin
         ELSEIF (par2phy(i).EQ.12) THEN
c***  this is bathymetry
            x(i)=zb(par2lay(i))
         ELSEIF (par2phy(i).EQ.13) THEN
c***  this is bathymetry range
            x(i)=rb(par2lay(i))
         ELSEIF (par2phy(i).EQ.11) THEN
c***  this is the EOF
            x(i)=aeof(par2lay(i)) 
         ELSE
            WRITE(*,*)'Getmodelreal: opt. not defined for parameter',i
         ENDIF
      ENDDO
      END
c*********************
      SUBROUTINE setmodelx(ixlay,ixpar,iz,theta)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './ramgeo/ram.h'
      INCLUDE 'comramgeo.h'     
      INTEGER ixlay,ixpar,iz
      REAL theta
      IF (ixpar.EQ.1) THEN
c***  this is the ocean sound speed
         geopar(ixlay,1,iz)=theta
      ELSEIF (ixpar.EQ.2) THEN
c***  this is the sediment sound speed
         geopar(ixlay,2,iz)=theta
      ELSEIF (ixpar.EQ.3) THEN
c***  this is the bottom density
         geopar(ixlay,3,iz)=theta
      ELSEIF (ixpar.EQ.4) THEN
c***  this is bottom attenuation
         geopar(ixlay,4,iz)=theta
      ELSEIF (ixpar.EQ.5) THEN
c***  this is the ocean sound speed depth
         geodep(ixlay,1,iz)=theta
      ELSEIF (ixpar.EQ.6) THEN
c***  this is the sediment sound speed depth
         geodep(ixlay,2,iz)=theta
      ELSEIF (ixpar.EQ.7) THEN
c***  this is the bottom density depth
         geodep(ixlay,3,iz)=theta
      ELSEIF (ixpar.EQ.10) THEN
         geodep(ixlay,4,iz)=theta
      ELSEIF (ixpar.EQ.8) THEN
c***  this is bottom attenuation depth
         zs=theta
      ELSEIF (ixpar.EQ.9) THEN
c***  this is max range
         rmin=theta
      ELSEIF (ixpar.EQ.12) THEN
c***  this is bathymetry
         zb(ixlay)=theta
      ELSEIF (ixpar.EQ.13) THEN
c***  this is bathymetry range
         rb(ixlay)=theta
      ELSEIF (ixpar.EQ.11) THEN
c***  this is the EOF
         aeof(ixlay) = theta
      ELSEIF (ixpar.EQ.11) THEN
c***  this is the EOF
         aeof(ixlay) = theta
      ELSE
         WRITE(*,*)'setmodelx: opt. not defined for parameter'
      ENDIF
      END

c*************
      SUBROUTINE setmodeleof(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './ramgeo/ram.h'
      INCLUDE 'comramgeo.h'
      INTEGER i
      REAL x(*)
      DO i=1,neofvar
         IF (par2phy_eof(i).EQ.1) THEN
c***  this is the ocean sound speed
            geopar(par2lay_eof(i),1,par3_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.2) THEN
c***  this is the sediment sound speed
            geopar(par2lay_eof(i),2,par3_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.3) THEN
c***  this is the bottom density
            geopar(par2lay_eof(i),3,par3_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.4) THEN
c***  this is bottom attenuation
            geopar(par2lay_eof(i),4,par3_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.5) THEN
c***  this is the ocean sound speed depth
            geodep(par2lay_eof(i),1,par3_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.6) THEN
c***  this is the sediment sound speed depth
            geodep(par2lay_eof(i),2,par3_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.7) THEN
c***  this is the bottom density depth
            geodep(par2lay_eof(i),3,par3_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.10) THEN
            geodep(par2lay_eof(i),4,par3_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.8) THEN
c***  this is bottom attenuation depth
            zs=x(i)
         ELSEIF (par2phy_eof(i).EQ.9) THEN
c***  this is max range
            rmin=x(i)
         ELSEIF (par2phy_eof(i).EQ.12) THEN
c***  this is bathymetry
            zb(par2lay_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.13) THEN
c***  this is bathymetry range
            rb(par2lay_eof(i))=x(i)
         ELSE
            WRITE(*,*)'setmodeleof: opt. not defined for parameter',i
         ENDIF
      ENDDO
      END
