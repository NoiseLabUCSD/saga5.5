C*************
      SUBROUTINE setmodel(iq)
      USE global
      INTEGER iq
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'
      REAL fval                 ! function for computation of real values
      INTEGER i
      DO i=1,nparm
         IF (par2phy(i).EQ.1) THEN
c***  this is the refractivity
            rf.refmsl(par2lay(i),par3(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.2) THEN
c***  this is the height corresponding to a refractivity point.
c     WRITE(*,*)'fval', fval(model(i,iq),i)
            rf.hmsl(par2lay(i),par3(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.3) THEN
c***  this is the terrain X-coordinates
            tr.terx(par2lay(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.4) THEN
c***  this is the terrain Y-coordinates
            tr.tery(par2lay(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.5) THEN
c***  this is the source-receiver distance
            vnp.rmax =fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.6) THEN
c***  this is the antenna height
            sv.antht = fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.7) THEN
c***  this is the capping
            cap_coef(par2lay(i)) = fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.9) THEN
c***  this is the source-receiver distance
            vnp.rmax =fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.12) THEN
c***  this is the baseheight
            rp.base(par2lay(i)) =fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.13) THEN
c***  this is the thickness
            rp.thick(par2lay(i)) =fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.14) THEN
c***  this is the offset
            rp.offset(par2lay(i)) =fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.15) THEN
c***  this is the m-defecit
            rp.mdef(par2lay(i)) =fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.16) THEN
c***  this is the noise-floor
            znoise =fval(model(i,iq),i)

         ELSEIF (par2phy(i).EQ.17) THEN
c***  this is the coeficients
            rp.coef(par2lay(i), par3(i)) =fval(model(i,iq),i)

         ELSEIF (par2phy(i).EQ.18) THEN
c***  clutter cross section
            ccs(par2lay(i)) =fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.19) THEN
c***  this is the rd factors
            rp.factor(par2lay(i), par3(i)) =fval(model(i,iq),i)
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
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'
      REAL fval                 ! function for computation of real values
      DO i=1,nparm
         IF (par2phy(i).EQ.1) THEN
c***  this is the refractivity
            rf.refmsl(par2lay(i),par3(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.2) THEN
c***  this is the height corresponding to a refractivity point.
c     WRITE(*,*)'fval', fval(modelb(i),i)
            rf.hmsl(par2lay(i),par3(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.3) THEN
c***  this is the terrain X-coordinates
            tr.terx(par2lay(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.4) THEN
c***  this is the terrain Y-coordinates
            tr.tery(par2lay(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.5) THEN
c***  this is the source-receiver distance
            vnp.rmax =fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.6) THEN
c***  this is the antenna height
            sv.antht = fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.7) THEN
c***  this is the capping coeff
            cap_coef(par2lay(i)) = fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.9) THEN
c***  this is the source-receiver distance
            vnp.rmax =fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.11) THEN
c***  this is the EOF
            aeof(par2lay(i)) = fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.12) THEN
c***  this is the baseheight
            rp.base(par2lay(i)) =fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.13) THEN
c***  this is the thickness
            rp.thick(par2lay(i)) =fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.14) THEN
c***  this is the offset
            rp.offset(par2lay(i)) =fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.15) THEN
c***  this is the m-defecit
            rp.mdef(par2lay(i)) =fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.16) THEN
c***  this is the noise-floor
            znoise =fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.17) THEN
c***  this is the coeficients
            rp.coef(par2lay(i), par3(i)) =fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.18) THEN
c***  clutter cross section
            ccs(par2lay(i)) =fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.19) THEN
c***  this is the rd factors
            rp.factor(par2lay(i), par3(i)) =fval(modelb(i),i)
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
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'
      INTEGER i
      REAL x(*)
      DO i=1,nparm
         IF (par2phy(i).EQ.1) THEN
c***  this is the refractivity
            rf.refmsl(par2lay(i),par3(i))=x(i)
         ELSEIF (par2phy(i).EQ.2) THEN
c***  this is the height corresponding to a refractivity point.
c     WRITE(*,*)'fval', fval(modelb(i),i)
            rf.hmsl(par2lay(i),par3(i))=x(i)
         ELSEIF (par2phy(i).EQ.3) THEN
c***  this is the terrain X-coordinates
            tr.terx(par2lay(i))=x(i)
         ELSEIF (par2phy(i).EQ.4) THEN
c***  this is the terrain Y-coordinates
            tr.tery(par2lay(i))=x(i)
         ELSEIF (par2phy(i).EQ.5) THEN
c***  this is the source-receiver distance
            vnp.rmax =x(i)
         ELSEIF (par2phy(i).EQ.6) THEN
c***  this is the antenna height
            sv.antht = x(i)
         ELSEIF (par2phy(i).EQ.7) THEN
c***  this is the capping coeff
            cap_coef(par2lay(i)) = x(i)
         ELSEIF (par2phy(i).EQ.9) THEN
c***  this is the source-receiver distance
            vnp.rmax =x(i)
         ELSEIF (par2phy(i).EQ.11) THEN
c***  this is the EOF
            aeof(par2lay(i)) = x(i)
         ELSEIF (par2phy(i).EQ.12) THEN
c***  this is the baseheight
            rp.base (par2lay(i))=x(i)
         ELSEIF (par2phy(i).EQ.13) THEN
c***  this is the thickness
            rp.thick(par2lay(i)) =x(i)
         ELSEIF (par2phy(i).EQ.14) THEN
c***  this is the offset
            rp.offset(par2lay(i)) =x(i)
         ELSEIF (par2phy(i).EQ.15) THEN
c***  this is the m-defecit
            rp.mdef(par2lay(i)) =x(i)
         ELSEIF (par2phy(i).EQ.16) THEN
c***  this is the noise-floor
            znoise =x(i)
         ELSEIF (par2phy(i).EQ.17) THEN
c***  this is the coeficients
            rp.coef(par2lay(i), par3(i)) =x(i)

         ELSEIF (par2phy(i).EQ.18) THEN
c***  clutter cross section
            ccs(par2lay(i)) =x(i)
         ELSEIF (par2phy(i).EQ.19) THEN
c***  this is the rd factors
            rp.factor(par2lay(i), par3(i)) =x(i)
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
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'
      INTEGER i
      REAL x(*)
      DO i=1,nparm
         IF (par2phy(i).EQ.1) THEN
c***  this is the refractivity
            x(i)=rf.refmsl(par2lay(i),par3(i))
         ELSEIF (par2phy(i).EQ.2) THEN
c***  this is the height corresponding to a refractivity point.
c     WRITE(*,*)'fval', fval(modelb(i),i)
            x(i)= rf.hmsl(par2lay(i),par3(i))
         ELSEIF (par2phy(i).EQ.3) THEN
c***  this is the terrain X-coordinates
            x(i)=tr.terx(par2lay(i))
         ELSEIF (par2phy(i).EQ.4) THEN
c***  this is the terrain Y-coordinates
            x(i)=tr.tery(par2lay(i))
         ELSEIF (par2phy(i).EQ.5) THEN
c***  this is the source-receiver distance
            x(i)=vnp.rmax 
         ELSEIF (par2phy(i).EQ.6) THEN
c***  this is the antenna height
            x(i)=sv.antht 
         ELSEIF (par2phy(i).EQ.7) THEN
c***  this is the capping coeff
            x(i)= cap_coef(par2lay(i))
         ELSEIF (par2phy(i).EQ.9) THEN
c***  this is the source-receiver distance
            x(i)=vnp.rmax 
         ELSEIF (par2phy(i).EQ.11) THEN
c***  this is the EOF
            x(i)=aeof(par2lay(i)) 
         ELSEIF (par2phy(i).EQ.12) THEN
c***  this is the baseheight
            x(i)=rp.base (par2lay(i))
         ELSEIF (par2phy(i).EQ.13) THEN
c***  this is the thickness
            x(i)=rp.thick(par2lay(i))
         ELSEIF (par2phy(i).EQ.14) THEN
c***  this is the offset
            x(i)=rp.offset(par2lay(i))
         ELSEIF (par2phy(i).EQ.15) THEN
c***  this is the m-defecit
            x(i)=rp.mdef(par2lay(i))
         ELSEIF (par2phy(i).EQ.16) THEN
c***  this is the noise-floor
            x(i)=znoise 
         ELSEIF (par2phy(i).EQ.17) THEN
c***  this is the coeficients
            x(i)=rp.coef(par2lay(i), par3(i)) 
         ELSEIF (par2phy(i).EQ.18) THEN
c***  clutter cross section
            x(i)= ccs(par2lay(i)) 
         ELSEIF (par2phy(i).EQ.19) THEN
c***  this is the rd factors
            x(i)=rp.factor(par2lay(i), par3(i))
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
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'     
      INTEGER ixlay,ixpar,iz
      REAL theta
      IF (ixpar.EQ.1) THEN
c***  this is the refractivity
         rf.refmsl(ixlay,iz)=theta
      ELSEIF (ixpar.EQ.2) THEN
c***  this is the height corresponding to a refractivity point.
c     WRITE(*,*)'fval', fval(modelb(i),i)
         rf.hmsl(ixlay,iz)=theta
      ELSEIF (ixpar.EQ.3) THEN
c***  this is the terrain X-coordinates
         tr.terx(ixlay)=theta
      ELSEIF (ixpar.EQ.4) THEN
c***  this is the terrain Y-coordinates
         tr.tery(ixlay)=theta
      ELSEIF (ixpar.EQ.5) THEN
c***  this is the source-receiver distance
         vnp.rmax =theta
      ELSEIF (ixpar.EQ.6) THEN
c***  this is the antenna height
         sv.antht=theta 
      ELSEIF (ixpar.EQ.7) THEN
c***  this is the capping coeff
         cap_coef(ixlay) =  theta
      ELSEIF (ixpar.EQ.9) THEN
c***  this is the source-receiver distance
         vnp.rmax =theta
      ELSEIF (ixpar.EQ.11) THEN
c***  this is the EOF
         aeof(ixlay) = theta
      ELSEIF (ixpar.EQ.12) THEN
c***  this is the baseheight
         rp.base(ixlay) =theta
      ELSEIF (ixpar .EQ.13) THEN
c***  this is the thickness
         rp.thick(ixlay) =theta
      ELSEIF (ixpar .EQ.14) THEN
c***  this is the offset
         rp.offset(ixlay) =theta
      ELSEIF (ixpar .EQ.15) THEN
c***  this is the m-defecit
         rp.mdef(ixlay) =theta
      ELSEIF (ixpar .EQ.16) THEN
c***  this is the noise-floor
         znoise =theta
      ELSEIF (ixpar.EQ.17) THEN
c***  this is the coeficients
         rp.coef(ixlay, iz) =theta
      ELSEIF (ixpar.EQ.18) THEN
c***  clutter cross section
         ccs(ixlay) = theta 
      ELSEIF (ixpar.EQ.19) THEN
c***  this is the rd factors
         rp.factor(ixlay, iz) =theta
      ELSE
         WRITE(*,*)'setmodelx: opt. not defined for parameter'
      ENDIF
      END

c*************
      SUBROUTINE setmodeleof(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'
      INTEGER i
      REAL x(*)
      DO i=1,neofvar
         IF (par2phy_eof(i).EQ.1) THEN
c***  this is the refractivity
            rf.refmsl(par2lay_eof(i),par3_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.2) THEN
c***  this is the height corresponding to a refractivity point.
            rf.hmsl(par2lay_eof(i),par3_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.3) THEN
c***  this is the terrain X-coordinates
            tr.terx(par2lay_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.4) THEN
c***  this is the terrain Y-coordinates
            tr.tery(par2lay_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.5) THEN
c***  this is the source-receiver distance
            vnp.rmax =x(i)
         ELSEIF (par2phy_eof(i).EQ.6) THEN
c***  this is the antenna height
            sv.antht = x(i)
         ELSEIF (par2phy_eof(i).EQ.7) THEN
c***  this is the capping
            cap_coef(par2lay_eof(i)) = x(i)
         ELSEIF (par2phy_eof(i).EQ.9) THEN
c***  this is the source-receiver distance
            vnp.rmax =x(i)
         ELSEIF (par2phy_eof(i).EQ.12) THEN
c***  this is the baseheight
            rp.base(par2lay_eof(i)) = x(i)
         ELSEIF (par2phy_eof(i).EQ.13) THEN
c***  this is the thickness
            rp.thick(par2lay_eof(i)) = x(i)
         ELSEIF (par2phy_eof(i).EQ.14) THEN
c***  this is the offset
            rp.offset(par2lay_eof(i)) = x(i)
         ELSEIF (par2phy_eof(i).EQ.15) THEN
c***  this is the m-defecit
            rp.mdef(par2lay_eof(i)) = x(i)
         ELSEIF (par2phy_eof(i).EQ.16) THEN
c***  this is the noise-floor
            znoise = x(i)
         ELSEIF (par2phy_eof(i).EQ.17) THEN
c***  this is the coeficients
            rp.coef(par2lay_eof(i), par3_eof(i)) =x(i)
         ELSEIF (par2phy_eof(i).EQ.18) THEN
c***  clutter cross section
            ccs(par2lay_eof(i)) = x(i)
         ELSEIF (par2phy_eof(i).EQ.19) THEN
c***  this is the rd factors
            rp.factor(par2lay_eof(i), par3_eof(i)) =x(i)
         ELSE
            WRITE(*,*)'setmodeleof: opt. not defined for parameter',i
         ENDIF
      ENDDO
      END
