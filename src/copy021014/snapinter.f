      
c------------------------------------------------------------------------
c     Map environment variables to SNAP
c------------------------------------------------------------------------

C*************
      subroutine setmodel(iq)
      USE global
      integer iq
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnap.h'
      real fval              ! function for computation of real values
      integer ii,i
      DO i=1,nparm
         IF (par2phy(i).eq.1) then
c***  this is the water depth
            hw=fval(model(i,iq),i)
c     write(*,*)i,iq,model(i,iq),fval(model(i,iq),i),hw
         elseif (par2phy(i).eq.2) then
c***  this is the water velocity
c     write(*,*)'fval', fval(model(i,iq),i)
            c0(par2lay(i))=fval(model(i,iq),i)
         elseif (par2phy(i).eq.3) then
c***  this is the sediment velocity
            c1(par2lay(i))=fval(model(i,iq),i)
         elseif (par2phy(i).eq.4) then
c***  this is the attenuation
            beta(par2lay(i))=fval(model(i,iq),i)
         elseif (par2phy(i).eq.5) then
c***  this is the roughness
            scatt(par2lay(i))=fval(model(i,iq),i)
         elseif (par2phy(i).eq.6) then
c***  this is the sediment density
            r1 = fval(model(i,iq),i)
         elseif (par2phy(i).eq.8) then
c***  source depth
            srd(par2lay(i))=fval(model(i,iq),i)
         elseif (par2phy(i).eq.9) then
c***  receiver range...
            rng(par2lay(i))= fval(model(i,iq),i)
         elseif (par2phy(i).eq.11) then
c***  this is the EOF
            aeof(par2lay(i)) = fval(model(i,iq),i)
         elseif (par2phy(i).eq.12) then
c***  this is the bottom P-velocity
            c2 = fval(model(i,iq),i)
         elseif (par2phy(i).eq.13) then
c***  this is the bottom S-velocity
            c2s = fval(model(i,iq),i)
         elseif (par2phy(i).eq.14) then
c***  this is the bottom density
            r2 = fval(model(i,iq),i)
         elseif (par2phy(i).eq.15) then
c***  this is the receiver depth
c     fldrd(1) = fval(model(i,iq),i)
c     fldrd(2) = fldrd(1)+fldrd(3)*(ndep-1)
            do ii=1,ndep
               rdep(ii)=rdref(ii)-rdref(1)+fval(model(i,iq),i)
            enddo
         elseif (par2phy(i).eq.16) then
c***  this is the velocity points in water
            z0(par2lay(i)) = fval(model(i,iq),i)
         elseif (par2phy(i).eq.17) then
c***  this is the sediment depth
            h1 = fval(model(i,iq),i)
         elseif (par2phy(i).eq.18) then
c***  this is the velocity points in sediment
            z1(par2lay(i)) = fval(model(i,iq),i)
         elseif (par2phy(i).eq.19) then
c***  this is tilt
            dtilt = fval(model(i,iq),i)
         elseif (par2phy(i).eq.20) then
c***  this is change of sound speed profile
c     write(*,*)'Using profile no',nint(fval(model(i,iq),i))
            do ii=1,nobspts
               c0(ii) = xobsssp(nint(fval(model(i,iq),i)),ii)
c     write(*,*)z0(ii),c0(ii)
            enddo
         elseif (par2phy(i).eq.21) then
c***  this is the arraysnape
            arrayshape(par2lay(i)) = fval(model(i,iq),i)
         elseif (par2phy(i).eq.22) then
c***  this is the source factor
            pfact(par2lay(i))=fval(model(i,iq),i)
         elseif (par2phy(i).eq.23) then
c***  this is the multistatic repeater
            multigeo(par2lay(i))=fval(model(i,iq),i)
         elseif (par2phy(i).eq.24) then
c***  this is for grain size     .cfh.
            phim=fval(model(i,iq),i)
         elseif (par2phy(i).eq.25) then
c***  this is for range offsets 
            dr_ind(par2lay(i))=fval(model(i,iq),i)
         elseif (par2phy(i).eq.26) then
c***  this is for time offsets 
            time_ind(par2lay(i))=fval(model(i,iq),i)
         ELSEIF (par2phy(i).EQ.28) THEN
c***  this is for time delay
            del_time = fval(model(i,iq),i)
         elseif (par2phy(i).eq.29) then
c***  this is for error variance .cfh.
            numh=fval(model(i,iq),i)
         else
            write(*,*) ' setmodel: option not defined for parameter',i
         endif
c     WRITE(prtfil,*)' setmodel:, i,fval',i,fval(model(i,iq),i)
      enddo
      end
c*************
      subroutine setmodelbest(modelb)
      USE global
      integer i,ii,modelb(*)
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnap.h'
      real fval              ! function for computation of real values
      do i=1,nparm
         if (par2phy(i).eq.1) then
c***  this is the water depth
            hw=fval(modelb(i),i)
         elseif (par2phy(i).eq.2) then
c***  this is the water velocity
c     write(*,*)'fval', fval(modelb(i),i)
            c0(par2lay(i))=fval(modelb(i),i)
         elseif (par2phy(i).eq.3) then
c***  this is the sediment velocity
            c1(par2lay(i))=fval(modelb(i),i)
         elseif (par2phy(i).eq.4) then
c***  this is the sediment attenuation
            beta(par2lay(i))=fval(modelb(i),i)
         elseif (par2phy(i).eq.5) then
c***  this is the roughness
            scatt(par2lay(i))=fval(modelb(i),i)
         elseif (par2phy(i).eq.6) then
c***  this is the sediment density
            r1 = fval(modelb(i),i)
         elseif (par2phy(i).eq.8) then
c***  source depth
            srd(par2lay(i))=fval(modelb(i),i)
         elseif (par2phy(i).eq.9) then
c***  receiver range...
            rng(par2lay(i))= fval(modelb(i),i)
         elseif (par2phy(i).eq.11) then
c***  this is the EOF
            aeof(par2lay(i)) = fval(modelb(i),i)
         elseif (par2phy(i).eq.12) then
c***  this is the bottom P-velocity
            c2 = fval(modelb(i),i)
         elseif (par2phy(i).eq.13) then
c***  this is the bottom S-velocity
            c2s = fval(modelb(i),i)
         elseif (par2phy(i).eq.14) then
c***  this is the bottom density
            r2 = fval(modelb(i),i)
         elseif (par2phy(i).eq.15) then
c***  this is the receiver depth
c     fldrd(1) = fval(modelb(i),i)
c     fldrd(2) = fldrd(1)+fldrd(3)*(ndep-1)
            do ii=1,ndep
               rdep(ii)=rdref(ii)-rdref(1)+fval(modelb(i),i)
            enddo
         elseif (par2phy(i).eq.16) then
c***  this is the velocity points in water
            z0(par2lay(i)) = fval(modelb(i),i)
         elseif (par2phy(i).eq.17) then
c***  this is the sediment depth
            h1 = fval(modelb(i),i)
         elseif (par2phy(i).eq.18) then
c***  this is the velocity points in sediment
            z1(par2lay(i)) = fval(modelb(i),i)
         elseif (par2phy(i).eq.19) then
c***  this is tilt
            dtilt = fval(modelb(i),i)
         elseif (par2phy(i).eq.20) then
c***  this is change of sound speed profile
            do ii=1,nobspts
               c0(ii) = xobsssp(nint(fval(modelb(i),i)),ii)
            enddo
         elseif (par2phy(i).eq.21) then
c***  this is the array shape
            arrayshape(par2lay(i)) = fval(modelb(i),i)
         elseif (par2phy(i).eq.22) then
c***  this is the source factor
            pfact(par2lay(i))=fval(modelb(i),i)
         elseif (par2phy(i).eq.23) then
c***  this is the multistatic repeater
            multigeo(par2lay(i))=fval(modelb(i),i)
         elseif (par2phy(i).eq.24) then
c***  this is for grain size
            phim=fval(modelb(i),i)
         elseif (par2phy(i).eq.25) then
c***  this is for range offsets
            dr_ind(par2lay(i))=fval(modelb(i),i)
         elseif (par2phy(i).eq.26) then
c***  this is for time offsets
            time_ind(par2lay(i))=fval(modelb(i),i)
         ELSEIF (par2phy(i).EQ.28) THEN
c***  this is for time delay
            del_time = fval(modelb(i),i)
         elseif (par2phy(i).eq.29) then
c***  this is for error variance
            numh=fval(modelb(i),i)
         else
            write(*,*) ' snapinter: option not defined'
         endif
      enddo
      end

      subroutine setmodelreal(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnap.h'
      integer i,ii
      real x(*)
      do i=1,nparm
         if (par2phy(i).eq.1) then
c***  this is the water depth
            hw=x(i)
         elseif (par2phy(i).eq.2) then
c***  this is the water velocity
c     write(*,*)'fval', fval(modelb(i),i)
            c0(par2lay(i))=x(i)
         elseif (par2phy(i).eq.3) then
c***  this is the sediment velocity
            c1(par2lay(i))=x(i)
         elseif (par2phy(i).eq.4) then
c***  this is the  attenuation
            beta(par2lay(i))=x(i)
         elseif (par2phy(i).eq.5) then
c***  this is the roughness
            scatt(par2lay(i))=x(i)
         elseif (par2phy(i).eq.6) then
c***  this is the sediment density
            r1 = x(i)
         elseif (par2phy(i).eq.8) then
c***  source depth
            srd(par2lay(i))=x(i)
         elseif (par2phy(i).eq.9) then
c***  receiver range...
            rng(par2lay(i))= x(i)
         elseif (par2phy(i).eq.11) then
c***  this is the EOF
            aeof(par2lay(i)) = x(i)
         elseif (par2phy(i).eq.12) then
c***  this is the bottom P-velocity
            c2 = x(i)
         elseif (par2phy(i).eq.13) then
c***  this is the bottom S-velocity
            c2s =x(i)
         elseif (par2phy(i).eq.14) then
c***  this is the bottom density
            r2 = x(i)
         elseif (par2phy(i).eq.15) then
c***  this is the receiver depth
            do ii=1,ndep
               rdep(ii)=rdref(ii)-rdref(1)+x(i)
            enddo
c     fldrd(1) = x(i)
c     fldrd(2) = fldrd(1)+fldrd(3)*(ndep-1)
         elseif (par2phy(i).eq.16) then
c***  this is the velocity points in water
            z0(par2lay(i)) = x(i)
         elseif (par2phy(i).eq.17) then
c***  this is the sediment depth
            h1 = x(i)
         elseif (par2phy(i).eq.18) then
c***  this is the velocity points in sediment
            z1(par2lay(i)) = x(i)
         elseif (par2phy(i).eq.19) then
c***  this is tilt
            dtilt = x(i)
         elseif (par2phy(i).eq.20) then
c***  this is change of sound speed profile
            do ii=1,nobspts
               c0(ii) = xobsssp(nint(x(i)),ii)
            enddo
         elseif (par2phy(i).eq.21) then
c***  this is the array shape
            arrayshape(par2lay(i)) = x(i)
         elseif (par2phy(i).eq.22) then
c***  this is the source factor
            pfact(par2lay(i))=x(i)
         elseif (par2phy(i).eq.23) then
c***  this is the multistatic repeater
            multigeo(par2lay(i))=x(i)
         elseif (par2phy(i).eq.24) then
c***  this is for grain size
            phim=x(i)
         elseif (par2phy(i).eq.25) then
c***  this is for range offsets
            dr_ind(par2lay(i))=x(i)
         elseif (par2phy(i).eq.26) then
c***  this is for timee offsets
            time_ind(par2lay(i))=x(i)
         ELSEIF (par2phy(i).EQ.28) THEN
c***  this is for time delay
            del_time = x(i)
         elseif (par2phy(i).eq.29) then
c***  this is for error variance
            numh=x(i)
         else
            write(*,*) ' snapinter: option not defined'
         endif
      enddo
      end
c     
c=======================================================================

      subroutine getmodelreal(x)

c=======================================================================

c------------------------------------------------------------------------
c     Map environment variables from ORCA to array x.
c------------------------------------------------------------------------
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnap.h'
      integer i
      real x(*)
      do i=1,nparm
         if (par2phy(i).eq.1) then
c***  this is the water depth
            x(i)=hw
         elseif (par2phy(i).eq.2) then
c***  this is the water velocity
            x(i)=c0(par2lay(i))
         elseif (par2phy(i).eq.3) then
c***  this is the sediment velocity
            x(i)= c1(par2lay(i))
         elseif (par2phy(i).eq.4) then
c***  this is the sediment attenuation
            x(i)=beta(par2lay(i))
         elseif (par2phy(i).eq.5) then
c***  this is the roughness
            x(i)=scatt(par2lay(i))
         elseif (par2phy(i).eq.6) then
c***  this is the sediment density
            x(i)=r1 
         elseif (par2phy(i).eq.8) then
c***  source depth
            x(i)=srd(par2lay(i))
         elseif (par2phy(i).eq.9) then
c***  receiver range...
            x(i)= rng(par2lay(i))
         elseif (par2phy(i).eq.11) then
c***  this is the EOF
            x(i)=aeof(par2lay(i)) 
         elseif (par2phy(i).eq.12) then
c***  this is the bottom P-velocity
            x(i) = c2 
         elseif (par2phy(i).eq.13) then
c***  this is the bottom S-velocity
            x(i) = c2s 
         elseif (par2phy(i).eq.14) then
c***  this is the bottom density
            x(i) = r2
         elseif (par2phy(i).eq.15) then
c***  this is the receiver depth
            x(i) = rdep(1) 
c     fldrd(2) = fldrd(1)+fldrd(3)*(ndep-1)
         elseif (par2phy(i).eq.16) then
c***  this is the velocity points in water
            x(i)= z0(par2lay(i))
         elseif (par2phy(i).eq.17) then
c***  this is the sediment depth
            x(i)=h1
         elseif (par2phy(i).eq.18) then
c***  this is the velocity points in sediment
            x(i)= z1(par2lay(i))
         elseif (par2phy(i).eq.19) then
c***  this is tilt
            x(i)=dtilt
         elseif (par2phy(i).eq.20) then
c***  this is change of sound speed profile
            x(i) = 10           ! This does realy not work
         elseif (par2phy(i).eq.21) then
c***  this is the array shape
            x(i)=arrayshape(par2lay(i)) 
         elseif (par2phy(i).eq.22) then
c***  this is the source factor
            x(i)=pfact(par2lay(i))
         elseif (par2phy(i).eq.23) then
c***  this is the multistatic repeater
            x(i)=multigeo(par2lay(i))
         elseif (par2phy(i).eq.24) then
c***  this is for grain size
            x(i)= phim  
         elseif (par2phy(i).eq.25) then
c***  this is for range offsets
            x(i)=dr_ind(par2lay(i))
         elseif (par2phy(i).eq.26) then
c***  this is for time offsets
            x(i)=time_ind(par2lay(i))
         ELSEIF (par2phy(i).EQ.28) THEN
c***  this is for time delay
            x(i)=del_time 
         elseif (par2phy(i).eq.29) then
c***  this is for error variance
            x(i)= numh
         else
            write(*,*) ' snapinter: option not defined'
         endif
      enddo
      end
c*********************
      subroutine setmodelx(ixlay,ixpar,theta)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnap.h'     
      integer ii,ixlay,ixpar
      real theta
      if (ixpar.eq.1) then
c***  this is the water depth
         hw=theta
      elseif (ixpar.eq.2) then
c***  this is the water velocity
c     write(*,*)'fval', theta
         c0(ixlay)=theta
      elseif (ixpar.eq.3) then
c***  this is the sediment velocity
         c1(ixlay)=theta
      elseif (ixpar.eq.4) then
c***  this is the attenuation
         beta(ixlay)=theta
      elseif (ixpar.eq.5) then
c***  this is the roughness
         scatt(ixlay)=theta
      elseif (ixpar.eq.6) then
c***  this is the sediment density
         r1 = theta
      elseif (ixpar.eq.8) then
c***  source depth
         srd(ixlay)=theta
      elseif (ixpar.eq.9) then
c***  receiver range...
         rng(ixlay)= theta
      elseif (ixpar.eq.11) then
c***  this is the EOF
         aeof(ixlay) = theta
      elseif (ixpar.eq.12) then
c***  this is the bottom P-velocity
         c2 = theta
      elseif (ixpar.eq.13) then
c***  this is the bottom S-velocity
         c2s =theta
      elseif (ixpar.eq.14) then
c***  this is the bottom density
         r2 = theta
      elseif (ixpar.eq.15) then
c***  this is the receiver depth
c     fldrd(1) = theta
         do ii=1,ndep
            rdep(ii)=rdref(ii)-rdref(1)+theta
         enddo
c     fldrd(2) = fldrd(1)+fldrd(3)*(ndep-1)
      elseif (ixpar.eq.16) then
c***  this is the velocity points in water
         z0(ixlay)= theta
      elseif (ixpar.eq.17) then
c***  this is the sediment depth
         h1=theta
      elseif (ixpar.eq.18) then
c***  this is the velocity points in sediment
         z1(ixlay)= theta
      elseif (ixpar.eq.19) then
c***  this is tilt
         dtilt = theta
      elseif (ixpar.eq.20) then
c***  this is change of sound speed profile
         do ii=1,nobspts
            c0(ii) = xobsssp(nint(theta),ii)
         enddo
      elseif (ixpar.eq.21) then
c***  this is the array shape
         arrayshape(ixlay) = theta
      elseif (ixpar.eq.22) then
c***  this is the source factor
         pfact(ixlay)=theta
      elseif (ixpar.eq.23) then
c***  this is the multistatic repeater
         multigeo(ixlay)=theta
      elseif (ixpar.eq.24) then
c***  this is for grain size
         phim=theta
      elseif (ixpar.eq.25) then
c***  this is the range offsets 
         dr_ind(ixlay)=theta
      elseif (ixpar.eq.26) then
c***  this is the time offsets 
         time_ind(ixlay)=theta
      ELSEIF (ixpar.EQ.28) THEN
c***  this is for time delay
         del_time =theta
      elseif (ixpar.eq.29) then
c***  this is for error variance
         numh=theta
      else
         write(*,*) ' snapinter: option not defined'
      endif
      end

c*************
      subroutine setmodeleof(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnap.h'
      integer i,ii
      real x(*)

c     write(*,*) ' snapinter: entering setmodeleof'
      do i=1,neofvar
         if (par2phy_eof(i).eq.1) then
c***  this is the water depth
            hw=x(i)
         elseif (par2phy_eof(i).eq.2) then
c***  this is the water velocity
c     write(*,*)'fval', fval(modelb(i),i)
            c0(par2lay_eof(i))=x(i)
         elseif (par2phy_eof(i).eq.3) then
c***  this is the sediment velocity
            c1(par2lay_eof(i))=x(i)
         elseif (par2phy_eof(i).eq.4) then
c***  this is the sediment attenuation
            beta(par2lay_eof(i))=x(i)
         elseif (par2phy_eof(i).eq.5) then
c***  this is the roughness
            scatt(par2lay_eof(i))=x(i)
         elseif (par2phy_eof(i).eq.6) then
c***  this is the sediment density
            r1 = x(i)
         elseif (par2phy_eof(i).eq.8) then
c***  source depth
            srd(par2lay_eof(i))=x(i)
         elseif (par2phy_eof(i).eq.9) then
c***  receiver range...
            rng(par2lay_eof(i))= x(i)
         elseif (par2phy_eof(i).eq.12) then
c***  this is the bottom P-velocity
            c2 = x(i)
         elseif (par2phy_eof(i).eq.13) then
c***  this is the bottom S-velocity
            c2s = x(i)
         elseif (par2phy_eof(i).eq.14) then
c***  this is the bottom density
            r2 = x(i)
         elseif (par2phy_eof(i).eq.15) then
c***  this is the receiver depth
c     fldrd(1) = x(i)
c     fldrd(2) = fldrd(1)+fldrd(3)*(ndep-1)
            do ii=1,ndep
               rdep(ii)=rdref(ii)-rdref(1)+x(i)
            enddo
         elseif (par2phy_eof(i).eq.16) then
c***  this is the velocity points in water
            z0(par2lay_eof(i)) = x(i)
         elseif (par2phy_eof(i).eq.17) then
c***  this is the sediment depth
            h1 = x(i)
         elseif (par2phy_eof(i).eq.18) then
c***  this is the velocity points in water
            z1(par2lay_eof(i)) = x(i)
         elseif (par2phy_eof(i).eq.19) then
c***  this is tilt
            dtilt = x(i)
         elseif (par2phy_eof(i).eq.21) then
c***  this is the array shape
            arrayshape(par2lay_eof(i)) = x(i)
         elseif (par2phy_eof(i).eq.22) then
c***  this is the source factor
            pfact(par2lay_eof(i))=x(i)
         elseif (par2phy_eof(i).eq.23) then
c***  this is the multistatic repeater
            multigeo(par2lay_eof(i))=x(i)
         elseif (par2phy_eof(i).eq.24) then
c***  this is for grain size
            phim=x(i)
         elseif (par2phy_eof(i).eq.25) then
c***  this is for range offsets
            dr_ind(par2lay_eof(i))=x(i)
         elseif (par2phy_eof(i).eq.26) then
c***  this is for time offsets
            time_ind(par2lay_eof(i))=x(i)
         ELSEIF (par2phy_eof(i).EQ.28) THEN
c***  this is for time delay
            del_time=x(i) 
         elseif (par2phy_eof(i).eq.29) then
c***  this is for error variance
            numh=x(i)
        else
            write(*,*) ' snapinter: option not defined'
         endif
      enddo
c     write(*,*) ' snapinter: exiting setmodeleof ...'
      end


