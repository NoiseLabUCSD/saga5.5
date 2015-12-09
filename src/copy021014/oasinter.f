      
c------------------------------------------------------------------------
c     Map environment variables to OASE
c------------------------------------------------------------------------

C*************
      subroutine setmodel(iq)
      USE global
      integer iq
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      real fval                 ! function for computation of real values
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnrd.f'     
      integer i,j,ii
      real  deltathick
      do i=1,nparm
         if (par2phy(i).eq.7) then
c***  this is the thickness of each layer
c     write(*,*)'fval', fval(model(i,iq),i)
            deltathick= 
     &           fval(model(i,iq),i)-
     &           (v(par2lay(i)+1,1)-v(par2lay(i),1))
c     write(*,*)'deltathick', deltathick 
            do j=par2lay(i)+1,nlay
               v(j,1)=v(j,1)+deltathick
            enddo
         elseif (par2phy(i).eq.8) then
c***  source depth...
            sdc(1)= fval(model(i,iq),i)
         elseif (par2phy(i).eq.9) then
c***  receiver range...
            xranges(1)= fval(model(i,iq),i)
         elseif (par2phy(i).eq.11) then
c***  this is the EOF
            aeof(par2lay(i)) = fval(model(i,iq),i)
         elseif (par2phy(i).eq.15) then
c***  receiver depth...
            do ii=ndep,1,-1
               rdep(ii)=rdep(ii)-rdep(1)+fval(model(i,iq),i)
            enddo
         else
            v(par2lay(i),par2phy(i))=fval(model(i,iq),i)
         endif
      enddo
c     write(*,*)'parameters'
c     write(*,*)iq,(par2phy(i),model(i,iq),i=1,nparm)      
c     write(*,'(8f10.3)')(fval(model(i,iq),i),i=1,nparm)  
c     i=1
c     DO j=0,7
c     write(*,*) 'fval',fval(j,i)
c     ENDDO
c     
c     write(*,'(8f10.3)')(v(par2lay(i),par2phy(i)),i=1,nparm)      
c     do i=1,nlay
c     write(*,'(6f12.3)')(v(i,j),j=1,6)
c     ENDDO
      end
c*************
      subroutine setmodelbest(modelb)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      real fval                 ! function for computation of real values
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnrd.f'     
      integer i,ii,modelb(*)
      real  deltathick
      do i=1,nparm
         if (par2phy(i).eq.7) then
c***  this is the thickness of each layer
c     write(*,*)'fval', fval(model(i,iq),i)
            deltathick= 
     &           fval(modelb(i),i)-(v(par2lay(i)+1,1)-v(par2lay(i),1))
c     write(*,*)'deltathick', deltathick 
            do j=par2lay(i)+1,nlay
               v(j,1)=v(j,1)+deltathick
            enddo
         elseif (par2phy(i).eq.8) then
c***  source depth...
            sdc(1)= fval(modelb(i),i)
         elseif (par2phy(i).eq.9) then
c***  receiver range...
            xranges(1)= fval(modelb(i),i)
         elseif (par2phy(i).eq.11) then
c***  this is the EOF
            aeof(par2lay(i)) = fval(modelb(i),i)
         elseif (par2phy(i).eq.15) then
c***  receiver depth...
            do ii=ndep,1,-1
               rdep(ii)=rdep(ii)-rdep(1)+fval(modelb(i),i)
            enddo
         else
            v(par2lay(i),par2phy(i))=fval(modelb(i),i)
         endif
      enddo
      end
      subroutine setmodelreal(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnrd.f'     
      integer i,ii
      real x(*)
      real  deltathick
      do i=1,nparm
         if (par2phy(i).eq.7) then
c***  this is the thickness of each layer
            deltathick= 
     &           x(i)-(v(par2lay(i)+1,1)-v(par2lay(i),1))
            do j=par2lay(i)+1,nlay
               v(j,1)=v(j,1)+deltathick
            enddo
         elseif (par2phy(i).eq.8) then
c***  source depth...
            sdc(1)= x(i)
         elseif (par2phy(i).eq.9) then
c***  receiver range...
            xranges(1)= x(i)
         elseif (par2phy(i).eq.11) then
c***  this is the EOF
            aeof(par2lay(i)) = x(i)
         elseif (par2phy(i).eq.15) then
c***  receiver depth...
            do ii=ndep,1,-1
               rdep(ii)=rdep(ii)-rdep(1)+x(i)
            enddo
         else
            v(par2lay(i),par2phy(i))=x(i)
         endif
      enddo
      end
c********************
      subroutine getmodelreal(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnrd.f'     
      integer i
      real x(*)
c     write(*,*)' oasinter,sdc',sdc(1)
      do i=1,nparm
         if ( par2phy(i).eq.7) then
            x(i)=(v(par2lay(i)+1,1)-v(par2lay(i),1))
         elseif ( par2phy(i).eq.8) then
            x(i)=sdc(1)
         elseif (par2phy(i).eq.9) then
c***  receiver range...
            x(i)=xranges(1)
         elseif (par2phy(i).eq.11) then
c***  this is the EOF
            x(i)=aeof(par2lay(i)) 
         elseif (par2phy(i).eq.15) then
c***  receiver depth...
            x(i)=rdep(1)
         else
            x(i)=v(par2lay(i),par2phy(i))
         endif
      enddo
      end
c*********************
      subroutine setmodelx(ixlay,ixpar,theta)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnrd.f'     
      integer ixlay,ixpar,ii
      real theta
      real  deltathick
      if (ixpar.eq.7) then
c***  this is the thickness of each layer
c     write(*,*)'fval', fval(model(i,iq),i)
         deltathick= 
     &        theta-(v(ixlay+1,1)-v(ixlay,1))
c     write(*,*)'deltathick', deltathick 
         do j=ixlay+1,nlay
            v(j,1)=v(j,1)+deltathick
         enddo
      elseif (ixpar.eq.8) then
c***  source depth...
         sdc(1)= theta
      elseif (ixpar.eq.9) then
c***  receiver range...
         xranges(1)= theta
      elseif (ixpar.eq.11) then
c***  this is the EOF
         aeof(ixlay) = theta
      elseif (ixpar.eq.15) then
c***  receiver depth...
         do ii=ndep,1,-1
            rdep(ii)=rdep(ii)-rdep(1)+theta
         enddo
      else
         v(ixlay,ixpar)=theta 
      endif
      
      end

      subroutine setmodeleof(x)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnrd.f'     
      integer i,ii
      real x(*)
      real  deltathick
      do i=1,neofvar
         if (par2phy_eof(i).eq.7) then
c***  this is the thickness of each layer
c     write(*,*)'fval', fval(model(i,iq),i)
            deltathick= 
     &           x(i)-(v(par2lay_eof(i)+1,1)-v(par2lay_eof(i),1))
c     write(*,*)'deltathick', deltathick 
            do j=par2lay_eof(i)+1,nlay
               v(j,1)=v(j,1)+deltathick
            enddo
         elseif (par2phy_eof(i).eq.8) then
c***  source depth...
            sdc(1)= x(i)
         elseif (par2phy_eof(i).eq.9) then
c***  receiver range...
            xranges(1)= x(i)
         elseif (par2phy_eof(i).eq.15) then
c***  receiver depth...
            do ii=ndep,1,-1
               rdep(ii)=rdref(ii)-rdref(1)+x(i)
            enddo
         else
            v(par2lay_eof(i),par2phy_eof(i))=x(i)
         endif
      enddo
      end
