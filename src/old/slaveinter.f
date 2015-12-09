C*************
      subroutine setmodel(iq)
      integer iq
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comslave.h'
      real fval              ! function for computation of real values
      integer ii,i,j
      real  deltathick
      do i=1,nparm
        if (par2phy(i).eq.1) then
c***    this is the slave parameters
           slavepar(par2lay(i))=fval(model(i,iq),i)
c           write(*,*)i,iq,model(i,iq),fval(model(i,iq),i),hw
        elseif (par2phy(i).eq.11) then
c***    this is the EOF
           aeof(par2lay(i)) = fval(model(i,iq),i)
        else
          write(*,*) ' setmodel: option not defined for parameter',i
        endif
c        WRITE(prtfil,*)' setmodel:, i,fval',i,fval(model(i,iq),i)
      enddo
c      write(*,*)'parameters'
c      write(*,'(8f10.3)')(v(par2lay(i),par2phy(i)),i=1,nparm)      
c      do i=1,nlay
c        write(*,'(6f12.3)')(v(i,j),j=1,6)
c      ENDDO
      end
c*************
      subroutine setmodelbest(modelb)
      integer i,ii,modelb(*)
      real  deltathick
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comslave.h'
      real fval              ! function for computation of real values
      do i=1,nparm
        if (par2phy(i).eq.1) then
c***    this is the slave parameters
           slavepar(par2lay(i))=fval(modelb(i),i)
c           write(*,*)i,iq,model(i,iq),fval(model(i,iq),i),hw
        elseif (par2phy(i).eq.11) then
c***    this is the EOF
           aeof(par2lay(i)) = fval(modelb(i),i)
        else
          write(*,*) ' snapinter: option not defined'
        endif
      enddo
      end

      subroutine setmodelreal(x)
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comslave.h'
      real fval              ! function for computation of real values
      integer i,ii
      real  x(*)
      do i=1,nparm
        if (par2phy(i).eq.1) then
c***    this is the slave parameters
           slavepar(par2lay(i))=x(i)
c           write(*,*)i,iq,model(i,iq),fval(model(i,iq),i),hw
        elseif (par2phy(i).eq.11) then
c***    this is the EOF
           aeof(par2lay(i)) = x(i)
        else
          write(*,*) ' snapinter: option not defined'
        endif
      enddo
      end
c********************
      subroutine getmodelreal(x)
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comslave.h'
      real fval              ! function for computation of real values
      integer i,ii
      real x(*)
      do i=1,nparm
        if (par2phy(i).eq.1) then
c***    this is the slave parameters
           x(i)=slavepar(par2lay(i))
c           write(*,*)i,iq,model(i,iq),fval(model(i,iq),i),hw
        elseif (par2phy(i).eq.11) then
c***    this is the EOF
           x(i)=aeof(par2lay(i)) 
        else
          write(*,*) ' snapinter: option not defined'
        endif
      enddo
      end
c*********************
      subroutine setmodelx(ixlay,ixpar,theta)
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comslave.h'     
      real fval              ! function for computation of real values
      integer ii,ixlay,ixpar
      real theta
        if (ixpar.eq.1) then
c***    this is the slave parameters
           slavepar(ixlay)=theta
c           write(*,*)i,iq,model(i,iq),fval(model(i,iq),i),hw
          elseif (ixpar.eq.11) then
c***    this is the EOF
           aeof(ixlay) = theta
        else
          write(*,*) ' snapinter: option not defined'
        endif
      end

c*************
      subroutine setmodeleof(x)
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comslave.h'
      real fval              ! function for computation of real values
      integer i
      real x(*)
      real  deltathick
      do i=1,neofvar
        if (par2phy_eof(i).eq.1) then
c***    this is the water depth
            slavepar(par2lay_eof(i))=x(i)
        else
          write(*,*) ' snapinter: option not defined'
        endif
      enddo
      end


