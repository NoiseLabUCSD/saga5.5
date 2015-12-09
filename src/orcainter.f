      
c------------------------------------------------------------------------
c     Map environment variables to ORCA
c------------------------------------------------------------------------

      subroutine  setmodel(iq)
      USE global
      use fgs_com
      use i_o_com
      implicit none
      INCLUDE 'comopt.h'
      include 'comforw.h'
      real fval                 ! function for computation of real values
      integer iq

      integer n,ipar,ilay,idx2,kd

      do i=1,nparm
         ipar=par2phy(i)
         ilay=par2lay(i)
         idx2=par3(i)

         select case( ipar )

      case( 1 )                 ! Water Depth
         zsvp(nsvp) = fval(model(i,iq),i)
         
      case( 2 )                 ! P Speed		 
         if( idx2==0 ) geob(1:2,1,ilay) =  fval(model(i,iq),i)
         if( idx2==1 ) geob(1,1,ilay) =  fval(model(i,iq),i)
         if( idx2==2 ) geob(2,1,ilay) =  fval(model(i,iq),i)

      case( 3 )                 ! S speed		 
         if( idx2==0 ) geob(1:2,2,ilay) = fval(model(i,iq),i)
         if( idx2==1 ) geob(1,2,ilay) = fval(model(i,iq),i)
         if( idx2==2 ) geob(2,2,ilay) =  fval(model(i,iq),i)

      case( 4 )                 ! P Attenuation
         if( idx2==0 ) geob(1:2,4,ilay) = fval(model(i,iq),i)
         if( idx2==1 ) geob(1,4,ilay) =  fval(model(i,iq),i)
         if( idx2==2 ) geob(2,4,ilay) = fval(model(i,iq),i)

      case( 5 )                 ! S Attenuation 
         if( idx2==0 ) geob(1:2,5,ilay) =  fval(model(i,iq),i)
         if( idx2==1 ) geob(1,5,ilay) =  fval(model(i,iq),i)
         if( idx2==2 ) geob(2,5,ilay) = fval(model(i,iq),i)

      case( 6 )                 ! Density	
         if( idx2==0 ) geob(1:2,3,ilay) =  fval(model(i,iq),i)
         if( idx2==1 ) geob(1,3,ilay) = fval(model(i,iq),i)
         if( idx2==2 ) geob(2,3,ilay) =  fval(model(i,iq),i)

      case( 7 )                 ! Layer Thickness	 
         hb(ilay) =	 fval(model(i,iq),i)
         
      case( 8 )                 ! Source Depth
         zsrc(1) =  fval(model(i,iq),i)

      case( 9 )                 ! Source Range
         rkm(1) = fval(model(i,iq),i)

      case(11)                  ! this is the EOF
         aeof(par2lay(i)) = fval(model(i,iq),i)

      case( 15 )                ! First Receiver Depth
         zrec(1) = fval(model(i,iq),i)

      case( 19 )                ! Vertical Array Tilt
         dtilt =  fval(model(i,iq),i)
         
      case( 20 )                ! Sound Speed In Water
         csvp(ilay) = fval(model(i,iq),i)

      case( 29 )                ! MH error
         numh = fval(model(i,iq),i)

      case default
      end select
c     write(*,*)'setmodel orcainiter,i,ipar,ilay,idx2 x(i)',
c     1            i,ipar,ilay,idx2, fval(model(i,iq),i)

      enddo

      return
      end

c*************

      subroutine  setmodelbest(modelb)
      USE global
      use fgs_com
      use i_o_com
      implicit none
      INCLUDE 'comopt.h'
      include 'comforw.h'
      integer iq
      real fval                 ! function for computation of real values
      integer modelb(*)
      integer n,ivar,ipar,ilay,idx2,kd

      do ivar=1,nparm

         ipar=par2phy(ivar)
         ilay=par2lay(ivar)
         idx2=par3(ivar)

         select case( ipar )

      case( 1 )                 ! Water Depth
         zsvp(nsvp) = fval(modelb(i),ivar)
         
      case( 2 )                 ! P Speed		 
         if( idx2==0 ) geob(1:2,1,ilay) =  fval(modelb(i),ivar)
         if( idx2==1 ) geob(1,1,ilay)   =  fval(modelb(i),ivar)
         if( idx2==2 ) geob(2,1,ilay)   =  fval(modelb(i),ivar)

      case( 3 )                 ! S speed		 
         if( idx2==0 ) geob(1:2,2,ilay) = fval(modelb(i),ivar)
         if( idx2==1 ) geob(1,2,ilay)   = fval(modelb(i),ivar)
         if( idx2==2 ) geob(2,2,ilay)   = fval(modelb(i),ivar)

      case( 4 )                 ! P Attenuation
         if( idx2==0 ) geob(1:2,4,ilay) = fval(modelb(i),ivar)
         if( idx2==1 ) geob(1,4,ilay)   = fval(modelb(i),ivar)
         if( idx2==2 ) geob(2,4,ilay)   = fval(modelb(i),ivar)

      case( 5 )                 ! S Attenuation 
         if( idx2==0 ) geob(1:2,5,ilay) = fval(modelb(i),ivar)
         if( idx2==1 ) geob(1,5,ilay)   = fval(modelb(i),ivar)
         if( idx2==2 ) geob(2,5,ilay)   = fval(modelb(i),ivar)

      case( 6 )                 ! Density	
         if( idx2==0 ) geob(1:2,3,ilay) = fval(modelb(i),ivar)
         if( idx2==1 ) geob(1,3,ilay)   = fval(modelb(i),ivar)
         if( idx2==2 ) geob(2,3,ilay)   = fval(modelb(i),ivar)

      case( 7 )                 ! Layer Thickness	 
         hb(ilay) =	 fval(modelb(i),ivar)
         
      case( 8 )                 ! Source Depth
         zsrc(1) =  fval(modelb(i),ivar)

      case( 9 )                 ! Source Range
         rkm(1) = fval(modelb(i),ivar)

      case( 11 )   
         aeof(par2lay(i)) = fval(modelb(i),ivar)

      case( 15 )                ! First Receiver Depth
         zrec(1) = fval(modelb(i),ivar)

      case( 19 )                ! Vertical Array Tilt
         dtilt =  fval(modelb(i),i)
         
      case( 20 )                ! Sound Speed In Water
         csvp(ilay) = fval(modelb(i),ivar)

      case( 29 )                ! MH error
         numh = fval(modelb(i),ivar)

      case default
      end select

      enddo

      return
      end

c*************
      subroutine  setmodelreal(x)
      USE global
      use fgs_com
      use i_o_com
      implicit none
      INCLUDE 'comopt.h'
      include 'comforw.h'

      integer n,idim,ipar,ilay,idx2,kd
      real    x(*) 

      do idim=1,nparm

         ipar=par2phy(idim)
         ilay=par2lay(idim)
         idx2=par3(idim)

         select case( ipar )

      case( 1 )                 ! Water Depth
         zsvp(nsvp) = x(idim)
         
      case( 2 )                 ! P Speed		 
         if( idx2==0 ) geob(1:2,1,ilay) = x(idim)
         if( idx2==1 ) geob(1,1,ilay) = x(idim)
         if( idx2==2 ) geob(2,1,ilay) = x(idim)

      case( 3 )                 ! S speed		 
         if( idx2==0 ) geob(1:2,2,ilay) = x(idim)
         if( idx2==1 ) geob(1,2,ilay) = x(idim)
         if( idx2==2 ) geob(2,2,ilay) = x(idim)

      case( 4 )                 ! P Attenuation
         if( idx2==0 ) geob(1:2,4,ilay) = x(idim)
         if( idx2==1 ) geob(1,4,ilay) = x(idim)
         if( idx2==2 ) geob(2,4,ilay) = x(idim)

      case( 5 )                 ! S Attenuation 
         if( idx2==0 ) geob(1:2,5,ilay) = x(idim)
         if( idx2==1 ) geob(1,5,ilay) = x(idim)
         if( idx2==2 ) geob(2,5,ilay) = x(idim)

      case( 6 )                 ! Density	
         if( idx2==0 ) geob(1:2,3,ilay) = x(idim)
         if( idx2==1 ) geob(1,3,ilay) = x(idim)
         if( idx2==2 ) geob(2,3,ilay) = x(idim)

      case( 7 )                 ! Layer Thickness	 
         hb(ilay) = x(idim)	
         
      case( 8 )                 ! Source Depth
         zsrc(1) = x(idim)

      case( 9 )                 ! Source Range
         rkm(1) = x(idim)
 
      case( 11 )                ! this is the EOF
         aeof(par2lay(idim)) = x(idim)

      case( 15 )                ! First Receiver Depth
         zrec(1) = x(idim)

      case( 19 )                ! Vertical Array Tilt
         dtilt = x(idim)
         
      case( 20 )                ! Sound Speed In Water
         csvp(ilay) = x(idim)

      case( 29 )                ! MH error
         numh =  x(idim)

      case default
      end select

      enddo

      return
      end

c     
c=======================================================================

      subroutine  getmodelreal(x)

c=======================================================================

c------------------------------------------------------------------------
c     Map environment variables from ORCA to array x.
c------------------------------------------------------------------------

      USE global
      use fgs_com
      use i_o_com
      implicit none
      INCLUDE 'comopt.h'
      include 'comforw.h'
      
      integer n,idim,ipar,ilay,idx2
      real    x(*)
      
      do idim=1,nparm
         ipar=par2phy(idim)
         ilay=par2lay(idim)
         idx2=par3(idim)

         select case( ipar )

      case( 1 )                 ! Water Depth
         x(idim) = zsvp(nsvp)
         
      case( 2 )                 ! P Speed		 
         if( idx2==0 ) x(idim) = geob(1,1,ilay)
         if( idx2==1 ) x(idim) = geob(1,1,ilay)
         if( idx2==2 ) x(idim) = geob(2,1,ilay)

      case( 3 )                 ! S Speed	
         if( idx2==0 ) x(idim) = geob(1,2,ilay)
         if( idx2==1 ) x(idim) = geob(1,2,ilay)
         if( idx2==2 ) x(idim) = geob(2,2,ilay)

      case( 4 )                 ! P Attenuation
         if( idx2==0 ) x(idim) = geob(1,4,ilay)
         if( idx2==1 ) x(idim) = geob(1,4,ilay)
         if( idx2==2 ) x(idim) = geob(2,4,ilay)
         
      case( 5 )                 ! S Attenuation
         if( idx2==0 ) x(idim) = geob(1,5,ilay)
         if( idx2==1 ) x(idim) = geob(1,5,ilay)
         if( idx2==2 ) x(idim) = geob(2,5,ilay)

      case( 6 )                 ! Density
         if( idx2==0 ) x(idim) = geob(1,3,ilay)
         if( idx2==1 ) x(idim) = geob(1,3,ilay)
         if( idx2==2 ) x(idim) = geob(2,3,ilay)

      case( 7 )                 ! Layer Thickness
         x(idim) = hb(ilay)

      case( 8 )                 ! Source Depth
         x(idim) = zsrc(1)

      case( 9 )                 ! Source Range
         x(idim) = rkm(1)

      case(11)                  ! this is the EOF
         x(idim) = aeof(par2lay(idim)) 

      case( 15 )                ! First Receiver Depth
         x(idim) = zrec(1)

      case( 19 )                ! Vertical Array Tilt
         x(idim) = dtilt

      case( 20 )                ! Sound Speed In Water
         x(idim) = csvp(ilay)

      case( 29 )                ! MH error
         x(idim) = numh

      case default
      end select
c     write(*,*)'orcainiter,idim,ipar,ilay,idx2 x(idim)',
c     1            idim,ipar,ilay,idx2, x(idim)
      enddo
      
      return
      end
c************************************
      subroutine  setmodelx( ilay,ipar,idx2,theta)
      USE global
      use fgs_com
      use i_o_com
      implicit none
      INCLUDE 'comopt.h'
      include 'comforw.h'

      integer n,idim,ipar,ilay,idx2,kd
      real    theta 

c     write(*,*)'setmodelx( ipar,ilay,idx2,theta)'
c     write(*,*) ipar,ilay,idx2,theta

      select case( ipar )

      case( 1 )                 ! Water Depth
         zsvp(nsvp) = theta
         
      case( 2 )                 ! P Speed		 
         if( idx2==0 ) geob(1:2,1,ilay) = theta
         if( idx2==1 ) geob(1,1,ilay) = theta
         if( idx2==2 ) geob(2,1,ilay) = theta

      case( 3 )                 ! S speed		 
         if( idx2==0 ) geob(1:2,2,ilay) = theta
         if( idx2==1 ) geob(1,2,ilay) = theta
         if( idx2==2 ) geob(2,2,ilay) = theta

      case( 4 )                 ! P Attenuation
         if( idx2==0 ) geob(1:2,4,ilay) = theta
         if( idx2==1 ) geob(1,4,ilay) = theta
         if( idx2==2 ) geob(2,4,ilay) = theta

      case( 5 )                 ! S Attenuation 
         if( idx2==0 ) geob(1:2,5,ilay) = theta
         if( idx2==1 ) geob(1,5,ilay) = theta
         if( idx2==2 ) geob(2,5,ilay) = theta

      case( 6 )                 ! Density	
         if( idx2==0 ) geob(1:2,3,ilay) = theta
         if( idx2==1 ) geob(1,3,ilay) = theta
         if( idx2==2 ) geob(2,3,ilay) = theta

      case( 7 )                 ! Layer Thickness	 
         hb(ilay) = theta	
         
      case( 8 )                 ! Source Depth
         zsrc(1) = theta

      case( 9 )                 ! Source Range
         rkm(1) = theta

      case(11)                  ! this is the EOF
         aeof(ilay) = theta 

      case( 15 )                ! First Receiver Depth
         zrec(1) = theta

      case( 19 )                ! Vertical Array Tilt
         dtilt = theta
         
      case( 20 )                ! Sound Speed In Water
         csvp(ilay) = theta

      case( 29 )                ! MH error
         numh =  theta
      case default
      end select

      return
      end
c     

c*************
      subroutine  setmodeleof(x)
      USE global
      use fgs_com
      use i_o_com
      implicit none
      INCLUDE 'comopt.h'
      include 'comforw.h'

      integer n,idim,ipar,ilay,idx2,kd
      real    x(*) 

      Do idim=1,neofvar

         ipar=par2phy_eof(idim)
         ilay=par2lay_eof(idim)
         idx2=par3_eof(idim)
C         write(*,*) ipar, ilay,idx2
         select case( ipar )

      case( 1 )                 ! Water Depth
         zsvp(nsvp) = x(idim)
         
      case( 2 )                 ! P Speed		 
         if( idx2==0 ) geob(1:2,1,ilay) = x(idim)
         if( idx2==1 ) geob(1,1,ilay) = x(idim)
         if( idx2==2 ) geob(2,1,ilay) = x(idim)

      case( 3 )                 ! S speed		 
         if( idx2==0 ) geob(1:2,2,ilay) = x(idim)
         if( idx2==1 ) geob(1,2,ilay) = x(idim)
         if( idx2==2 ) geob(2,2,ilay) = x(idim)

      case( 4 )                 ! P Attenuation
         if( idx2==0 ) geob(1:2,4,ilay) = x(idim)
         if( idx2==1 ) geob(1,4,ilay) = x(idim)
         if( idx2==2 ) geob(2,4,ilay) = x(idim)

      case( 5 )                 ! S Attenuation 
         if( idx2==0 ) geob(1:2,5,ilay) = x(idim)
         if( idx2==1 ) geob(1,5,ilay) = x(idim)
         if( idx2==2 ) geob(2,5,ilay) = x(idim)

      case( 6 )                 ! Density	
         if( idx2==0 ) geob(1:2,3,ilay) = x(idim)
         if( idx2==1 ) geob(1,3,ilay) = x(idim)
         if( idx2==2 ) geob(2,3,ilay) = x(idim)

      case( 7 )                 ! Layer Thickness	 
         hb(ilay) = x(idim)	
         
      case( 8 )                 ! Source Depth
         zsrc(1) = x(idim)

      case( 9 )                 ! Source Range
         rkm(1) = x(idim)

      case( 15 )                ! First Receiver Depth
         zrec(1) = x(idim)

      case( 19 )                ! Vertical Array Tilt
         dtilt = x(idim)
         
      case( 20 )                ! Sound Speed In Water
         csvp(ilay) = x(idim)
         csvp(ilay) = x(idim)

      case( 29 )                ! MH error
         numh =  x(idim)

      case default
      end select
      enddo
C      write(*,*) (geob(idim,1,1),idim=1,2)
C      write(*,*) (geob(idim,1,2),idim=1,2)

      return
      end

