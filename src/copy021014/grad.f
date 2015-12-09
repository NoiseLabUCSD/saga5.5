c****************************************************************************
c     Subroutines for performing gauss-Newton inversion
c     
c     PETER GERSTOFT, 1992 
c****************************************************************************
      SUBROUTINE gaunew(theta,maxiter)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comgrad.h'
C     
C**   local variables
      INTEGER J
      REAL fitn
c-----PARAMETER vector
      REAL theta(mpar)
      REAL thetain(mpar)
      INTEGER maxiter
      COMPLEX a1(mpar*mpar),a2(mpar*mpar)
c---  initialize the starting point
      IF (iopt(4).EQ.2) THEN
         WRITE(*,*)'The starting point  is..'
         DO j=1,nparm
            thetain(j)=theta(j) 
            WRITE(*,*)' ',j,theta(j)
         ENDDO
      ENDIF
C     

      IF (iopt(22).EQ.1) THEN 
         IF (iopt(13).EQ.1) THEN
            IF (nparm.GE.2) THEN
               CALL crrao(thetain,nparm,a1) ! nparm=2 for hess. pl.
            ENDIF

            WRITE(*,*)' calling crline'
            DO j=1,nparm 
               thetain(j)=theta(j)
               WRITE(*,*)' ',j,theta(j)
            ENDDO
            CALL setmodelreal(theta)
            CALL crline(thetain,1,a1,a2) ! nparm=2 for hess. pl.
         ELSE
            WRITE(*,*)' routines found in the directory  old'
c     CALL unceir(nx,ndep,theta,nparm)
c     CALL plhess(nx,ndep,theta,2) ! nparm=2 for hess. pl.
         ENDIF
         RETURN
      ENDIF


      IF (iopt(25).EQ.1) THEN
c     WRITE(*,*) 'gnmin w/o phase'
         CALL gnmin(nx,ndep,theta,maxiter)
      ELSE
         WRITE(*,*) 'gnmin w phase is not developped'
      ENDIF

      IF (iopt(4).EQ.2) THEN
         CALL cost(fitn)
         WRITE(*,*)'The best  energy',fitn          
         DO j=1,nparm
            WRITE(prtfil,*)' ',j,theta(j)
            WRITE(*,*)' ',j,theta(j)
         ENDDO
      ENDIF 

      END

c***************************************
      SUBROUTINE jacobi(theta)  
c-    This computes the forward dirivatives of the response.
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comgrad.h'

      REAL theta(Mpar),xtemp		 ! the current parameter vector
      COMPLEX*16 chelp(Mpar)
      INTEGER i,j,index,ix,ifrq,jdep
c---  set the model and CALL the modelling routine
      REAL*8 rnorm
      CALL setmodelrealjac(theta)
      CALL forw2grad

c
c---  move the gradient  for shape functions
c       
      IF (iopt(17).EQ.1) THEN
         IF (nfrq.GT.1) THEN
            CALL movegrad
         ELSE
            STOP 'Shapefunctions for gradient only for one freq'
         ENDIF
      ENDIF

c Frequecy loop:

       DO 200 ifrq=1,nfrq      
cc
c---- IF we have to USE weighting
      IF (iopt(3).GE.3) THEN
            DO ix=1,nx            
               DO jdep=1,ndep
                  j=(ifrq-1)*ndep+jdep
                  index=(j-1)*nx
                  resp(ix+index)=resp(ix+index)*weight(ix+index)
                  DO i=1,nthet
                     grad(ix+index,i)=
     .                    grad(ix+index,i)*weight(ix+index) ! normalized
               ENDDO
            ENDDO
         ENDDO
      ENDIF
c
c---  we are ONLY using amplitude information
c
      IF (iopt(25).EQ.1) THEN     ! normalizing to real
        DO i=1,nthet
            DO ix=1,nx            
               DO jdep=1,ndep
                  j=(ifrq-1)*ndep+jdep
                  index=(j-1)*nx
                  grad(ix+index,i)=
     1                 (REAL(grad(ix+index,i))*REAL(resp(ix+index))
     2                 +AIMAG(grad(ix+index,i))*AIMAG(resp(ix+index)))
     3                 /ABS(resp(ix+index))
            ENDDO
         ENDDO
        ENDDO
c---    for response
            DO ix=1,nx            
               DO jdep=1,ndep
                  j=(ifrq-1)*ndep+jdep
                  index=(j-1)*nx
           resp(ix+index)=ABS(resp(ix+index))
          ENDDO
        ENDDO
      ENDIF                             ! end transformation to real
c
c---- for unit-vector optimization
c
      IF ((iopt(13).EQ.-1).OR.(iopt(27).EQ.1)) THEN
c---  for response
         rnorm=0.
         DO ix=1,nx            
            DO jdep=1,ndep
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               rnorm = rnorm 
     1              +REAL(resp(ix+index))**2 +AIMAG(resp(ix+index))**2
            ENDDO
         ENDDO
c 
         rnorm=1./SQRT(rnorm)
c---- normalize gradient and resp
         DO ix=1,nx            
            DO jdep=1,ndep
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               resp(ix+index)=resp(ix+index)*rnorm
               DO i=1,nthet
                  grad(ix+index,i)=grad(ix+index,i)*rnorm
               ENDDO
            ENDDO
         ENDDO
c---    axilulary matrix
         DO i=1,nthet
            chelp(i)=0.
            DO ix=1,nx            
               DO jdep=1,ndep
                  j=(ifrq-1)*ndep+jdep
                  index=(j-1)*nx
                  chelp(i)=chelp(i)
     .                 +CONJG(resp(ix+index))*grad(ix+index,i)
               ENDDO
            ENDDO
         ENDDO
c----- 
         DO ix=1,nx            
            DO jdep=1,ndep
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nthet
                  grad(ix+index,i)=
     1                 grad(ix+index,i)-resp(ix+index)*dreal(chelp(i))
               ENDDO
            ENDDO
         ENDDO
      ENDIF                     ! unit normalization
c-----------------
      DO i=1,nthet
c>>>>>>>>>>>>>> This normalisation is based on Lines and Treitel 
c        xtemp=0
c        DO k=1,ndep
c          index=(k-1)*nx
c          DO j=1,nx
c            xtemp=xtemp+ABS(grad(j+index,i))**2
c          ENDDO
c        ENDDO
c        xtemp=SQRT(xtemp)
c        IF (xtemp.EQ.0)xtemp=1
c>>>>>>>>>>>>>>>>>>>>>>>> END lines and Treitel
         xtemp=theta(i)
         xtemp=1
         jacscale(i)=xtemp
         DO ix=1,nx            
            DO jdep=1,ndep
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               grad(ix+index,i)=  grad(ix+index,i)*xtemp ! normalized
            ENDDO
         ENDDO
      ENDDO
 200  CONTINUE
c
      END

c**************************************************************
      SUBROUTINE movegrad
c-    This moves the gradient
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comgrad.h'
      REAL xhelp
      INTEGER ngrad2,ngrad3,i,i2,i3,i4,j,k,index
c---  set the model and CALL the modelling routine
c---- we normalize the gradient

      ngrad3=ngrad-nforw_eofvar
      ngrad2=ngrad+ nforw_eof
c---- reset
      DO k=1,ndep
         index=(k-1)*nx
         DO j=1,nx
            DO i=ngrad+1,ngrad2
               grad(j+index,i)=0 ! normalized
            ENDDO
         ENDDO
      ENDDO
c     
c---  we construct the EOF gradients after the first computed
c
      DO i=ngrad+1,ngrad2
         i3=par2lay(forw2eof(i-ngrad)) ! this tells which col to use
         DO i2=1,neofvar
            xhelp=eofcoef(i2,i3)
            IF (xhelp.EQ.0) GOTO 111
            i4=i2+ngrad3
            DO k=1,ndep
               index=(k-1)*nx
               DO j=1,nx
                  grad(j+index,i)=
     1            grad(j+index,i)+grad(j+index,i4)*xhelp  ! normalized
               ENDDO
            ENDDO
 111        CONTINUE
         ENDDO
      ENDDO
c
c---  put the physical paprameters on correct spot.
c
      DO i=ngrad3,1,-1
         i4=forw2parm(i)
         DO k=1,ndep
            index=(k-1)*nx
            DO j=1,nx
               grad(j+index,i4)=grad(j+index,i) ! normalized
            ENDDO
         ENDDO
      ENDDO

c
c---  Now we move them back
c
      DO i=ngrad+1,ngrad2
         i4=forw2eof(i-ngrad)
         DO k=1,ndep
            index=(k-1)*nx
            DO j=1,nx
               grad(j+index,i4)=grad(j+index,i) ! normalized
            ENDDO
         ENDDO
      ENDDO

      END

C***********************************************************
      SUBROUTINE map2forw(index,nparml,nparmu)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comgrad.h'
      INCLUDE 'comoas.h'
      INTEGER maplpflag
      COMMON /map2forwfl/maplpflag
      INTEGER index,iforw,iforwe,iparm,ieof,ilay,il,j,ii
      INTEGER nparml,nparmu
      iforw=0
      iforwe=0
      ii=0
      maplpflag= maplpflag+1
c     
c---  For local method on all paprameters
c     
      IF (index.EQ.1) THEN 
         DO iparm=1,nparm
            IF (par2phy(iparm).NE.11) THEN
               iforw=iforw+1
               par2phy_forw(iforw)=par2phy(iparm) 
               par2lay_forw(iforw)=par2lay(iparm) 
               forw2parm(iforw)=iparm ! points to the EOF_coefs
            ELSE
               iforwe=iforwe+1
               forw2eof(iforwe)=iparm ! points to the EOF_coefs
            ENDIF  
            ii=ii+1
            forw2opt(ii)=iparm  ! points to the EOF_coefs
         ENDDO
         nforw_eof=iforwe
         nthet = iforw + iforwe
         DO ieof=1,neofvar
            iforw=iforw+1
            IF (iforw.GT.mpar) STOP 'map2forw: increase mpar'
            par2phy_forw(iforw)=par2phy_eof(ieof) 
            par2lay_forw(iforw)=par2lay_eof(ieof) 
         ENDDO
         nforw_eofvar=neofvar
         ngrad=iforw
      ELSE
c     
c---  ONLY a local method on a fraction of the parameters
c     
         DO iparm=nparml,nparmu
            IF (par2phy(iparm).NE.11) THEN
               iforw=iforw+1
               par2phy_forw(iforw)=par2phy(iparm) 
               par2lay_forw(iforw)=par2lay(iparm) 
               forw2parm(iforw)=iparm ! points to the NON-EOF_coefs
            ELSE
               iforwe=iforwe+1
               forw2eof(iforwe)=iparm ! points to the EOF_coefs
            ENDIF  
            ii=ii+1
            forw2opt(ii)=iparm  ! points to the EOF_coefs
         ENDDO
         nforw_eof=iforwe
         nthet = iforw + iforwe
         DO ieof=1,neofvar
            DO iparm=1,nforw_eof
               IF  (forw2eof(iparm).EQ.ieof) GOTO 50
            ENDDO
            GOTO 51
 50         CONTINUE
            iforw=iforw+1
            IF (iforw.GT.mpar) STOP 'map2forw: increase mpar'
            par2phy_forw(iforw)=par2phy_eof(ieof) 
            par2lay_forw(iforw)=par2lay_eof(ieof) 
 51         CONTINUE
         ENDDO
         nforw_eofvar=neofvar
         ngrad=iforw
      ENDIF            

c***  How many parameters in each lay
      DO il=1,nlay
         Nparmlay(il)=0
      ENDDO
      DO iparm=1,ngrad
         ilay=par2lay_forw(iparm)
         Nparmlay(ilay)=Nparmlay(ilay)+1
         lay2parm(ilay,Nparmlay(ilay))=iparm
      ENDDO
      IF (maplpflag.LT.5) THEN
         WRITE(prtfil,*)'  ngrad= ',ngrad
         WRITE(prtfil,*)'  nthet= ',nthet
         WRITE(prtfil,*)'  nforw_eof = ',nforw_eof
         WRITE(prtfil,*)'  nforw_eofvar = ',nforw_eofvar
         DO iforw=1,ngrad
            WRITE(prtfil,*)'grad ', par2phy_forw(iforw),
     &           par2lay_forw(iforw),
     &           forw2opt(iforw) 
         ENDDO
         IF ((ngrad+nforw_eof).GT.mpar) THEN
            WRITE(*,*)' mpar not large enough for grad'
            WRITE(*,*)'ngrad+nforw_eof,mpar=',ngrad+nforw_eof,mpar
            STOP
         ENDIF
         WRITE(prtfil,*)' nparmlay for each layer'
         DO il=1,nlay
            IF (Nparmlay(il).GE.1)
     &           WRITE(prtfil,*)' lay',il,' points',
     &           (lay2parm(il,j),j=1,Nparmlay(il))
         ENDDO    
         WRITE(prtfil,*)
      ENDIF
      
      END
