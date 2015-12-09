      SUBROUTINE forwardmodel(iopt,mopt)
      INTEGER  mopt,i,iopt(mopt) 
      DO i=1,mopt
         iopt(i)=0.
      ENDDO
      iopt(1)=2
      iopt(30)=9                ! 9 is for ramgeo.
      iopt(12)=1   
      END
c     
c     Initialize the parameters, acoustic field, and matrices.
c     
      SUBROUTINE input
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE  'comramgeo.h'
      INCLUDE 'ramgeo/ram.h'
      CHARACTER*80 dumch,  dumch2
      EQUIVALENCE (dumch,opt)
      INTEGER i,id,ic,irange,itype,j,ierr,idep
      REAL*8   zdep,zval,zvalold
      

      ierr=0
      nfrq=1
c---  READ the GA-parameters
      CALL readinputstart
      
      READ(1,'(a)') dumch2       ! This line has no effect in snap
c     CALL ramoption(opt)
      WRITE(*,*)'reading nfrq..'
      READ(1,*) nfrq
 92   WRITE(*,'(a,i3,a,i4)') ' nfrq=', nfrq
      IF (ABS(nfrq).GT.mfreq) STOP ' nfreq > mfreq' ! test
      
      IF (nfrq.GT.0) THEN
         READ(1,*)(frq(I),i=1,nfrq)
         WRITE(*,'(a,10f10.3)')' frequencies:',(frq(I),i=1,nfrq)
      ELSE
         STOP 'ramgeo: not yet implemented'
      ENDIF
      
      READ(1,*) zs
      pi=4.0*ATAN(1.0)
      omega=2.0*pi*freq
      READ(1,*)rmin,rmax,drng
 102  WRITE(*,*)'min,max range, ',rmin,rmax
      WRITE(*,*)'Hydrophone range spacing',drng
      IF (rmin.GT.rmax) THEN
         STOP 'rmin>rmax'
      ENDIF
      READ(1,*) dr0
      WRITE(*,*)'PE model dr',dr0
      dr=dr0
      rarray=rmax-rmin
      dr2=drng
      ndr=1
      READ(1,*)zmax,dz,ndz
      READ(1,*)rdmin,rdmax
      IF (rdmin.GT.rdmax) THEN
         STOP 'rdmin>rdmax'
      ENDIF
      zr=rdmin
      zmplt=rdmax
      irdmin=INT(rdmin/dz+0.5)
      rdmin=dz*irdmin
      
      rdmax=dz* INT(rdmax/dz+0.5)
      ndep=INT((rdmax-rdmin)/(dz*ndz))+1
      WRITE(*,*)'rdmin,rdlow,ndep',rdmin,rdlow,ndep
c     equidistant receiver depths
      DO 920 idep=1,ndep
 920     rdep(idep)=(idep-1)*(dz*ndz) + rdmin

         READ(1,*)c0,np,ns,rs0
         nz=zmax/dz-0.5
         nzplt=ndep

c---- READ the bathymetry
         i=1
 1       READ(1,*)rb(i),zb(i)
         IF(rb(i).LT.0.0)go to 2
         i=i+1
         go to 1
 2       CONTINUE
         IF (i.GT.1) THEN
            IF (rb(i-1).GT.200*rmax) THEN
               i=i-1
            ELSE
               rb(i)=200.0*rmax
               zb(i)=zb(i-1)
            ENDIF
         ENDIF
         nbathy=i
c---- check dimensions
c         IF(nz+2.GT.mz)THEN
c            WRITE(*,*)'   Need to increase parameter mz to ',nz+2
c            STOP
c         END IF
c         IF(np.GT.mp)THEN
c            WRITE(*,*)'   Need to increase parameter mp to ',np
c            STOP
c         END IF
         IF(i.GT.mr)THEN
            WRITE(*,*)'   Need to increase parameter mr to ',i
            STOP
         END IF
         rangeprof(1)=0
         DO irange=1,mrange
            DO itype=1,4
               ic=1  
               zval=0
               DO id=1,mpardep
                  zvalold=zval
                  READ(1,*,err=20) zdep,zval
                  geodep(ic,itype,irange)=zdep
                  geopar(ic,itype,irange)=zval
                  igeodep(ic,itype,irange)= zdep/dz+1.5
                  
                  IF (zdep.EQ.zmax) THEN
                     igeodep(ic,itype,irange)= zdep/dz+1.5-1
                     geodep(ic,itype,irange)=zdep-dz
                  ENDIF

                  IF ((id.EQ.1).AND. (igeodep(id,itype,irange).NE.1)) 
     &                 THEN
                     igeodep(ic+1,itype,irange)= 
     &                    igeodep(ic,itype,irange)
                     geopar(ic+1,itype,irange)=
     &                    geopar(ic,itype,irange)
                     igeodep(ic,itype,irange)=1
                     ic=ic+1
                  ENDIF


                  IF (geodep(ic,itype,irange).LT.0) THEN
                     geodep(ic,itype,irange)=(nz+1)*dz
                     geopar(ic,itype,irange)=zvalold
                     igeodep(ic,itype,irange)=nz+2
                     GOTO 22
                  ELSEIF (igeodep(ic,itype,irange).LT.(nz+2)) THEN
                     IF (ic.GE.2) THEN
                        IF (igeodep(ic,itype,irange).EQ.
     &                       igeodep(ic-1,itype,irange)) then
                           igeodep(ic,itype,irange)=
     &                          igeodep(ic-1,itype,irange)+1
                        ENDIF
                     ENDIF
                     ic=ic+1
                  ENDIF
               ENDDO
 22            ndepprof(irange,itype)=ic
               IF (ic.GT.Mpardep ) THEN
                  WRITE(*,*) 'Need to increase parameter',
     #                 ' Mpardep to ',(ic+1),' or more'
                  STOP
               ENDIF
c     234567890123456789012345678901234567890123456789012345678901234567890
            ENDDO
            READ(1,*,err=20)rangeprof(irange+1)
            IF (rangeprof(irange+1).LT.0) GOTO 20
         ENDDO
 20      Nprof=irange
         rangeprof(Nprof+1)=2*rmax
         IF ((Nprof+1).GT.Mrange) THEN
            WRITE(*,*) 'Need to increase parameter Mrange to ',
     &           (Nprof+1)
            STOP
         ENDIF

         CALL writeenvmodel
         IF (drng.GT.0) THEN
            nx=INT((rmax-rmin)/drng +1)
            nx=MAX(1,nx)
         ELSE
            nx=1
         ENDIF
         WRITE(*,*)'nx=',nx
         DO i=1,nx
            xranges(i)=rmin+drng*(i-1)
         ENDDO


         IF (nx.LE.0) THEN
            WRITE(*,*)
            WRITE(*,*)' rmin, rmax,dr,ndr=',rmin, rmax,dr,ndr
            WRITE(*,*) ' The field must be output',
     &           'for more than one range point'
            WRITE(*,*)' Change range parameters'
            STOP
         ENDIF

         ncurv=nfrq*ndep        ! iopt(5)=1,3,7
         
c*****************Copied from SUBROUTINE "read_ran_in" **************
         IF ((iopt(5).EQ.4) .OR. (iopt(5).EQ.6) 
     &        .OR. (iopt(13).EQ.1)) THEN
            ncurv=nfrq*nx
         ENDIF
         IF (iopt(13).GT.0) THEN ! for covariance matrix
            nbart=nfrq*nx
            WRITE(prtfil,*)'number of bartlett estimators',nbart
            IF (nbart.GT.mfreq) THEN
               WRITE(*,*)'nbart > mfreq'
               ierr=1
            ENDIF
            ncov_siz=nbart*ndep*ndep
            WRITE(*,*)' Elements in cov-matrix:',ncov_siz
            IF (ncov_siz.GT.mobs) THEN
               WRITE(*,*)'ncov_siz > mobs,  mobs=',mobs
               ierr=1
            ENDIF
         ELSEIF (iopt(5).EQ.5) THEN
            ncurv=nx*ndep       ! broadband matched filter
         ENDIF

         WRITE(prtfil,*)' ncurv=',ncurv,'nbart=',nbart
         IF ((nx.EQ.1).AND.(ndep.EQ.1).AND.(nfrq.EQ.1)) THEN
            WRITE(*,*) ' nx, ndep and nfrq must not all be 1'
            ierr=1
         ENDIF
         IF ((iopt(5).EQ.7).AND. (nx.EQ.1)) THEN
            WRITE(*,*)' *****   Error   ******'
            WRITE(*,*)' Only one phone specified in range '
            WRITE(*,*)' at least two phones required ',
     &           'for this inversion'
            ierr=1
         ENDIF
         IF ((iopt(5).EQ.4).AND. (ndep.EQ.1)) THEN
            WRITE(*,*)' *****   Error   ******'
            WRITE(*,*)' Only one phone specified in depth'
            WRITE(*,*)' at least two phones required for ',
     &           'this inversion'
            ierr=1
         ENDIF
         IF ((iopt(5).EQ.5).AND. (nfrq.EQ.1)) THEN
            WRITE(*,*)' *****   Error   ******'
            WRITE(*,*)' Only one frequency specified'
            WRITE(*,*)'at least two frequencies required ',
     &           'for this inversion'
            ierr=1
         ENDIF
         IF (nx*ndep.GT.mobs) THEN
            WRITE(*,*) ' nx*ndep >mobs, nx,ndep,mobs',nx,ndep,mobs
            ierr=1
         ENDIF
         IF (nx.GT.mx) THEN
            WRITE(*,*) 'nx, mx=',nx,mx
            STOP 'nx is greater than mx'
         ENDIF

c---  is this what we optimize over


c---  should the  EOF be READ ? 

         IF (iopt(17).EQ.1) THEN
            CALL eofinit
         ENDIF

c***  create text strings for physical parameters
c     123456789012345678901234567890
         phystxt(1) = 'Ocean sound speed (m/s)$'
         phystxt(2) = 'Bottom sound speed (m/s)$'
         phystxt(3) = 'Bottom density (g/cm^3)$'
         phystxt(4) = 'Bottom attenuation (dB)$'
         phystxt(5) = 'Depth (m)$'
         phystxt(6) = 'Depth (m)$'
         phystxt(7) = 'Depth (m)$'
         phystxt(8) = 'Depth (m)$'
         phystxt(9) = 'Range (km)$'
         phystxt(11)= 'Shape coefficient$'
         phystxt(12)= 'Bathymetry (m)$'
         phystxt(13)= 'Bathymetry range (m)$'
c***  

         DO 8 j=1,mphys
            phystxt2(j)='                                             '
 6          DO 7 I=40,1,-1
               IF(phystxt(j)(I:I).NE.'$') GO TO 7
               phystxt2(j)(1:I-1)=phystxt(j)(1:I-1)
c     WRITE(*,*)phystxt2(j)
               GO TO 8
 7          CONTINUE
 8       CONTINUE
         
c---- READ the optimization PARAMETER
c         WRITE(*,*) 'hi peter'
         CALL readoptparm

         ichangedepth=0
         DO i=1,nparm
            IF (par2phy(i).EQ.5 .OR.par2phy(i).EQ.6 .OR. 
     &           par2phy(i).EQ.7 .OR.par2phy(i).EQ.10 .OR. 
     &           par2phy(i).EQ.11) THEN
               ichangedepth=1
            ENDIF
            IF (fmin(i).GT.fmax(i))THEN
               WRITE(*,*)' *** fmin > fmax for parm',i
               ierr=1
            ENDIF
            IF (par2phy(i).LT.1 .OR. par2phy(i).GT.19 
     &           )THEN
               WRITE(*,*)' *** par2phy not correct, parm',i
               ierr=1
            ELSE
               IF (par3(i).GT.Nprof) THEN
                  WRITE(*,*)' *** Optimzation variable #:',i
                  WRITE(*,*)' This points past the last profile'
               ENDIF
            ENDIF
            IF (par2phy(i).EQ.11 )THEN
               IF (par2lay(i).GT.neof) THEN
                  WRITE(*,*)' *** Optimzation variable #:',i
                  WRITE(*,*)' *** The shapecoffient number is not',
     &                 ' defined'
                  WRITE(*,*)' par2lay(i)', par2lay(i)
                  ierr=1
               ENDIF
            ENDIF
            IF (par2phy(i).EQ.5 )THEN
               IF (par2lay(i).GT.ndepprof(par3(i),1)) THEN
                  WRITE(*,*)' *** Optimization variable #:',i
                  WRITE(*,*)' The depth point is below last point ', 
     &                 ndepprof(par3(i),1),
     &                 'depth points, the pointer is too deep'
                  WRITE(*,*)' par2lay(i)', par2lay(i)
                  ierr=1
               ENDIF
            ENDIF
            IF (par2phy(i).EQ.6 )THEN
               IF (par2lay(i).GT.ndepprof(par3(i),2)) THEN
                  WRITE(*,*)' *** Optimization variable #:',i
                  WRITE(*,*)' The depth point is below last point ', 
     &                 ndepprof(par3(i),2),
     &                 'depth points, the pointer is too deep'
                  WRITE(*,*)' par2lay(i)', par2lay(i)
                  ierr=1
               ENDIF
            ENDIF
            IF (par2phy(i).EQ.7 )THEN
               IF (par2lay(i).GT.ndepprof(par3(i),3)) THEN
                  WRITE(*,*)' *** Optimization variable #:',i
                  WRITE(*,*)' The depth point is below last point ', 
     &                 ndepprof(par3(i),3),
     &                 ', the pointer is too deep'
                  WRITE(*,*)' par2lay(i)', par2lay(i)
                  ierr=1
               ENDIF
            ENDIF
            IF (par2phy(i).EQ.10 )THEN
               IF (par2lay(i).GT.ndepprof(par3(i),4)) THEN
                  WRITE(*,*)' *** Optimization variable #:',i
                  WRITE(*,*)' The depth point is below last point ', 
     &                 ndepprof(par3(i),4),
     &                 ', the pointer is too deep'
                  WRITE(*,*)' par2lay(i)', par2lay(i)
                  ierr=1
               ENDIF
            ENDIF
c     IF (par2phy(i).EQ.1 )THEN
c     IF (fmin(i).LT.z0(nd0-1)) THEN
c     ENDIF
c     ENDIF
c     IF (par2phy(i).EQ.17)THEN
c     ENDIF

         ENDDO

c***  errors ?
         IF (ierr.EQ.1)STOP 'stopped in input'
c***  CLOSE input unit
         CLOSE(1)                                    
         END    

      
      SUBROUTINE forwinit
c     
c     This version of ram handles multiple sediment layers that 
c     parallel the bathymetry.
c     
c     ******************************************************************
c     ***** Range-dependent Acoustic Model, Version 1.4g, 05-Jul-00 ****
c     ******************************************************************
c     
c     This code was developed by Michael D. Collins at the Naval
c     Research Laboratory in Washington, DC. It solves range-dependent 
c     ocean acoustics problems WITH the split-step Pade algorithm 
c     [M. D. Collins, J. Acoust. Soc. Am. 93, 1736-1742 (1993)]. A 
c     user's guide and updates of the code are available via anonymous 
c     ftp from ram.nrl.navy.mil. 
c     
c     Version 1.5g contains a correction to a bug in the dimension of
c     quantities passed to subroutines fndrt and guerre that Laurie
c     Fialkowski noticed. 
c     
c     Version 1.4g contains a correction to a minor bug in subroutine
c     guerre that Dave King noticed (amp1 and amp2 were declared
c     twice) and a few other minor improvements. 
c     
c     Version 1.3g contains a new root-finding subroutine.
c     
c     Version 1.2g contains a minor modification. The output to tl.grid
c     is no longer zeroed out along the ocean bottom. This was done in
c     previous versions so that the ocean bottom would be highlighted
c     in graphical displays. The graphics codes ramclr, ramctr, and
c     ramcc read in the bathymetry from ram.in and plot the ocean
c     bottom directly. 
c     
c     Version 1.1g contains two improvements:
c     
c     (1) An improved self starter. Stability is improved by using the 
c     factor (1-X)**2 instead of (1+X)**2 to smooth the delta function. 
c     The factor (1+X)**2 is nearly singular for some problems involving
c     deep water and/or weak attenuation. Numerical problems associated 
c     with this singularity were detected by Eddie Scheer of Woods Hole 
c     Oceanographic Institute. 
c     
c     (2) Elimination of underflow problems. A very small number is 
c     added to the solution in subroutine solve to prevent underflow,
c     which can adversely affect run time on some computers. This
c     improvement was suggested by Ed McDonald of the SACLANT Undersea
c     Research Centre. 
c     
c     
c     mr=bathymetry points, mz=depth grid, mp=pade terms.
c     
c     parameter (mr=100,mz=20000,mp=10)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      include 'comramgeo.h'
      include 'ramgeo/ram.h'
      real cw(mz),cb(mz),rhob(mz),attn(mz),alpw(mz),
     >     alpb(mz),f1(mz),f2(mz),f3(mz),ksqw(mz),tlg(mz)
      complex ci,ksq(mz),ksqb(mz),u(mz),
     >     v(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),s1(mz,mp),
     >     s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)

c      real, dimension(:),allocatable :: cw,cb,rhob,attn,alpw,
c     >     alpb,f1,f2,f3,ksqw,tlg
c      complex ci
c      complex, dimension(:),allocatable :: ksq,ksqb,u,v,pd1,pd2
c      complex, dimension(:,:),allocatable :: r1,r2,r3,s1,s2,s3
      real k0,r,eps,dir,eta,rp,rs
      integer ir,lz,mdr,ib,ifrq
      integer iz,index,id,idep,irange
      complex phasefactor 
      flagpu=-3
c      write(*,*)nz,nprof
c      mz=nz+4
c      mp=nprof+1
c      allocate(cw(mz))
c      allocate(cb(mz))
c      allocate(rhob(mz))
c      allocate(attn(mz))
c      allocate(alpw(mz))
c      allocate(alpb(mz))
c      allocate(f1(mz))
c      allocate(f2(mz))
c      allocate(f3(mz))
c      allocate(ksqw(mz))
c      allocate(tlg(mz))
c      allocate(ksq(mz))
c      allocate(ksqb(mz))
c      allocate(u(mz))
c      allocate(v(mz))
c      allocate(pd1(mp))
c      allocate(pd2(mp))
c      allocate(r1(mz,mp))
c      allocate(r2(mz,mp))
c      allocate(r3(mz,mp))
c      allocate(s1(mz,mp))
c      allocate(s2(mz,mp))
c      allocate(s3(mz,mp))

      entry forw2

c      write(*,*)'updated'
c     
c---- If writing out trf function
c     
      if (iWriteTrf.ge.100 ) then
         rmin=0
         rmax=5000
         nzplt=2
c         nxplt=2
c     
         rangeprof(Nprof+1)=2*rmax    
c     

         rb(nbathy)=200*rmax
         zmplt=180
         nfrq=1
         frq(1)=250
c     frq(1)=300
         write(*,*)' Running ramgeo For optimised env at',
     &        frq(1),' Hz'
c     Fmin, Fmax, DF, nfrq: ', 
c     .        Fmindum, Fmaxdum, DF_temp,nfrq
         dr=5
         dr0=5
         dr2=dr0
         ndr=5
         ndz=20    
      ELSE
c     flagpu=1
      ENDIF

c     ENDIF
      lwascomp=1
c     WRITE(*,*) 'has entered snap'
      flagpu=1 +flagpu

c     debugging
c     
      IF (iopt(6).EQ.1) flagpu=-10 
c     
c---- USE of EOF
c     

      IF (iopt(17).EQ.1) THEN
         CALL eofval
         IF (lwascomp.EQ.-1) THEN
            RETURN
         ENDIF
      ENDIF
c     make sure last bathymetry point is correct
      zb(nbathy)=   zb(nbathy-1)
c     
c---- for modifying depths in sediment
c     
      CALL changedepth

      IF (iWriteTrf.LT.1 ) THEN
         IF (nx.EQ.1) THEN
c     ndr=INT(rmax/dr+0.5)
            ndr=1
            rmax=rmin
         ELSE
            rmax=rmin+rarray
            dr=rmin/INT(rmin/dr0)
c     WRITE(*,*)' new dr0',dr,rmax,rmin,rarray
         ENDIF
      ENDIF

      
c     WRITE(*,*)' new dr0',dr,rmax,rmin,rarray,ndr,nx

      IF (flagpu.LT.0) CALL writeenvmodel
      

      DO ifrq=1,nfrq
         omega=2.*pi*frq(ifrq)
c       WRITE(*,*)'omega,frq', omega,frq(ifrq),ndz,r,ndr
         irange=0
         dr=dr0

         CALL setup(mr,mz,nz,mp,np,ns,mdr,ndr,ndz,iz,nzplt,lz,
     >        ib,ir,dir,dr,
     >        dz,pi,eta,eps,omega,rmax,c0,k0,ci,r,rp,rs,rb,zb,
     >        cw,cb,rhob,
     >        attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,
     >        r3,s1,s2,s3,
     >        pd1,pd2,tlg)


         IF (iWriteTrf.gE.100 ) THEN
c            ndr=1
            mdr=0
            WRITE(16,*)dr*ndr,dz*ndz,INT(rmax/(ndr*dr)+0.5)-1,lz
         ENDIF
c     
c     March the acoustic field out in range.
c     
c     WRITE(*,*)' before  updat',r,rmin,rmax,mdr,dr,dr0,dr2
 1       CALL updat(mr,mz,nz,mp,np,iz,ib,dr,dz,eta,omega,
     >        rmax,c0,k0,ci,r,
     >        rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb,ksq,
     >        ksqw,ksqb,f1,f2,f3,
     >        r1,r2,r3,s1,s2,s3,pd1,pd2)

c      WRITE(*,*)' after  updat',r,rmin,rmax,mdr,dr,dr0,dr2
         CALL solvegeo(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
         r=r+dr
c      WRITE(*,*)' range',r,rmin,rmax,mdr,dr,dr0,dr2
c     DO i=1,nz+2 
c     WRITE(91,*)f1(i),f2(i),f3(i)
c     ENDDO
c     WRITE(99,*)r,-20.0*alog10(cabs(f3(ir)*u(ir)))+10.0*alog10(r)
c     1    ,f3(ir),u(ir),ir,dir
c     
         IF (iWriteTrf.gE.100) THEN
c       WRITE(*,*)'am i here, r=',r,mdr,ndr
            CALL outpt(mz,mdr,ndr,ndz,iz,nzplt,lz,ir,dir,eps,
     &           r,f3,u,tlg)
c       WRITE(*,*)'after outpt, r=',r
         ELSEIF (r.GE.rmin) THEN
c     WRITE(*,*)'why do I not come here',ndep,ndr,mdr
            mdr=mdr+1
            IF (mdr.EQ.ndr) THEN
               mdr=0
               phasefactor=cexp(CMPLX(0.,omega/c0*r))/SQRT(r)

               irange=irange+1      
               DO id=1,ndep     ! irin
                  idep=irdmin+(id-1)*ndz
                  index=(id +((ifrq-1))*ndep-1)*nx
c            WRITE(prtfil,*)'resp',resp(irange+index),irange,index,idep,ndz
                  resp(irange+index)=u(idep)*f3(idep)*phasefactor
               ENDDO            ! depth
            ENDIF
         ENDIF      
c     
         IF(r.LT.rmax)go to 1
c     WRITE(95,*)' range',r
c     DO i=1,nz+2 
c     WRITE(95,*)u(i)
c     ENDDO
c      WRITE(*,*)' hello a'
         IF ((iopt(6).EQ.1) .AND. (iWriteTrf.lt.1)) THEN
            DO id=1,ndep        ! irin
               idep=irdmin+(id-1)*ndz
               index=(id +((ifrq-1))*ndep-1)*nx
c               write(*,*)'id,idep,index,ndep,nx',id,idep,index,ndep,nx
c               write(*,*)'irange,irange+index',irange,irange+index
               WRITE(prtfil,*)'resp ',resp(irange+index),
     &              irange,index,idep
               WRITE(prtfil,*)'data ',DATA(irange+index)
            ENDDO               ! depth
         ENDIF
      ENDDO                     !frequency
      RETURN
      END
c     
c---- 
c     
      SUBROUTINE writeenvmodel
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE  'comramgeo.h'
      INCLUDE 'ramgeo/ram.h'
      INTEGER i,id,irange,itype


      WRITE(prtfil,*)'Frequency:',(frq(i),i=1,nfrq)
      WRITE(prtfil,*)'zs=', zs
      pi=4.0*ATAN(1.0)
      omega=2.0*pi*freq
      WRITE(prtfil,*)'rmin,rmax,dr,ndr',rmin,rmax,dr,ndr
      WRITE(prtfil,*)'zmax,dz,ndz',zmax,dz,ndz
      WRITE(prtfil,*)'rdmin,rdmax,ndz',rdmin,rdmax,ndz
      WRITE(prtfil,*)'  c0,np,ns,rs',c0,np,ns,rs0
      WRITE(prtfil,*)'bathymetry'
      DO i=1,nbathy
         WRITE(prtfil,*)rb(i),zb(i)
      ENDDO
      DO irange=1,Nprof
         WRITE(prtfil,*)'Profile at range:',rangeprof(irange)
         DO itype=1,4
            IF (itype.EQ.1) THEN
               WRITE(prtfil,*)' Ocean Sound Speed'
            ELSEIF (itype.EQ.2) THEN
               WRITE(prtfil,*)' Bottom Sound Speed'
            ELSEIF (itype.EQ.3) THEN
               WRITE(prtfil,*)' Bottom Density'
            ELSEIF (itype.EQ.4) THEN
               WRITE(prtfil,*)' Bottom Attenuation'
            ENDIF
            WRITE(prtfil,*)'Number of points:',ndepprof(irange,itype)

            DO id=1,ndepprof(irange,itype)
c     WRITE(*,*)'id,itype,irange',id,itype,irange
               WRITE(prtfil,*) 
     #              geodep(id,itype,irange),geopar(id,itype,irange),
     #              igeodep(id,itype,irange)
            ENDDO
         ENDDO
      ENDDO
      END
c     
c---- 
c     
      SUBROUTINE changedepth
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE  'comramgeo.h'
      INCLUDE 'ramgeo/ram.h'
      INTEGER id,irange,itype,ic
      REAL zval,zvalold,zdep


      DO irange=1,Nprof
         DO itype=1,4
            ic=1  
            zval=0
            DO id=1,mpardep
               zvalold=zval
               zval= geopar(ic,itype,irange)
               zdep= geodep(ic,itype,irange)
               igeodep(ic,itype,irange)= zdep/dz+1.5

               IF ((id.EQ.1).AND. (igeodep(id,itype,irange).NE.1)) THEN
                  igeodep(ic+1,itype,irange)=igeodep(ic,itype,irange)
                  geopar(ic+1,itype,irange)=geopar(ic,itype,irange)
                  igeodep(ic,itype,irange)=1
                  ic=ic+1
               ENDIF

               IF (geodep(ic,itype,irange).EQ.((nz+1)*dz)) THEN
c     geodep(ic,itype,irange)=(nz+1)*dz
                  geopar(ic,itype,irange)=zvalold
c     igeodep(ic,itype,irange)=nz+2
                  GOTO 22
               ELSEIF (igeodep(ic,itype,irange).LT.(nz+2)) THEN
                  IF (ic.GE.2) THEN
                     IF (igeodep(ic,itype,irange).EQ.
     #                    igeodep(ic-1,itype,irange)) then
                        igeodep(ic,itype,irange)=
     &                       igeodep(ic-1,itype,irange)+1
                     ENDIF
                  ENDIF
                  ic=ic+1
               ENDIF
            ENDDO
 22         IF (ic.GT.Mpardep ) THEN
               CALL  writeenvmodel
               WRITE(*,*) 'from changedepth:'
               WRITE(*,*) 'Need to increase parameter',
     #              ' Mpardep to more than',(ic+1)
               WRITE(*,*)' ..the environment model is written to *.out'
               STOP
            ENDIF
c     234567890123456789012345678901234567890123456789012345678901234567890
         ENDDO
      ENDDO
      END


c*********************************************
      SUBROUTINE writeTL()
      USE global
      IMPLICIT NONE
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comramgeo.h'

      INTEGER ifreq,i,id,index,mx1,nf1
c     pg      INTEGER*4 jj0,jsrc,jrec,j,mx1
      REAL*4 rspace,freqs,dt1
      LOGICAL bintrf
      INTEGER luttrf
      data luttrf/16/
      CHARACTER*30 outroot      !not used
      bintrf=.FALSE.

c     bintrf=.TRUE.

c     
c     rspace = range increment for receiver position
c     pln      rkm(1)=rkm(1)*1.0e-3
      IF (nx.EQ.1) THEN
         rspace=1
      ELSE
         rspace=dr
      END IF
      mx1=1 !nf1+nfrq-1
      dt1=1.E0 !/fsbb
      freqs=1  !(fmaxdum-fmindum)/2.+fmindum
c     sd=zsr(mzsrc(1))
      nf1=1
      if (iforwt==1) then
      write(*,*)'Calling trfhead, np,ndep=',nx,ndep
      CALL trfhead(outroot,title,rdep(1),rdep(ndep),
     &     dr0,rspace,nfrq,nf1,mx1,dt1,freqs,sd,
     &     bintrf,ndep,1,nx)
c      CALL trfhead(outroot,title,rdep(1),rdep(ndep),
c    &     rng(1),rspace,nfftbb,nf1,mx1,dt1,freqs,sd,
c    &     bintrf,ndep,1,np)
      endif
c                  write(*,*)'writing trf..'
                  WRITE(luttrf,*)
                  WRITE(luttrf,'(1000000f8.3)')
     &    ((((20*log10(abs(resp(i+(id +((ifreq-1))*ndep-1)*nx))))
     &        ,id=1,ndep),i=1,nx), ifreq=1,nfrq)
      RETURN
      END  !wrtetrf

ccc
 
      SUBROUTINE writecomplex
      USE global
      IMPLICIT NONE
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comramgeo.h'

      INTEGER ifreq,i,id,index,mx1,nf1
c     pg      INTEGER*4 jj0,jsrc,jrec,j,mx1
      REAL*4 rspace,freqs,dt1
      LOGICAL bintrf
      INTEGER luttrf
      data luttrf/16/
      CHARACTER*30 outroot      !not used
      bintrf=.FALSE.

      IF (nx.EQ.1) THEN
         rspace=1
      ELSE
         rspace=dr
      END IF
      mx1=1 !nf1+nfrq-1
      dt1=1.E0 !/fsbb
      freqs=1  !(fmaxdum-fmindum)/2.+fmindum
c     sd=zsr(mzsrc(1))
      nf1=1
      if (iforwt==1) then
      write(*,*)'Calling trfhead, np,ndep=',nx,ndep
      CALL trfhead(outroot,title,rdep(1),rdep(ndep),
     &     dr0,rspace,nfrq,nf1,mx1,dt1,freqs,sd,
     &     bintrf,ndep,1,nx)
c      CALL trfhead(outroot,title,rdep(1),rdep(ndep),
c    &     rng(1),rspace,nfftbb,nf1,mx1,dt1,freqs,sd,
c    &     bintrf,ndep,1,np)
      endif

                   WRITE(luttrf,*)
                  WRITE(luttrf,'(1000000f15.8)')
     &    ((( REAL(resp(i+(id +((ifreq-1))*ndep-1)*np)), 
     &        imag(resp(i+(id +((ifreq-1))*ndep-1)*np))
     &        ,id=1,ndep),i=1,np), ifreq=1,nfrq)
                  WRITE(34)
     &    ((( REAL(resp(i+(id +((ifreq-1))*ndep-1)*np)), 
     &        imag(resp(i+(id +((ifreq-1))*ndep-1)*np))
     &        ,id=1,ndep),i=1,np), ifreq=1,nfrq)
      RETURN
      END  !writecomplex
ccc




