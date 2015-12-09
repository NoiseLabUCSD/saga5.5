      subroutine unceir(my,ny,theta,ntheta)

      
      INCLUDE 'comopt.h'
      INCLUDE 'comgrad.h'
      INCLUDE 'comforw.h'
c-----dimensions of data
      INTEGER my,ny,ntheta
c-----parameter vector
      REAL theta(mpar)
c-----pseudo-Hessian
      REAL zx(mx,mpar)   ! matrix to find SVD of 
      REAL ssvd(mpar+1),esvd(mpar),usvd(mx,mpar ),vsvd(mpar,mpar)
      REAL work2(mx)             
      INTEGER nsvd   ! number of observations
c-----old parameter vector	$(forde) snapinit$(for) $(object)snapinit$(obj)

      REAL V    ! current and old value of square error
c-----counters
      INTEGER j,k,n,m,index
      LOGICAL linfo  ! info control
c-----information returned
      INTEGER info
c-----check dimensions
      nsvd=ny*my
      if(nsvd.gt.mx) stop 'gnmin: Nsvd > mx; increase mx'
      IF (nthet.gt.mpar) stop ' gnmin: increase mpar  '
c-----calculate gradient
      CALL jacobi(theta)
      IF (iopt(4).ne.2) then
        linfo=.false.
      ELSE
       linfo=.true.
      ENDIF

       V=0
        DO n=1,ny
          index=(n-1)*my
          DO m=1,my
            resp(m+index)=(resp(m+index))-(data(m+index))
            V=V+(resp(m+index))**2
         ENDDO 
        ENDDO 
       IF (linfo) WRITE(*,*)'the start fitness is',V

c
c-----  computation of SVD
c
        DO j=1,ntheta
          k=0
          DO n=1,ny
            index=(n-1)*my
            DO m=1,my
              k=k+1
              zx(k,j)=grad(m+index,j)             
            ENDDO 
          ENDDO 
        ENDDO 
    
        call ssvdc(zx,mx,nsvd,ntheta,ssvd,esvd,usvd,mx,vsvd,mpar,
     &           work2,21,info)
      if (info.ne.0) then
        write(*,*)'********* info from zsvdc ************',info
      endif
c
c***** plot the singular vectors
c
      call plotsvd(vsvd,mpar,ntheta)

c      do j=1,ntheta
c        WRITE(prtfil,'(a,3e14.6)')'svd',ssvd(j),(vsvd(i,j),i=1,2)
c      enddo  
      
      return
      end
c******************************************************
