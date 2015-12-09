      MODULE fgs_com
          
      IMPLICIT NONE

      REAL, DIMENSION(:), ALLOCATABLE :: rdarr,rrarr
!     pg      REAL, DIMENSION(:), ALLOCATABLE :: frq,sdev,rdarr,rrarr
!     pg      REAL, DIMENSION(:), ALLOCATABLE :: mtrue,maxlim,minlim                      
!     pg       INTEGER, DIMENSION(:), ALLOCATABLE :: midx1,midx2
      INTEGER, PARAMETER :: v_arr=1,h_arr=2,nZL=1,ntZL=1

      INTEGER nran,ndim,nmod,arraytype
      INTEGER nhyd,ocount,i,j,k               
     		
      REAL  dtilt,rrva,rrhv,valength	
      LOGICAL optc,optp,tilt                    
      END MODULE fgs_com
