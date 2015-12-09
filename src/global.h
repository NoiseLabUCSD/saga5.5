      MODULE global
      IMPLICIT NONE
      INTEGER,   DIMENSION(:,:), ALLOCATABLE ::  allmodel,model
      REAL   ,   DIMENSION(:),   ALLOCATABLE ::  allfit
      INTEGER,   DIMENSION(:,:), ALLOCATABLE ::  parents
      COMPLEX,   DIMENSION(:),   ALLOCATABLE ::  weight, DATA, cov
      COMPLEX,   DIMENSION(:),   ALLOCATABLE ::  Aun
      SAVE
      END MODULE global

