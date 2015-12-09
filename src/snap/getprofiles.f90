SUBROUTINE getprofiles 
!PROGRAM getprofiles 
  IMPLICIT NONE
  INTEGER, PARAMETER           :: MAX_SIZE = 500
  ! maximum vector size for z_sed
  DOUBLE PRECISION, DIMENSION(1:MAX_SIZE) :: z1, c1p, c1s, r1,&
       alp2betap, alp2betas
  DOUBLE PRECISION             :: cw_bot, phim, max_dep, d1, d2
  INTEGER                      :: ndep, i
!----------------------
! OPEN OUTPUT FILE
!----------------------
  OPEN(10, FILE='temp.m', status='unknown')
  OPEN(11, FILE='temp.dat', status='old')

  READ(11,*)  phim, cw_bot, ndep, max_dep

  z1(1) = 0.
  d1 = -1.301
  d2 = LOG10(max_dep)
  DO i = 2, ndep-1 
     z1(i) = 10 ** (d1+ (i-2.) * (d2 - d1) / (ndep-1.) )
  ENDDO
  z1(ndep) = max_dep

  CALL PHI2PARM(cw_bot, ndep, z1, phim, c1p, c1s, r1, alp2betap, alp2betas) 

  WRITE(10,29)
  DO i = 1, ndep
     WRITE(10,19) z1(i), c1p(i), c1s(i), r1(i), alp2betap(i), alp2betas(i)
  END DO
  WRITE(10,39) 

19 FORMAT( 6E15.5,';')
29 FORMAT('profs = [')
39 FORMAT('];')
  CLOSE(10)
  CLOSE(11)
!CONTAINS
END SUBROUTINE getprofiles 
!END PROGRAM getprofiles 
! ----------------------------------------------------------
! DOUBLE PRECISION FUNCTION  alpp_sur(alpp_0, phim):
! ----------------------------------------------------------
DOUBLE PRECISION  FUNCTION ALPP_SUR(phim)
  IMPLICIT  NONE 
  DOUBLE PRECISION, INTENT(IN) :: phim 
  DOUBLE PRECISION             :: alpp_0

  IF ((phim > 0) .AND. (phim < 2.5)) THEN
     alpp_0 = 0.23 + 0.0268 * phim
  ELSE IF ((phim >= 2.5 ) .AND. (phim < 4.12)) THEN
     alpp_0 = -0.1516607 + 0.1794643 * phim
  ELSE IF ((phim >= 4.12) .AND. (phim < 4.82)) THEN
     alpp_0 =  2.556396  - 0.4778311 * phim
  ELSE IF ((phim >= 4.82) .AND. (phim < 5.35)) THEN
     alpp_0 =  1.206895  - 0.1978517 * phim
  ELSE IF ((phim >= 5.35) .AND. (phim < 6)) THEN
     alpp_0 =  0.6295265 - 0.08993236*phim
  ELSE IF ((phim >= 6) .AND. (phim < 7.1)) THEN
     alpp_0 =  0.2974684 - 0.03458935 * phim
  ELSE IF (phim >= 7.1) THEN
     alpp_0 =  0.1299449 - 0.01099449 * phim
  END IF
  alpp_sur = alpp_0

END FUNCTION alpp_sur

! ----------------------------------------------------------
! SUBROUTINE  PHI2PARM():
! ----------------------------------------------------------
SUBROUTINE PHI2PARM(cw_bot, ndep, z1, phim, c1p, c1s, r1, alp2betap, alp2betas) 
  IMPLICIT  NONE
  INTEGER, INTENT(IN)      :: ndep
  DOUBLE PRECISION, INTENT(IN)         :: cw_bot, phim
  DOUBLE PRECISION, DIMENSION(1:ndep), INTENT(IN) :: z1
  DOUBLE PRECISION, DIMENSION(1:ndep), INTENT(OUT):: c1p, c1s, r1, &
       alp2betap, alp2betas
  DOUBLE PRECISION         :: alpp_0, SSR, fac1 = 0.0 , fac2 = 0.0
  DOUBLE PRECISION, DIMENSION(1:ndep)  :: c1p_1, c1s_1, r1_1, alp1p_1, alp1s_1
  DOUBLE PRECISION, DIMENSION(1:ndep)  :: c1p_2, c1s_2, r1_2, alp1p_2, alp1s_2
  DOUBLE PRECISION, DIMENSION(1:ndep)  :: alp1p, alp1s
  ! these two variables are for converting alpha_f to alpha_wavelength 
  DOUBLE PRECISION         :: alpp_sur
  INTEGER                  :: i                            ! index

    c1p_1 = 0.0
    c1p_2 = 0.0
    c1s_1 = 0.0
    c1s_2 = 0.0
    r1_1  = 0.0
    r1_2  = 0.0
    alp1p_1 = 0.0 
    alp1p_2 = 0.0 
    alp1s_1 = 0.0 
    alp1s_2 = 0.0 
  SSR = 1.18 - 0.034 * phim + 0.0013*phim ** 2.   ! sound speed ratio
  alpp_0 = alpp_sur(phim)   ! find surficial compressional attenuation

! determine which subroutine is used
  IF ((phim < 0) .AND. (phim >= 10)) THEN
     WRITE(*,*) 'phim=',phim
     STOP 'phim too small or too large'
     RETURN
  ELSE IF (phim <= 3.25) THEN   ! sand: coarse-grained sediment (phim <= 3.25)
     CALL COARSE_GRAINED(phim, cw_bot, SSR, alpp_0, z1, ndep, &
          c1p, c1s, r1, alp1p, alp1s)
  ELSE IF (phim >= 5.75) THEN   ! silt: fine-grained sediment (phim >= 5.75)
     CALL FINE_GRAINED(phim, cw_bot, SSR, alpp_0, z1, ndep, &
          c1p, c1s, r1, alp1p, alp1s)
  ELSE IF ((phim > 3.25) .AND. (phim < 5.75 )) THEN
     CALL FINE_GRAINED(phim, cw_bot, SSR, alpp_0, z1, ndep, &
          c1p_1, c1s_1, r1_1, alp1p_1, alp1s_1)
     CALL COARSE_GRAINED(phim, cw_bot, SSR, alpp_0, z1, ndep, &
          c1p_2, c1s_2, r1_2, alp1p_2, alp1s_2)
     fac1 = (phim - 3.25)/(5.75-3.25)
     fac2 = 1-(phim - 3.25)/(5.75-3.25)
     c1p = c1p_1*fac1 + c1p_2*fac2
     c1s = c1s_1*fac1 + c1s_2*fac2
     r1  = r1_1*fac1 + r1_2*fac2
     alp1p = alp1p_1*fac1 + alp1p_2*fac2
     alp1s = alp1s_1*fac1 + alp1s_2*fac2
  END IF

! Convert the units of attenuations to dB/wavelength
  alp2betap = 0.001 * c1p * alp1p 
  alp2betas = 0.001 * c1s * alp1s  

END SUBROUTINE PHI2PARM

! ----------------------------------------------------------
! SUBROUTINE  FINE_GRAINED():
! ----------------------------------------------------------
SUBROUTINE FINE_GRAINED(phim, cw_bot, SSR, alpp_0, z1, ndep,&
     c1p, c1s, r1, alp1p, alp1s)
  IMPLICIT  NONE
  INTEGER, INTENT(IN)      :: ndep  
  DOUBLE PRECISION, DIMENSION(1:ndep), INTENT(IN) ::  z1
  DOUBLE PRECISION, DIMENSION(1:ndep)  ::  x
  DOUBLE PRECISION, INTENT(IN)         :: cw_bot, phim, alpp_0, SSR
  DOUBLE PRECISION, DIMENSION(1:ndep), INTENT(OUT):: c1p, c1s, r1,alp1p,alp1s 
  DOUBLE PRECISION         :: new_phim
  INTEGER                  :: i  

  IF (phim < 5.75) THEN
     new_phim = 5.75
  ELSE
     new_phim = phim
  END IF

! comp. sound speed profile  (m/s)  
  c1p = SSR * cw_bot + 0.712*z1;

! P attenuation profile (dB/m/kHz)    
  alp1p(1) = alpp_0
  alp1p(2:ndep) = alpp_0 + 9.088d-5 * (z1(2:ndep)) - &
      2.285d-7 * ((z1(2:ndep))**2.) + 1.336d-10 * ((z1(2:ndep))**3.)

! shear sound speed profile  (m/s)  
  c1s(1) = c1p(1)/13.994                ! surficial shear s 
  x(1)   = 1.  
  x(2:ndep) = 1.063 - 0.01106486*LOG10(z1(2:ndep)) - &
       0.03604141*2*LOG10(z1(2:ndep))
  c1s(2:ndep) = c1p/(10**x(2:ndep))     ! shear s w/ depth
  
! density profile  (g/cm^3)  
  r1 = (22.85 - new_phim)/10.275 + 0.001395*z1 - 6.17e-7 * (z1**2.)
  
! shear attenuation profile (dB/m/kHz)  
  alp1s(1) = 17.3
  alp1s(2:ndep) = alp1p(2:ndep) * (alp1s(1)/alp1p(1))

END SUBROUTINE FINE_GRAINED 

! ----------------------------------------------------------
! SUBROUTINE  COARSE_GRAINED():
! ----------------------------------------------------------
SUBROUTINE COARSE_GRAINED(phim, cw_bot, SSR, alpp_0, z1, ndep, &
     c1p, c1s, r1, alp1p, alp1s)
  IMPLICIT  NONE
  INTEGER, INTENT(IN)      :: ndep  
  DOUBLE PRECISION, DIMENSION(1:ndep), INTENT(IN):: z1
  DOUBLE PRECISION, INTENT(IN)         :: cw_bot, phim, alpp_0, SSR
  DOUBLE PRECISION, DIMENSION(1:ndep), INTENT(OUT):: c1p, c1s, r1, alp1p, alp1s
  DOUBLE PRECISION         :: new_phim
  INTEGER                  :: i  

  IF (phim > 3.25) THEN
     new_phim = 3.25
  ELSE
     new_phim = phim
  END IF
  
! comp. sound speed profile (m/s) 
  c1p(1) = SSR * cw_bot
  c1p(2:ndep) = c1p(1) * ((z1(2:ndep)/0.05)**0.015)      
  
! P attenuation profile (dB/m/kHz)    
  alp1p(1) = alpp_0
  alp1p(2:ndep) = alpp_0 * ((z1(2:ndep)/0.05)**(-1./6.))  

! shear sound speed profile  (m/s)  
  c1s(1) = c1p(1)/31.4
  c1s(2:ndep) = c1s(1) * ((z1(2:ndep)/0.05)**0.28)       ! shear s w/ depth
  
! density constant w/ depth  (g/cm^3) 
  r1 = (22.85-new_phim)/10.275

! shear attenuation profile (dB/m/kHz)  
  alp1s(1) = 13.2
  alp1s(2:ndep) = alp1s(1) * ((z1(2:ndep)/0.05)**(-1./6.))

END SUBROUTINE COARSE_GRAINED

!END PROGRAM GetProfiles  

