C
C     ARRAYS FOR DEPTH-DEPENDENT GREEN'S FUNCTIONS
C     AND WORKING ARRAYS
C
      COMPLEX CFF(NP,3)
      COMPLEX CFFS(NP)
      COMMON /STRDIS/ CFF,CFFS
      COMPLEX CBUF(NP),CFILE(ISIZE)
      REAL ARG(NP),FAC(NP)
      COMMON /BUFF1/ ARG,FAC
      COMMON /BUFF2/ CBUF
      COMMON /BUFF3/ CFILE
c************ This matrix contains the wave number integrands
C      complex  wavenoint(1,1)
      complex  wavenoint(np,nrd),facsqrt(np)
      common /wavenoint/wavenoint,facsqrt
