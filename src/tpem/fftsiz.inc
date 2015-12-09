c MXNFFT:  Maximum power of 2 for transform size
c MAXPTS:  Maximum size of arrays for the real and imaginary fields
      
      integer*4 maxpts, mxnfft
      parameter ( mxnfft = 12)
      parameter ( maxpts = 2**mxnfft ) 
C REMEMBERER ALSO TO UPDATE THESE VARIABLE IN TPEM.INC
