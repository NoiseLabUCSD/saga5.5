      function besl0(x)
c
c: this function approximates a zero'th order modified bessel
c: function
c
      implicit integer*4(i-n)
      t=x/3.75
      if(abs(t) .le. 1.) then 
         tsq=t*t
         besl0=1.0 + tsq*(3.5156229 + tsq*(3.0899424 + tsq*(1.2067492 
     .       + tsq*(.2659732 + tsq*(.0360768 +tsq*.0045813)))))
      else
         besl0=.39894228 + (((((((.00392377/t - .01647633)/t
     .       + .02635537)/t - .02057706)/t + .00916281)/t
     .       - .00157565)/t + .00225319)/t + .01328592)/t
         besl0=besl0/(sqrt(x)*exp(-x))
      endif
c
      return
      end 
