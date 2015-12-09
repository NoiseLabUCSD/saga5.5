      function divrat(u,v,w,du,dv,dw)
c
c: this function takes the derivative of (uv)/(w) given their
c: derivative.
c
      implicit integer*4(i-n)
      divrat=(w*(u*dv+v*du) - u*v*dw)/w**2
c
      return
      end 
