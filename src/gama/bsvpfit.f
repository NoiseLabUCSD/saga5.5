      subroutine bsvpfit(z1,z2,c1,c2,b1,b2)
c
c: this subroutine finds the slope b2, given the depths z1 and z2,
c: the sound speeds c1 and c2, and the slope b1 at z1.  The 
c: criterion is based on making the radicand for the g2 parameter
c: of the Weinberg profile positive, but close to zero.
c
      implicit integer*4(i-n)
      dz=z2 - z1
      dc=c2 - c1
      b2=.999*(dc*(c1+c2)/(2.*dz))**2 / (c1*c2*b1)
c
      return
      end
