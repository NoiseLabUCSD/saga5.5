      subroutine lnfit(k,g1,g2,a,b)
c
c: this subroutine does a linear fit g=af + b to the two data points
c: (bdf(k-1),g1) and (bdf(k),g2).
      implicit integer*4(i-n)
      include 'common/bdcom'
c
      a=(g2-g1)/(bdf(k)-bdf(k-1))
c     a=(g2-g1)/(log10(bdf(k))-log10(bdf(k-1)))
      b=g2 - a*bdf(k)
c     b=g2 - a*log10(bdf(k))
c
      return
      end
