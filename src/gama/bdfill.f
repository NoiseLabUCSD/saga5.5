      subroutine bdfill(kl,xl,nr,kbdf,gslbd,tbd,ebd,trmbd,trpbd)
c
c: this subroutine fills beam displaced eigenrays characteristics
c: from kl and xl.
c
      implicit integer*4(i-n)
      integer*4 kl(2,2500)
      real xl(7,2500),gslbd(20),tbd(20),ebd(20),trmbd(20),trpbd(20)
c
      gslbd(kbdf)=abs(xl(2,nr))
      tbd(kbdf)=xl(3,nr)
      ebd(kbdf)=xl(4,nr)
      trmbd(kbdf)=xl(5,nr)
      trpbd(kbdf)=xl(7,nr)
      kl(1,nr)=kl(1,nr) - kbdf*1000000
c
      return
      end
