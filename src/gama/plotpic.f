      subroutine plotpic(range,jjnr)
c
c: this subroutine plots a picture of the ocean and the layered bottom.
c
      implicit integer*4(i-n)
      include 'common/depth'
      include 'common/srlox'
      include 'common/pii'
      include 'common/pic'
      include 'common/charcom'
      real range(nrtot)
      integer*4 jjnr(nrtot)
c
      rmax=0.
      do 10 jr=1,nrtot
         rmax=max(rmax,range(jr))
10    continue
c
      write(55,200) 0.,ctab,0.,ctab,2
      write(55,200) 0.,ctab,-zlev(nlev+4),ctab,0
      write(55,200) .001*rmax,ctab,-zlev(nlev+4),ctab,0
      write(55,200) .001*rmax,ctab,0.,ctab,0
      write(55,200) 0.,ctab,0.,ctab,0
200   format(f9.4,a1,f8.2,a1,i1)
      do 20 j=4,nlev+3
         write(55,200) 0.,ctab,-zlev(j),ctab,0
         write(55,200) .001*rmax,ctab,-zlev(j),ctab,0
         write(55,200) 0.,ctab,-zlev(j),ctab,0
20    continue
c
      do 28 jzr=1,nzr
         zr=xyz(kxy(jzr)+1,3)
         write(55,200) 0.,ctab,-zr,ctab,3
         write(55,200) 0.,ctab,-zs,ctab,0
         nnxy=nxy(jzr)
         do 29 jxy=1,nnxy
            do 30 kr=1,nrangx
               mr=jjnr(mr1(jzr) + nnxy*(kr-1) + jxy - 1) 
               write(55,200) .001*range(mr),ctab,-zs,ctab,4
30          continue
29       continue
         write(55,200) 0.,ctab,-zs,ctab,0
28    continue
c
      return
      end 
