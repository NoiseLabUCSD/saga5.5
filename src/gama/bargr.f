      subroutine bargr(nr,jcaus,t,e,db,a,phr)
c
c: this subroutine plots a bar graph for each desired range.  the
c: graph gives the arrival time, total attenuation in db,
c: and frequency-dependent attenuation factor for each eigenray found.
c
      implicit integer*4(i-n)
      include 'common/vhfacs' 
      include 'common/paths'
      include 'common/pathchr'
      include 'common/pathway'
      include 'common/pii'
      include 'common/svp'
c
      character*1 cq(0:3)
      data cq/' ','$','&','!'/
c
c: determine sign of receiver grazing angle (+ from above, - from below:
      sg=1.
      if(ldir(icpth) .eq. 'u') sg=-1.
      sipt=(sg*acos(a*cr)*piedeg - silo)*vlfac2
      tpt=(t-tlo)*hfac
c: db=999.99 means that the magnitude of the ray is essentially zero. 
      if(db .ne. 999.99) then 
         dbpt=(db-dblo)*vhfac 
      else
         dbpt=0.
      endif
      ept=e*vlfac
cxx   call symbol(tpt,dimv/2.+dbpt,.12,15,0.0,-1) 
cxx   call symbol(tpt,ept,.12,15,0.0,-2)
cxx   call symbol(tpt,sipt,.12,4,0.0,-1)
c: 8-25-89: put a triangle symbol for the total phase of the ray: 
c     phpt=(phr + 180.)*dimv/360.
c     call symbol(tpt,phpt,.12,2,0.,-1) 
c
cxx   call pltline(tpt,dimv/2.+dbpt+.08,-.12)
cxx   if(nr .lt. 10) then
cxx      write(8,110) nr
cxx   elseif(nr .lt. 100) then
cxx      write(8,130) nr
cxx   elseif(nr .ge. 100) then
cxx      write(8,130) nr
cxx   endif
cxx   call pltline(tpt,ept-.24,-.16)
cxx   write(8,150) cq(jcaus)
150   format(a1)
c
110   format(i1)
120   format(i2)
130   format(i3)
c
      return
      end 
