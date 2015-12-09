      subroutine traduc(icpth,fd,ktop,kbot,cz)
c
c: this subroutine translates the user's ocean path specification
c: into data the program needs to find rays and plot them.  the
c: ocean is divided into three sections (top, middle, and bottom)
c: by the source and receiver depths.  an ocean ray path is 
c: characterized by its initial direction from the source (fdir=
c: 'u' [up] or 'd' [down]), the number of times the ray enters the
c: top section (#t), and the number of times the ray enters the
c: bottom section (#b), and the ray type ('r'=refracting in the
c: bottom ocean section; 'w'=waterborne and totally reflected at
c: ocean bottom; 'p'=bottom penetrating; 'a'=all).  the array cpth is 
c: such that ncpth(icpth,j) is the number of legs the icpth'th ocean
c: path has in the j'th ocean section.
c: the ocean paths j=196,197,198 have a single leg in one of the 3
c: sections and are used to calculate the r'-a diagrams for each of
c: the 3 ocean sections.
c
      implicit integer*4(i-n)
      include 'common/pathchr'
      include 'common/paths'
      include 'common/depth'
      character*1 fd,cz
      data nbmax/0/ 
c
      data ncpth(196,1)/1/
      data ncpth(196,2)/0/
      data ncpth(196,3)/0/
      data ncpth(197,1)/0/
      data ncpth(197,2)/1/
      data ncpth(197,3)/0/
      data ncpth(198,1)/0/
      data ncpth(198,2)/0/
      data ncpth(198,3)/1/
c
      j=2 
      if(zs .lt. zr) j=1
c
      if((cz .ne. 'r') .and. (cz .ne. 'p') .and. (cz .ne. 'w')
     .   .and. (cz .ne. 'a')) call cperr(icpth,fd,ktop,kbot,cz)
      czprop(icpth)=cz
c
      if((fd .ne. 'u') .and. (fd .ne. 'd')) call
     .      cperr(icpth,fd,ktop,kbot,cz)
      kdiff=ktop-kbot
      if((fd .eq. 'u') .and. (kdiff .ne. 0) .and. (kdiff .ne. 1))
     .      call cperr(icpth,fd,ktop,kbot,cz)
      if((fd .eq. 'd') .and. (kdiff .ne. 0) .and. (kdiff .ne. -1))
     .      call cperr(icpth,fd,ktop,kbot,cz)
      ncpth(icpth,2)=(max0(ktop,kbot)-1)*2 + 1
      if((kdiff .eq. 0) .and. (((fd .eq. 'u') .and. (j .eq. 2))
     .      .or. ((fd .eq. 'd') .and. (j .eq. 1)))) ncpth(icpth,2)=
     .      ncpth(icpth,2) + 2
      if(ncpth(icpth,2) .lt. 0) then
c        print *,'direct ocean path included.'
         ncpth(icpth,2)=1
         if(fd .eq. 'u') then 
            fd='d'
         else
            fd='u'
         endif
      endif
      ncpth(icpth,1)=2*ktop
      ncpth(icpth,3)=2*kbot
      fdir(icpth)=fd
      ldir(icpth)='d'
      if((kbot .gt. ktop) .or. ((kbot .eq. ktop) .and. (fd .eq. 'u')))
     .   ldir(icpth)='u'
      nbmax=max0(nbmax,kbot)
c
      return
      end 
