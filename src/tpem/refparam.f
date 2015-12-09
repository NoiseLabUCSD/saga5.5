
      subroutine rdpara( rf,rp,sv, range,rmax,iperror)

      include 'tpem.inc'

      common / parinit / rv2, refdum(mxlvls), htdum(mxlvls),
     +  profint(0:maxpts), ht(0:maxpts), is
       
      TYPE(refractivity) rf
      TYPE(refparam) rp
      TYPE(systemvar) sv 
     
      real coef(mpcoef)
      real base,thick,offset,mdef
      real xx,pol(10)
      real*8 rmax

      integer nmarkov,mmarkov
      parameter( nmarkov=100, mmarkov=6)
      real markov( nmarkov,mmarkov)
      real remain,dum
      integer ip      

      common /markovchain/markov

c If there is a range-dependent refrac prof then interpolate horizontally
c using the two surrounding param at range RANGE with all duplicate levels.
  
c      xx=2*range/rmax -1
c      do i=1,rp.npoly
c         pol(i)=xx**(i-1)
c      enddo
c      pol(1)=1
c      pol(2)=   xx
c      pol(3)= 2*xx**2- 1
c      pol(4)= 4*xx**3- 3*xx
c      pol(5)= 8*xx**4- 8*xx**2+1
c      pol(6)=16*xx**5-20*xx**3+5*xx

       if (markov(1,1).ne.1) then
          open(unit=31,file='markov.in')
c     read(31,*) nmarkov,mmarkov  ! number of coef and 
          do ip=1,nmarkov
             read(31,*) dum,(markov(ip,i),i=2,mmarkov)
             markov(ip,1)=1
c             write(*,*)'Markov parameters'
c             write(*,*) dum,(markov(ip,i),i=2,mmarkov)
          enddo
          close(31)
       endif
       xx=range/rmax*(nmarkov-1)+1
       ip=xx
       remain=xx-ip
       if (ip.eq.nmarkov) then
          ip =ip-1
          remain=1
       endif
       do i=1,5
         pol(i)=markov(ip,i)*(1-remain) +markov(ip+1,i)*remain
       enddo
c      write(*,*)'range,xx,remain,ip',range,xx,remain,ip
c      write(*,*)pol(1),pol(2),pol(3),pol(4)
c      write(*,*)markov(ip,1),markov(ip+1,1)
      base=0.
      do i=1,rp%npoly
         base=base+pol(i)*rp%factor(i,1)
      enddo
c         write(*,*) ' Base =',base
      if (base.lt.0) then
c         write(*,*) ' Base =',base
c         write(*,*) ' Base has been truncated '
         base=0
      endif
      
         thick = rp%thick(1) 
         offset= rp% offset(1)
         mdef  = rp%mdef(1)
c        write(90,*)range,base
         do i=1,mpcoef
           coef(i)= rp%coef(i,1)
         enddo
            call paramprofile(base,thick,
     1           offset,mdef,coef,sv%antht,iperror)
           
      end    
c
c***************************
c

      subroutine paraminter( rf,rp,sv, range,iperror)

      include 'tpem.inc'

      common / parinit / rv2, refdum(mxlvls), htdum(mxlvls),
     +  profint(0:maxpts), ht(0:maxpts), is
       
      TYPE(refractivity) rf
      TYPE(refparam) rp
      TYPE(systemvar) sv       
      integer iperror
      save j, rv1
      data j, rv1 / 0, 0. /
      real coef(mpcoef)
      real base,thick,offset,mdef
c One-line interpolation function
c      pint( p1, p2 ) = p1 + fv * ( p2 - p1 ) 
      
c If there is a range-dependent refrac prof then interpolate horizontally
c using the two surrounding param at range RANGE with all duplicate levels.
  
      if( rf%nprof .gt. 1 ) then
 10      if ( (range.gt.rv2) .and. (is.lt. rf%nprof)) then
            j = is
            IS=IS+1
            rv1=rv2
            rv2=rf%rngprof(IS)
            goto 10
         end if
      

         if (range.gt.rv2) then
             fv=1
             fv2=0
         else
            FV =(range-rv1)/(rv2-rv1)  ! left point
            FV2=(rv2-range)/(rv2-rv1)  ! right point
         endif
c         write(*,*)'range,j,is',range,j,is,fv,fv2
         base  = rp%base(j)   *fv2+rp%base(is)  *fv
         thick = rp%thick(j)  *fv2+rp%thick(is) *fv
         offset= rp% offset(j)*fv2+rp%offset(is)*fv
         mdef  = rp%mdef(j)   *fv2+rp%mdef(is)  *fv
         do i=1,mpcoef
           coef(i)= rp%coef(i,j)*fv2+ rp%coef(i,is)  *fv
         enddo

            call paramprofile(base,thick,
     1           offset,mdef,coef,sv%antht,iperror)
           
      end if 
      end    
c####################################
      subroutine  paramprofile(baseheight,thick,offset,mdef,coef,
     &               zs,iperror)
c used for defining profile      
      real coef(*),zs
      real c0,c1,c2,delta,baseheight,thick,offset,Mdef,zup
      integer iperror
      include 'tpem.inc'

      common / pevar / wl, fko, delz, n, ln, zmax, n34, con
      common / parinit / rv2, refdum(mxlvls), htdum(mxlvls),
     +  profint(0:maxpts), ht(0:maxpts), is
      integer  lvlep_start,itpem_opt(40)
      common  /sagatpem/lvlep_start,itpem_opt
      integer iprof
      real mexcess,metcon,basehvar(10000)
      common/ mexc/ mexcess,metcon,basehvar,iprof
c local variables
      integer iz,i,izh
      real zr,mr,zloc,z(0:maxpts)
      real ahelp
c      write(*,*)'entering  paramprofile'
      iperror=0
c      c1 =0.118    !0.13         ! slope in evaporation and mixed layer adiabatic
c      c2 = 0.118       ! slope above inversion layer
c      delta=50        ! evaporation duct
c      z0 =0.00015      ! roughness
      c0=0.13
      c1=coef(1)
      c2=coef(2)
      delta=coef(3)
c      baseheight=250  !
c      thick=50
c      offset=339
c      Mdef= 49.5
c       write(*,*)'***************************************'
c       write(99,*)'baseheight, thick, offset, Mdef'
c       write(99,*)baseheight, thick, offset, Mdef
c       write(99,*)'slopes in mixed and top layer,delta',c1,c2,delta
c       write(99,*)'delz',delz
c      zmax=500      ! maximum heigth of PE

c       write(99,*)'delz',delz

      zup=zmax-baseheight-thick
      if (zup.lt.0 ) then
        if (baseheight.gt.zmax ) then
            baseheight=zmax
            thick=0
         else
            thick= zmax-baseheight
         endif
      endif
c      dz=0.7  ! pe discretization
c evaporation duct
      ahelp=1/(1-c1/c0)
      if ((ahelp.gt.0).and. (ahelp.lt.2)) then
        zd=ahelp*delta
      else
        zd=2*delta
      endif
        zd=min(zd, baseheight)
c      write(*,*)'Evaporation duct ',zd
      iz=((zd+delz/2)/delz)
      do i=0,iz
        zloc=delz*(i+0.01)
        z(i)=zloc
        profint(i)=c0*(z(i)-zd-delta*log(z(i)/zd))
      enddo
      if (iz.eq.0) then
        profint(iz)=0
      endif
      Mr = profint(iz) 
      zr = z(iz)

c mixed layer
      izh=(( baseheight-zd+delz/2)-(zr-zd))/delz
      do i=1,izh
         ih=i+iz
         zloc=zr-zd+delz*i
         z(ih)=zd+zloc
         profint(ih)=c1*z(i)
      enddo
      iz = iz+izh
      Mr = profint(iz) 
      zr = z(iz)
c      write(99,*)'mixed layer',zr,c1
c      write(*,*)'mixed layer',zr,c1
c
c inversion layer
c
c      zh = baseheight+thick+delz/2
      dM=-Mdef/(thick)
      izh=((thick+delz/2)-(zr-baseheight))/delz
c      write(*,*)'izh,Mr,zr,zr-baseheight'
c      write(*,*)izh,Mr,zr,zr-baseheight
      do i=1,izh
         ih=i+iz
         zloc =zr-baseheight+delz*i
         z(ih) = baseheight+zloc
         profint(ih)=Mr+dM*zloc
      enddo
      iz = iz+izh
      Mr = profint(iz)
      zr = z(iz)
c      write(99,*)'inversion layer top,slope',zr,dM
c      write(*,*)'inversion layer top,slope',zr,dM

C top layer

      izh=N-iz
      do i=1,izh
         ih=i+iz
         zloc =zr-baseheight-thick+delz*i
         z(ih) = baseheight+thick+zloc
         profint(ih)=Mr+c2*zloc
      enddo

       refmin=1000
       do i=0,N
          if ( profint(i).lt.refmin) then
             refmin=profint(i)
             iref=i
          endif
       enddo
        
       iprof= iprof+1
       basehvar(iprof)=baseheight
c        write(*,*)'iprof',iprof
c       write(*,*)'pg iref, refmin',iref, refmin,zs,profint(0)
c       write(*,*)'pg2',  profint(int(zs/delz+1)),int(zs/delz+1),delz
         refmin=refmin-0.118*z(iref)
c             write(*,*)'pg iref, refmin',iref, refmin,
c     1        profint(int(zs/delz+1))-profint(iref)
c
c  Meteorological constrain
c
       metcon=min(refmin, metcon)
       if ((itpem_opt(8).eq.1).and. (refmin.lt. -60)) then 
             write(*,*)'pg  refmin',iref, refmin,
     1        profint(int(zs/delz+1))-profint(iref)
         write(7,*)' profile rejected, refmin less than -60', 
     1               refmin, refmin,iref
         write(*,*)' profile rejected, refmin less than -60', 
     1               refmin, refmin,iref
         iperror=-1
       endif
c
c  Ref coef  constrain
c
      mexcess=max(mexcess,profint(int(zs/delz+1))-profint(iref))
      if ((itpem_opt(9).eq.1).and. 
     1   (mexcess.gt. 14)) then 
             write(*,*)'pg Mexcess',iref, refmin,
     1      mexcess
         write(7,*)' profile rejected Mexcess greter than 14',
     1             refmin,iref
         write(*,*)' profile rejected Mexcess greter than 14',
     1             refmin,iref
         iperror=-1
       else 
c         write(*,*)'accepted'
       endif

c add offset
       do i=0,N
         profint(i)= (profint(i)+offset)*con
       enddo
c write out
c       do i=0,N,N
c         write(99,*)i,z(i),profint(i)/con        
c       enddo

         
       end


 
