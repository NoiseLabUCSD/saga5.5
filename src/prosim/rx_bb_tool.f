c:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************

      subroutine lname80(chstr,leng)
c
c: Subroutine to determine the length of a character
c: string without blanks.  Starts at END of string.
c: chstr -- the string (length 64)
c: leng -- the length without blanks
c
      character*80 chstr
c
      do 40 i=80,1,-1
       	 if(chstr(i:i) .ne. ' ') then
            leng = i
            return
       	 endif
40    continue
      leng=0
c
      return
      end
ccc
      subroutine uni_space(n,x,fac)
c
      implicit none
      include 'Parms_com'
      integer*4 n,j
      real*4 x(NVRMAX),fac,xfac
c
      if(n .lt. 0) then
         n=iabs(n)
         xfac=(x(2) - x(1))/max(1,n-1)
         do j=2,n
            x(j)=x(1) + float(j-1)*xfac
         enddo
         if(fac .ne. 1.e0) then
            do j=1,n
               x(j)=fac*x(j)
            enddo
         endif
      endif
c
      return
      end
ccc
      subroutine hpsort(n,ra)
c
      implicit none
      integer*4 n
      integer*4 i,ir,j,l
      real*4 ra(n)
c: Sorts an array ra(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted 
c: rearrangement.
c
      real*4 rra
      if(n .lt. 2) return
c: The index l will be decremented from its initial value down to 1
c: during the "hiring" (heap creation) phase.  Once it reaches 1, the
c: index ir will be decremented from its initial value down to 1
c: during the "retirement and promotion" (heap selection) phase.
      l=n/2 + 1
      ir=n
10    continue
         if(l .gt. 1) then
            l=l-1
            rra=ra(l)
         else
            rra=ra(ir)
            ra(ir)=ra(1)
            ir=ir-1
            if(ir .eq. 1) then
               ra(1)=rra
               return
            endif
         endif
         i=l
         j=l+l
20       if(j .le. ir) then
            if(j .lt. ir) then
               if(ra(j) .lt. ra(j+1)) j=j+1
            endif
            if(rra .lt. ra(j)) then
               ra(i)=ra(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         ra(i)=rra
      goto 10
c
      end
ccc
      subroutine hpsort_indx(n,ra,indx)
c
      implicit none
      integer*4 n,indx(n)
      real*8 ra(n)
      integer*4 i,ir,j,l,iia
c: Sorts an array ra(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted
c: rearrangement.
c
      real*8 rra
c
      do j=1,n
         indx(j)=j
      enddo
c
      if(n .lt. 2) return
c: The index l will be decremented from its initial value down to 1
c: during the "hiring" (heap creation) phase.  Once it reaches 1, the
c: index ir will be decremented from its initial value down to 1
c: during the "retirement and promotion" (heap selection) phase.
      l=n/2 + 1
      ir=n
10    continue
         if(l .gt. 1) then
            l=l-1
            rra=ra(l)
            iia=indx(l)
         else
            rra=ra(ir)
            iia=indx(ir)
            ra(ir)=ra(1)
            indx(ir)=indx(1)
c
            ir=ir-1
            if(ir .eq. 1) then
               ra(1)=rra
               indx(1)=iia
               return
            endif
         endif
         i=l
         j=l+l
20       if(j .le. ir) then
            if(j .lt. ir) then
               if(ra(j) .lt. ra(j+1)) j=j+1
            endif
            if(rra .lt. ra(j)) then
               ra(i)=ra(j)
               indx(i)=indx(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         ra(i)=rra
         indx(i)=iia
      goto 10
c
      end
ccc
      subroutine r_hpsort_indx(n,ra,indx)
c
      implicit none
      integer*4 n,indx(n)
      real*4 ra(n)
      integer*4 i,ir,j,l,iia
c: Sorts an array ra(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted
c: rearrangement.
c
      real*4 rra
c
      do j=1,n
         indx(j)=j
      enddo
c
      if(n .lt. 2) return
c: The index l will be decremented from its initial value down to 1
c: during the "hiring" (heap creation) phase.  Once it reaches 1, the
c: index ir will be decremented from its initial value down to 1
c: during the "retirement and promotion" (heap selection) phase.
      l=n/2 + 1
      ir=n
10    continue
         if(l .gt. 1) then
            l=l-1
            rra=ra(l)
            iia=indx(l)
         else
            rra=ra(ir)
            iia=indx(ir)
            ra(ir)=ra(1)
            indx(ir)=indx(1)
c
            ir=ir-1
            if(ir .eq. 1) then
               ra(1)=rra
               indx(1)=iia
               return
            endif
         endif
         i=l
         j=l+l
20       if(j .le. ir) then
            if(j .lt. ir) then
               if(ra(j) .lt. ra(j+1)) j=j+1
            endif
            if(rra .lt. ra(j)) then
               ra(i)=ra(j)
               indx(i)=indx(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         ra(i)=rra
         indx(i)=iia
      goto 10
c
      end
ccc
      subroutine mem_lim(n,nlim,eline,le,vname,lv,lname,ll,iibad,
     .   iistop)
c
      implicit none
      include 'Parms_com'
c
      integer n,nlim,le,lv,ll,iibad,iistop
      character*128 eline,vname,lname
c
      if(n .gt. nlim) then
         print *,' '
         print *,eline(1:le)
         print *,'VARIABLE NAME = ',vname(1:lv),'; LIMIT = ',nlim,
     .      '; LIMIT NAME = ',lname(1:ll)
         print *,'ENTERED OR COMPUTED VALUE FOR THIS RUN = ',n
cc       do j=1,nlay+1
cc          print *,j,1,h(j),(geo(1,jj,j),jj=1,5)
cc          print *,j,2,h(j),(geo(2,jj,j),jj=1,5)
cc       enddo
         if(iistop .eq. 1) stop
         iibad=1
      endif
c
      return
      end
ccc
      subroutine var_lim(n,nlim,eline,le,vname,lv,ibad,istop)
c
      implicit none
      include 'Parms_com'
c
      integer n,nlim,le,lv,istop,ibad
      character*128 eline,vname
c
      ibad = 0
      if(n .le. nlim) then
         print *,' '
         print *,eline(1:le)
         print *,'VARIABLE NAME = ',vname(1:lv),'; LIMIT = ',nlim
         print *,'ENTERED OR COMPUTED VALUE FOR THIS RUN = ',n
         if(istop .eq. 1) stop
         ibad = 1
      endif
c
      return
      end

ccc
      SUBROUTINE zroot8_int(a,m,xroot,xlo,xhi)
      implicit none
      INTEGER m,MAXM
      REAL*8 EPS,xlo,xhi
      COMPLEX*16 a(m+1),xroot,roots(20)
      PARAMETER (EPS=1.d-6,MAXM=101)
CU    USES laguer8
      INTEGER j,jj,its,nf
      COMPLEX*16 ad(MAXM),x,b,c
      do 11 j=1,m+1
        ad(j)=a(j)
11    continue
      nf=0
      x=xroot
      do 13 j=m,1,-1
cc      x=cmplx(0.,0.)
        call laguer8(ad,j,x,its)
c: Check if root found is in interval and essentially real:
        if(dreal(x) .ge. xlo .and. dreal(x) .le. xhi .and.
     .     dabs(dimag(x)) .lt. EPS) then
c: Polish if necessary
           if(j .lt. m) then
              call laguer8(a,m,x,its)
           endif
           xroot=dcmplx(dreal(x),0.d0)
      nf=nf + 1
           if(j .lt. m) print *,'root not first: ',m-j+1
c          return
        endif
cc      if(abs(dimag(x)).le.2.*EPS**2*abs(dreal(x)))
cc   .     x=dcmplx(real(x),0.d0)
        roots(j)=x
        b=ad(j+1)
        do 12 jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
12      continue
        x=dcmplx(0.5d0,0.d0)
13    continue
cc    print *,'no roots found between xlo and xhi: ',xlo,xhi,roots
      if(nf .gt. 1) print *,'>1 root found: ',roots,xroot
      return
      END
ccc
      SUBROUTINE laguer8(a,m,x,its)
      INTEGER m,its,MAXIT,MR,MT
      REAL*8 EPSS
      COMPLEX*16 a(m+1),x
      PARAMETER (EPSS=2.e-7,MR=8,MT=10,MAXIT=MT*MR)
      INTEGER iter,j
      REAL*8 abx,abp,abm,err,frac(MR)
      COMPLEX*16 dx,x1,b,d,f,g,h,sq,gp,gm,g2
      SAVE frac
      DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
      do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=abs(b)
        d=dcmplx(0.d0,0.d0)
        f=dcmplx(0.d0,0.d0)
        abx=abs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=abs(b)+abx*err
11      continue
        err=EPSS*err
        if(abs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.*f/b
          sq=cdsqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          abp=abs(gp)
          abm=abs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.) then
            dx=m/gp
          else
            dx=cdexp(dcmplx(dlog(1.d0+abx),dble(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1) return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
12    continue
      pause 'too many iterations in laguer'
      return
      END
ccc
