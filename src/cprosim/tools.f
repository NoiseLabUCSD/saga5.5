      subroutine lname64(chstr,leng)
c
c: Subroutine to determine the length of a character
c: string without blanks
c: chstr -- the string (length 64)
c: leng -- the length without blanks
c
      character*64 chstr,blank
      character*1 char1
      data blank/'                                                  
     .              '/
c
      do 40 i=1,64
         char1 = chstr(i:i)
       	 if(char1 .eq. ' ') then
            leng = i-1
            chstr(i:64)=blank(i:64)
            return
       	 endif
40    continue
      leng=64
c
      return
      end
ccc
      subroutine lname80(chstr,leng)
c
c: Subroutine to determine the length of a character
c: string without blanks.  Starts at END of string.
c: chstr -- the string (length 80)
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
      subroutine usage(progname,lprog)
c
c: This subroutine updates a file containing a list of the times
c: certain programs have been run and the users that ran them.
c
      integer kdat(3),ktim(3)
      character*20 progname
      character*64 home
c
      open(83,file='/usr/people/westwood/.usage',status='unknown',
     .     form='formatted')
c     .   access='append')
      call getenv('HOME',home)
c      call idate(kdat(1),kdat(2),kdat(3))
c     nhq 20151005
      call idate(kdat)
      call itime(ktim)
      write(83,100) progname(1:20),home(1:20),
     .   kdat(2),kdat(1),kdat(3),(ktim(j),j=1,3)
100   format(a20,'; ',a20,'; ',2(i2.2,'/'),i4,'; ',i2.2,2(':',i2.2))
      close(83)
c
      return
      end
ccc
      subroutine hpsort_re_im(n,ra_re,ra_im,indx)
c
      implicit none
      integer*4 n,indx(n)
      integer*4 i,ir,j,l,iia
      real*8 ra_re(n),ra_im(n)
c: Sorts an array ra_re(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra_re is replaced on output by its sorted 
c: rearrangement. ra_im is rearranged in the same manner as ra_re is.
c
      real*8 rra,rrb
c
      do j=1,n
         indx(j)=j
      enddo
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
            rra=ra_re(l)
            rrb=ra_im(l)
            iia=indx(l)
         else
            rra=ra_re(ir)
            rrb=ra_im(ir)
            iia=indx(ir)
            ra_re(ir)=ra_re(1)
            ra_im(ir)=ra_im(1)
            indx(ir)=indx(1)
            ir=ir-1
            if(ir .eq. 1) then
               ra_re(1)=rra
               ra_im(1)=rrb
               indx(1)=iia
               return
            endif
         endif
         i=l
         j=l+l
20       if(j .le. ir) then
            if(j .lt. ir) then
               if(ra_re(j) .lt. ra_re(j+1)) j=j+1
            endif
            if(rra .lt. ra_re(j)) then
               ra_re(i)=ra_re(j)
               ra_im(i)=ra_im(j)
               indx(i)=indx(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         ra_re(i)=rra
         ra_im(i)=rrb
         indx(i)=iia
      goto 10
c
      end
ccc
      subroutine hpsort(n,ra)
c
      implicit none
      integer*4 n
      integer*4 i,ir,j,l
      real*8 ra(n)
c: Sorts an array ra(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted 
c: rearrangement.
c
      real*8 rra
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
      subroutine hpsort_r4(n,ra)
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
      do l=1,n
         write(6,*)l,ra(l)
      end do
c
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
      subroutine hpsort_indx_c16(n,ra,indx)
c
      implicit none
      integer*4 n,indx(n)
      complex*16 ra(n)
      integer*4 i,ir,j,l,iia
c: Sorts real part of complex*16 array ra(1:n) into DESCENDING order 
c: using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted
c: rearrangement.
c
      complex*16 rra
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
               if(dble(ra(j)) .gt. dble(ra(j+1))) j=j+1
            endif
            if(dble(rra) .gt. dble(ra(j))) then
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
      subroutine hpsort_i4(n,ra)
c
      implicit none
      integer*4 n
      integer*4 ra(n)
      integer*4 i,ir,j,l
c: Sorts an integer*4 array ra(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted
c: rearrangement.
c
      integer*4 rra
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
         else
            rra=ra(ir)
            ra(ir)=ra(1)
c
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
      subroutine hpsort_i4_indx(n,ra,indx)
c
      implicit none
      integer*4 n,indx(n)
      integer*4 ra(n)
      integer*4 i,ir,j,l
c: Sorts an integer*4 array ra(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted
c: rearrangement.
c
      integer*4 rra,iia
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
      subroutine openfftout(nfile,fftroot,lfft,fftfile,nfft)
c
c: This subroutine opens an FFT file with root name fftroot and 
c: suffix _fft for output.  Opens as direct access file.
c
      character*64 fftroot,fftfile
      logical qopen
c
c: open input _fft file:
      fftfile(1:lfft+4)=fftroot(1:lfft)//'_fft'
      lfft=lfft + 4
      inquire(nfile,opened=qopen)
      if(qopen) close(nfile)
c: reopen as direct access file with record length = 22+nfft words:
can   lenrec=4*(20+(nfft+2))
c: IRIS uses words rather than bytes:
      lenrec=(20+(nfft+2))
      open(nfile,file=fftfile(1:lfft),status='unknown',
     .   form='unformatted',
     .   access='direct',recl=lenrec)
c
      return
      end
ccc
      subroutine openfftout2(nfile,fftroot,lroot,fftfile,lfft,nfft,
     .   fs,fmin,fmax,xh,NRECL)
c
c: This subroutine opens an FFT file with root name fftroot and 
c: suffix _fft for output.  Opens as direct access file.
c: For fmin=-999., it uses the old FFT format, where all nff2/2 + 1
c: frequency bins are output. Otherwise, only bins from fmin to fmax 
c: are output.  xh(7)=1 means one frequency band from xh(15)*df to 
c: xh(16)*df is included in the file, where df=fs/nfft.  The total 
c: number of bins present is xh(16)-xh(15)+1.
c
      implicit none
      integer*4 nfile,lroot,lfft,nfft,NRECL,nf1,nf2,nffth1,lenrec,
     .   total_bins
      real*4 fs,fmin,fmax,xh(20),df
      character*64 fftroot,fftfile
c
      xh(5)=nfft
      xh(6)=fs
c: open input _fft file:
      fftfile(1:lroot+4)=fftroot(1:lroot)//'_fft'
      lfft=lroot + 4
c: Open as direct access file with record length = 22+nfft words:
c: SUN, Alliant has record length in bytes:
      df=fs/nfft
      nf1=nint(fmin/df) + 1
      nf2=nint(fmax/df) + 1
      nffth1=nfft/2 + 1
c: Set fmin=-999. to output all bins in the old way:
      if(fmin .eq. -999. .or. (nf1 .eq. 1 .and. nf2 .eq. nffth1)) then
         total_bins=nffth1
         xh(7)=0.
      else
         if(nf1 .lt. 0 .or. nf2 .gt. nffth1 .or. nf1 .gt. nf2) then
            print *,'Bad fmin,fmax,fs in openfftout2: ',fmin,fmax,fs
            stop
         endif
         total_bins=nf2-nf1+1
         xh(7)=1.
         xh(15)=nf1-1
         xh(16)=nf2-1
      endif
c: For direct acces files, NRECL should be 4 for lenrec in bytes
c: (SUN), 1 for lenrec in words (IRIS).  Set in Parms_com.
      lenrec=NRECL*(20 + 2*total_bins)
      open(nfile,file=fftfile(1:lfft),status='unknown',
     .   form='unformatted',
     .   access='direct',recl=lenrec)
c
      return
      end
ccc
      subroutine mem_lim(n,nlim,eline,le,vname,lv,lname,ll,iibad,
     .   iistop)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
c      include 'i_o_com'
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
      subroutine uni_space(n,x,fac)
c
      implicit none
      integer*4 n,j
      real*4 x(10000),fac,xfac
c
      if(abs(n) .gt. 10000) then
         write(6,*)'Stop in uni_space; array index exceeded'
      end if
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
      subroutine suffix(root,lroot,j,n,suf,lsuf,root2,lroot2)
c
      implicit none
      integer*4 lroot,j,n,lsuf,lroot2,ln
      character*128 root,suf,root2,nsuf
c
      if(n .eq. 1) then
         lroot2=lroot
         root2(1:lroot)=root(1:lroot)
      else
         if(n .lt. 10) then
            write(nsuf(1:1),'(i1)') j
            ln=1
         elseif(n .lt. 100) then
            write(nsuf(1:2),'(i2.2)') j
            ln=2
         elseif(n .lt. 1000) then
            write(nsuf(1:3),'(i3.3)') j
            ln=3
         elseif(n .lt. 10000) then
            write(nsuf(1:4),'(i4.4)') j
            ln=4
         else
            print *,'error in suffix: ',j,n
         endif
         lroot2=lroot + lsuf + ln
         root2(1:lroot2)=root(1:lroot)//suf(1:lsuf)//nsuf(1:ln)
      endif
c
      return
      end
ccc
      subroutine ludcmp_r8(a,n,np,indx,d)
c
c: Modified for real*8 matrices and vectors.
c: To solve A x = b the first time:
c: call ludcmp(a,n,np,indx,d)  [a destroyed here]
c: call lubksb(a,n,np,indx,b)  [b replaced by answer x]
c: Answer x will be returned in b.
c: To solve A x = b for the same A, but different b, just call
c: call lubksb(a,n,np,indx,b) 
c: with a and indx as computed by ludcmp the first time. Again, b is replaced
c: by the answer x.
c
      INTEGER n,np,indx(n),NMAX
      real*8 a(np,np),sum,cdum,TINY
      REAL*8 d
      PARAMETER (NMAX=4,TINY=(1.0d-20,0.d0))
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
c: Take out check for zeros:
cc      if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            cdum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=cdum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j) .eq. 0.d0) a(j,j)=TINY
        if(j.ne.n)then
          cdum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*cdum
18        continue
        endif
19    continue
      return
      END
ccc
      subroutine lubksb_r8(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      real*8 a(np,np),b(n),sum
      INTEGER i,ii,j,ll
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
ccc
      subroutine ludcmp(a,n,np,indx,d)
c
c: Modified for complex*16 matrices and vectors.
c: To solve A x = b the first time:
c: call ludcmp(a,n,np,indx,d)  [a destroyed here]
c: call lubksb(a,n,np,indx,b)  [b replaced by answer x]
c: Answer x will be returned in b.
c: To solve A x = b for the same A, but different b, just call
c: call lubksb(a,n,np,indx,b) 
c: with a and indx as computed by ludcmp the first time. Again, b is replaced
c: by the answer x.
c
      INTEGER n,np,indx(n),NMAX
      complex*16 a(np,np),sum,cdum,TINY
      REAL*8 d
      PARAMETER (NMAX=4,TINY=(1.0d-20,0.d0))
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (cdabs(a(i,j)).gt.aamax) aamax=cdabs(a(i,j))
11      continue
c: Take out check for zeros:
cc      if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*cdabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            cdum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=cdum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j) .eq. dcmplx(0.d0,0.d0)) a(j,j)=TINY
        if(j.ne.n)then
          cdum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*cdum
18        continue
        endif
19    continue
      return
      END
ccc
      subroutine lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      complex*16 a(np,np),b(n),sum
      INTEGER i,ii,j,ll
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
ccc
      SUBROUTINE zroots(a,m,roots,polish)
      INTEGER m,MAXM
      REAL EPS
      COMPLEX a(m+1),roots(m)
      LOGICAL polish
      PARAMETER (EPS=1.e-6,MAXM=101)
CU    USES laguer
      INTEGER i,j,jj,its
      COMPLEX ad(MAXM),x,b,c
      do 11 j=1,m+1
        ad(j)=a(j)
11    continue
      do 13 j=m,1,-1
        x=cmplx(0.,0.)
        call laguer(ad,j,x,its)
        if(abs(aimag(x)).le.2.*EPS**2*abs(real(x))) x=cmplx(real(x),0.)
        roots(j)=x
        b=ad(j+1)
        do 12 jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
12      continue
13    continue
      if (polish) then
        do 14 j=1,m
          call laguer(a,m,roots(j),its)
14      continue
      endif
      do 16 j=2,m
        x=roots(j)
        do 15 i=j-1,1,-1
          if(real(roots(i)).le.real(x))goto 10
          roots(i+1)=roots(i)
15      continue
        i=0
10      roots(i+1)=x
16    continue
      return
      END
ccc
      SUBROUTINE laguer(a,m,x,its)
      INTEGER m,its,MAXIT,MR,MT
      REAL EPSS
      COMPLEX a(m+1),x
      PARAMETER (EPSS=2.e-7,MR=8,MT=10,MAXIT=MT*MR)
      INTEGER iter,j
      REAL abx,abp,abm,err,frac(MR)
      COMPLEX dx,x1,b,d,f,g,h,sq,gp,gm,g2
      SAVE frac
      DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
      do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=abs(b)
        d=cmplx(0.,0.)
        f=cmplx(0.,0.)
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
          sq=sqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          abp=abs(gp)
          abm=abs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.) then
            dx=m/gp
          else
            dx=exp(cmplx(log(1.+abx),real(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1)return
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
        if(dble(x) .ge. xlo .and. dble(x) .le. xhi .and.
     .     dabs(dimag(x)) .lt. EPS) then
c: Polish if necessary
           if(j .lt. m) then
              call laguer8(a,m,x,its)
           endif
           xroot=dcmplx(dble(x),0.d0)
      nf=nf + 1
           if(j .lt. m) print *,'root not first: ',m-j+1
c          return
        endif
cc      if(abs(dimag(x)).le.2.*EPS**2*abs(dble(x))) 
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
      subroutine cdhankel(z,tol,H0)
c
c: Computes the zero'th order hankel function of the first kind of the
c: complex argument z using a tolerance of tol.
c
      implicit none
      complex*16 z,H0,zsq,J0,Y0,Y0sum,term,termx
      real*8 tol,gamma,ser,tw_o_pie
      integer*4 j,k,sg
      data gamma/.5772156659015d0/,tw_o_pie/0.63661977236758d0/
c
      zsq=z*z
      k=2 
      ser=1.d0
      term=zsq/4.d0
      Y0sum=term
      sg=1
      J0=1.d0 - term
c
      do j=2,100
         k=k+2
         sg=-sg
         ser=ser + 1.d0/dfloat(j)
         term=term*zsq/dfloat(k*k)
         termx=term*ser
         Y0sum=Y0sum + sg*termx
         J0=J0 - sg*term
         if(cdabs(termx/Y0sum) .lt. tol) then
            Y0=tw_o_pie*((cdlog(0.5d0*z) + gamma)*J0 + Y0sum)
            H0=J0 + dcmplx(-dimag(Y0),dble(Y0))
cc          print *,'j = ',j,z,H0
            return
         endif
      enddo
      print *,'cdhankel failed to converge: ',z,termx,Y0sum
c
      return
      end
ccc
      function magsq(z)
c
      implicit none
      complex*16 z
      real*8 magsq
c
      magsq=dble(z)*dble(z) + dimag(z)*dimag(z)
c
      return
      end
ccc
      function magsq_c8(z)
c
      implicit none
      complex*8 z
      real*4 magsq_c8
c
      magsq_c8=real(z)*real(z) + imag(z)*imag(z)
c
      return
      end
ccc
      function dot(z1,z2)
c
      implicit none
      complex*16 z1,z2
      real*8 dot
      dot=dble(z1)*dble(z2) + dimag(z1)*dimag(z2)
c
      return
      end
ccc
      subroutine hunt(xx,n,x,jlo)
      integer*4 jlo,n
      real*8 x,xx(n)
      integer*4 inc,jhi,jm
      logical ascnd
c
c: EKW FIX (always ascending for us, but tricked when n=1):
cc    ascnd=xx(n) .gt. xx(1)
      ascnd=xx(n) .ge. xx(1)
      if(jlo .le. 0 .or. jlo .gt. n) then
         jlo=0
         jhi=n+1
         goto 3
      endif
      inc=1
      if(x .ge. xx(jlo) .eqv. ascnd) then
1        jhi=jlo+inc
         if(jhi .gt. n)then
            jhi=n+1
         elseif(x .ge. xx(jhi) .eqv. ascnd)then
            jlo=jhi
            inc=inc+inc
            goto 1
         endif
      else
         jhi=jlo
2        jlo=jhi-inc
         if(jlo .lt. 1) then
            jlo=0
         elseif(x .lt. xx(jlo) .eqv. ascnd) then
            jhi=jlo
            inc=inc+inc
            goto 2
         endif
      endif
3     if(jhi-jlo .eq. 1) return
      jm=(jhi+jlo)/2
      if(x .gt. xx(jm) .eqv. ascnd) then
         jlo=jm
      else
         jhi=jm
      endif
      goto 3
c
      end
ccc
      subroutine hunt_re_im(xx_re,n,x,jlo)
      integer*4 jlo,n
      real*8 x,xx_re(n)
      integer*4 inc,jhi,jm
      logical ascnd
c
      ascnd=xx_re(n) .gt. xx_re(1)
      if(jlo .le. 0 .or. jlo .gt. n) then
         jlo=0
         jhi=n+1
         goto 3
      endif
      inc=1
      if(x .ge. xx_re(jlo) .eqv. ascnd) then
1        jhi=jlo+inc
         if(jhi .gt. n)then
            jhi=n+1
         elseif(x .ge. xx_re(jhi) .eqv. ascnd)then
            jlo=jhi
            inc=inc+inc
            goto 1
         endif
      else
         jhi=jlo
2        jlo=jhi-inc
         if(jlo .lt. 1) then
            jlo=0
         elseif(x .lt. xx_re(jlo) .eqv. ascnd) then
            jhi=jlo
            inc=inc+inc
            goto 2
         endif
      endif
3     if(jhi-jlo .eq. 1) return
      jm=(jhi+jlo)/2
      if(x .gt. xx_re(jm) .eqv. ascnd) then
         jlo=jm
      else
         jhi=jm
      endif
      goto 3
c
      end
