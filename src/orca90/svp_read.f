      subroutine svp_read

c	Converted for use of ORCA as a subroutine
c	Variables are either hard-coded herein or
c	set before call to orca
c	Unused lines have been commented out by 'co'
c	Note that several input file options will not be recognized,
c	a stop command will be issued.

      use parms_com
      use i_o_com
c: Local variables:
      integer*4 nline,j,j1,j2,iiblug(-4:4),iierr,ndel
      real*8 zdel(5*NLMAX)
      character*64 eline
      data eline/'INVALID INPUT IN SVP FILE: '/
c: For ktb,ktt<0, read in two more values (p- and s-wave attn freq exponents):
cc    data iiblug/4,4,4,2,0,0,2,2,2/
c: Read in all four even for ii=-1:
c: This not used for version 3.0 and up:
      data iiblug/4,4,4,4,0,0,2,2,2/
c
co      if(iiwrite .ne. 0) write(2,95)
95    format(/'### SVP FILE INFORMATION ###')
co      open(10,file=svp_file(1:lsvp),err=500) 
co      rewind(10)
      iierr=0
      nline=0

c------------------------------------------------
c	Variables now set here 
c------------------------------------------------

c	Line 1
	svp_ver		= 3.0	! should have been 2.0
	svp_title	= '   '

c	Line 2 - Upper halfspace
co        (geot(2,j,j1),j=1,5),(bpt(j,j1),j=3,4)
	geot(2,1,1) = 343.0		! using parameters for air
	geot(2,2,1) = 0.0
	geot(2,3,1) = 0.00121	! density of air
	geot(2,4,1) = 0.0
	geot(2,5,1) = 0.0
	bpt(3:4,1)  = 1.0d0
c	Line 3 - Sound speed in water
co		nsvp,ctol
	ctol = 0
co		zsvp(1),csvp(1),rho_svp,alpha_svp
	rho_svp		= 1.0
	alpha_svp	= 0.0
co        zsvp(j),csvp(j)

c	Line 4 - Bottom layering
co      nlayb
co      ktb(j),hb(j),((geob(j1,j2,j),j1=1,2),j2=1,5),(bpb(j1,j),j1=1,4)
c	Only linear profiles allowed
	do j=1,nlayb
		ktb(j)		= 1
		bpb(1:4,j)	= 1.0d0
	enddo
c	Line 5 - Lower Halfspace
co         (geob(1,j,j1),j=1,5),(bpb(j,j1),j=3,4)
		bpb(3:4,nlayb+1) = 1.0d0
c	Line 6 - Top Layers
co      read(10,*,end=510,err=510) nlayt
	nlayt = 0
c	Line 7 - Top Layer Parameters
co            read(10,*,end=510,err=510) ktt(j),ht(j),((geot(j1,j2,j),
co     .         j1=1,2),j2=1,5),(bpt(j1,j),j1=1,4)
c	------

co      call star2(10,nline,2,svp_file,1,iiwrite)
co      read(10,*,end=510,err=508) svp_ver,svp_title
      call check_val_r4(svp_ver,0.9e0,3.0e0,eline,27,nline,
     .   'SVP version',11,iierr)
490   continue
c
c: Upper halfspace:
      j1=1
co      call star2(10,nline,2,svp_file,1,iiwrite)
      if(svp_ver .gt. 2.e0) then
co         read(10,*,end=510,err=510) (geot(2,j,j1),j=1,5),
co     .      (bpt(j,j1),j=3,4)
      else
         read(10,*,end=510,err=510) (geot(2,j,j1),j=1,5)
         bpt(3:4,j1)=1.d0
      endif
      call check_val_r8(bpt(3,j1),.5d0,3.d0,eline,27,nline,
     .   'p-wave freq exp in upper h-space',35,iierr)
      call check_val_r8(bpt(4,j1),.5d0,3.d0,eline,27,nline,
     .   's-wave freq exp in upper h-space',35,iierr)
      if(geot(2,1,j1) .gt. 0.d0) then
c: c_hsp < 0 means use previous layer as halfspace:
         ht(j1)=1.d+20
         do j=1,5
            geot(1,j,j1)=geot(2,j,j1)
         enddo
         call svp_check_val_lay(2,0.d0,geot(1,1,j1),'upper h-space',13,
     .      eline,nline,iierr)
      endif
c
c: Sound speed profile in ocean:
co      call star2(10,nline,2,svp_file,1,iiwrite)
      if(svp_ver .gt. 0.9d0) then
co         read(10,*,end=510,err=510) nsvp,ctol
         call check_val_r8(ctol,0.d0,50.d0,eline,27,nline,'ctol',
     .      4,iierr)
      else
         read(10,*,end=510,err=510) nsvp
         ctol=0.d0
      endif
      call check_val_i(nsvp,2,5*NLMAX,eline,27,nline,'nsvp',4,iierr)
c
co      call star2(10,nline,2,svp_file,nsvp,iiwrite)
co      read(10,*,end=510,err=510) zsvp(1),csvp(1),rho_svp,alpha_svp
      call check_val_r8(zsvp(1),0.d0,0.d0,eline,27,nline,
     .   'zsvp(1) MUST BE 0',17,iierr)
      call check_val_r8(csvp(1),1.d-10,1.d+10,eline,27,nline,
     .   'csvp(1)',7,iierr)
      call check_val_r8(rho_svp,1.d-100,1.d+100,eline,27,nline,
     .   'rho_svp on first line',21,iierr)
      do j=2,nsvp
co         read(10,*,end=510,err=510) zsvp(j),csvp(j)
         call check_val_r8(zsvp(j),zsvp(j-1),1.d+10,eline,27,nline,
     .      'zsvp(j)',7,iierr)
         call check_val_r8(csvp(j),1.d-10,1.d+10,eline,27,nline,
     .      'csvp(j)',7,iierr)
      enddo
c
c: Check if the number of layers in ocean SVP can be reduced:
      ndel=0
      if(ctol .gt. 0.d0) then
         call svp_fit(nsvp,zsvp,csvp,ctol,ndel,zdel)
      endif
c
c: Bottom layering:
co      call star2(10,nline,2,svp_file,1,iiwrite)
co      read(10,*,end=510,err=510) nlayb
      call check_val_i(nlayb+1,1,NLMAX,eline,27,nline,'nlayb+1',
     .   7,iierr)
co      call star2(10,nline,2,svp_file,nlayb,iiwrite)
      do j=1,nlayb
         if(svp_ver .gt. 2.e0) then
co            read(10,*,end=510,err=510) ktb(j),hb(j),((geob(j1,j2,j),
co     .         j1=1,2),j2=1,5),(bpb(j1,j),j1=1,4)
         else
            read(10,*,end=510,err=510) ktb(j),hb(j),((geob(j1,j2,j),
     .         j1=1,2),j2=1,5),(bpb(j1,j),j1=1,iiblug(ktb(j)))
         endif
         if(j .eq. 1) then
c: Check for negative h, meaning two-way travel time, and negative c,
c: meaning csed/cwater ratio:
            if(geob(1,1,1) .lt. 0.d0) then
               cs_cw_rat=-geob(1,1,1)
               geob(1,1,1)=cs_cw_rat*csvp(nsvp)
			print *,'cb(1) = ',geob(1,1,1)
            endif
            if(hb(1) .lt. 0.d0) then
c: Nominal average gradient (since we don't know it for sure):
               gbar=0.75
               tau_2way=-hb(1)
               hb(1)=geob(1,1,1)*(exp(gbar*tau_2way/2.) - 1.)/gbar
			print *,'hb(1) = ',hb(1)
            endif
         endif
         call check_val_i(ktb(j),-4,4,eline,27,nline,'Profile Type',12,
     .      iierr)
         if(iabs(ktb(j)) .gt. 1) then
            call check_val_r8(bpb(2,j),.1d0,20.d0,eline,27,nline,
     .         'ctol for blug layer',19,iierr)
         endif
         call svp_check_val_lay(1,hb(j),geob(1,1,j),'bottom layer',12,
     .      eline,nline,iierr)
         call svp_check_val_lay(2,hb(j),geob(1,1,j),'bottom layer',12,
     .      eline,nline,iierr)
         call zero_sh(geob(1,2,j),geob(1,5,j),j,'top   ','bottom')
         call zero_sh(geob(2,2,j),geob(2,5,j),j,'bottom','bottom')
      enddo
c
c: Lower halfspace:
      j1=nlayb+1
co      call star2(10,nline,2,svp_file,1,iiwrite)
      if(svp_ver .gt. 2.e0) then
co         read(10,*,end=510,err=510) (geob(1,j,j1),j=1,5),
co     .      (bpb(j,j1),j=3,4)
      else
         read(10,*,end=510,err=510) (geob(1,j,j1),j=1,5)
         bpb(3:4,j1)=1.d0
      endif
      call check_val_r8(bpb(3,j1),.5d0,3.d0,eline,27,nline,
     .   'p-wave freq exp in lower h-space',35,iierr)
      call check_val_r8(bpb(4,j1),.5d0,3.d0,eline,27,nline,
     .   's-wave freq exp in lower h-space',35,iierr)
      if(geob(1,1,j1) .gt. 0.d0) then
c: c_hsp < 0 means use previous layer as halfspace:
c: Set thickness of halfspaces to large numbers (for use in zmx_init):
         hb(j1)=1.d+20
         do j=1,5
            geob(2,j,j1)=geob(1,j,j1)
         enddo
         call svp_check_val_lay(1,0.d0,geob(1,1,j1),'lower h-space',
     .      13,eline,nline,iierr)
         call zero_sh(geob(1,2,j1),geob(1,5,j1),j1,'top   ','bottom')
      endif
c
c: Layering above ocean:
c: NOTE: TOP LAYERS ASSUMED TO BE GIVEN FROM UPPER HALFSPACE TO OCEAN.
co      call star2(10,nline,2,svp_file,1,iiwrite)
co      read(10,*,end=510,err=510) nlayt
      call check_val_i(nlayt+1,1,NLMAX,eline,27,nline,'nlayt+1',
     .   7,iierr)
co      call star2(10,nline,2,svp_file,nlayt,iiwrite)
      do j=2,nlayt+1
         if(svp_ver .gt. 2.e0) then
co            read(10,*,end=510,err=510) ktt(j),ht(j),((geot(j1,j2,j),
co     .         j1=1,2),j2=1,5),(bpt(j1,j),j1=1,4)
         else
            read(10,*,end=510,err=510) ktt(j),ht(j),((geot(j1,j2,j),
     .         j1=1,2),j2=1,5),(bpt(j1,j),j1=1,iiblug(ktt(j)))
         endif
         if(j .eq. 1) then
c: Check for negative h, meaning two-way travel time, and negative c,
c: meaning csed/cwater ratio:
            if(geot(1,1,1) .lt. 0.d0) then
               cs_cw_rat=-geot(1,1,1)
               geot(1,1,1)=cs_cw_rat*csvp(1)
      print *,'ct(1) = ',geot(1,1,1)
            endif
            if(ht(1) .lt. 0.d0) then
c: Nominal average gradient (since we don't know it for sure):
               gbar=0.75
               tau_2way=-ht(1)
               ht(1)=geot(1,1,1)*(exp(gbar*tau_2way/2.) - 1.)/gbar
      print *,'ht(1) = ',ht(1)
            endif
         endif
         call check_val_i(ktt(j),-4,4,eline,27,nline,'Profile Type',12,
     .      iierr)
         if(iabs(ktt(j)) .gt. 1) then
            call check_val_r8(bpt(2,j),.1d0,20.d0,eline,27,nline,
     .         'ctol for blug layer',19,iierr)
         endif
         call svp_check_val_lay(1,ht(j),geot(1,1,j),'top layer',9,
     .      eline,nline,iierr)
         call svp_check_val_lay(2,ht(j),geot(1,1,j),'top layer',9,
     .      eline,nline,iierr)
         call zero_sh(geot(1,2,j),geot(1,5,j),j,'top   ','top   ')
         call zero_sh(geot(2,2,j),geot(2,5,j),j,'bottom','top   ')
      enddo
co      close(10)
c
      nlay=nsvp+nlayb+1+nlayt+1
      call check_val_i(nlay,2,NLMAX,eline,27,nline,'nlay',4,iierr)
c
      if(iierr .ne. 0) then
         print *,'STOPPING BECAUSE ERROR(S) IN INPUT SVP FILE FOUND.'
         stop
      endif
c
      if(ndel .gt. 0 .and. iiwrite .ne. 0) then
co         write(2,410) ndel,ndel+nsvp,(zdel(j),j=1,ndel)
410      format(/'== CTOL ALLOWED DELETION OF',i4,' OUT OF',i4,
     .      ' SVP LAYERS.'/'== DEPTHS DELETED = ',7(f7.2,1x)/
     .      40(3x,9(1x,f7.2)/))
      endif
c
      RETURN

500   print *,'Error opening SVP file ',svp_file
      stop
c: Old SVP format with version number:
508   svp_ver=0.9
      backspace(10)
      read(10,200,end=510,err=510) svp_title
200   format(a)
      goto 490
510   print *,'Endo or error reading SVP file ',svp_file,
     .   ' at line ',nline
      stop
      end
ccc
      subroutine zero_sh(cs,as,j,ch_tb1,ch_tb2)
c
      implicit none
      real*8 cs,as
      integer*4 j
      character*6 ch_tb1,ch_tb2
c
      if(cs .eq. 0.d0 .and. as .ne. 0.d0) then
         write(2,125) ch_tb1,ch_tb2,j
125      format('/*** SHEAR ATTENUATION SET TO 0 for 0 shear speed at ',
     .         a6,' of ',a6,' layer # ',i3,'***'/)
         as=0.d0
      endif
c
      return
      end
ccc
      subroutine svp_fit(nsvp,zsvp,csvp,ctol,ndel,zdel)
c
      implicit none
      integer*4 nsvp,j1,j2,j,jj,ndel,nd
      real*8 zsvp(nsvp),csvp(nsvp),ctol,zdel(nsvp),c1,c1_2,c2,cfit,gfac
c
      j1=1
      c1=csvp(j1)
      c1_2=1/(c1*c1)
15    continue
      j2=j1 + 2
      if(j2 .gt. nsvp) return
20    c2=csvp(j2)
      if(c1 .ne. c2) then
         gfac=(1/(c2*c2) - c1_2)/(zsvp(j2) - zsvp(j1))
      else
         cfit=c1
      endif
      do j=j1+1,j2-1
         if(c1 .ne. c2) then
            cfit=1.d0/sqrt(c1_2 + gfac*(zsvp(j) - zsvp(j1)))
         endif
c: Check fit:
         if(abs(cfit-csvp(j)) .gt. ctol) then
c: Fit bad, so check if any points can be deleted:
            nd=j2-j1-2
            if(nd .gt. 0) then
               do jj=j1+1,j2-2
                  ndel=ndel + 1
                  zdel(ndel)=zsvp(jj)
               enddo
               do jj=j1+1,nsvp-nd
                  csvp(jj)=csvp(jj+nd)
                  zsvp(jj)=zsvp(jj+nd)
               enddo
               nsvp=nsvp - nd
               j2=j2 - nd
            endif
            if(j2 .eq. nsvp) return
            j1=j2-1
            c1=csvp(j1)
            c1_2=1/(c1*c1)
            goto 15
         endif
      enddo
      if(j2 .lt. nsvp) then
         j2=j2 + 1
         goto 20
      else
         nd=j2-j1-1
         do jj=j1+1,j2-1
            ndel=ndel + 1
            zdel(ndel)=zsvp(jj)
         enddo
         do jj=j1+1,nsvp-nd
            csvp(jj)=csvp(jj+nd)
            zsvp(jj)=zsvp(jj+nd)
         enddo
         nsvp=nsvp - nd
         return
      endif
c
      end
ccc
      subroutine svp_check_val_lay(ii,h,geo,char_lay,nch,eline,nline,
     .   iierr)
c
      implicit none
      real*8 h,geo(2,5)
      integer*4 ii,nch,nline,iierr
      character*64 char_lay,eline
c
c: Make negative h mean two-way travel time:
      call check_val_r8(h,-10.0d0,1.5d4,eline,27,nline,
     .   'h '//char_lay(1:nch),2+nch,iierr)
c: Make negative c mean sound speed ratio:
      call check_val_r8(geo(ii,1),-10.d0,1.d+10,eline,27,nline,
     .   'cp '//char_lay(1:nch),3+nch,iierr)
      call check_val_r8(geo(ii,2),0.0d0,1.d+10,eline,27,nline,
     .   'cs '//char_lay(1:nch),3+nch,iierr)
      call check_val_r8(geo(ii,3),1.d-50,1.d+50,eline,27,nline,
     .   'rho '//char_lay(1:nch),4+nch,iierr)
      call check_val_r8(geo(ii,4),-200.d0,999.d0,eline,27,nline,
     .   'ap '//char_lay(1:nch),3+nch,iierr)
      call check_val_r8(geo(ii,5),-200.d0,200.d0,eline,27,nline,
     .   'as '//char_lay(1:nch),3+nch,iierr)
c
      return
      end
