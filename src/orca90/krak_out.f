      subroutine krak_out
c
c: Writes out a Kraken .env file corresponding to the current ORCA run.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 j,nmesh,mesh,meshfac,nmedia
      real*8 hwat,cphhigh,f_cw,f_rat
c
      if(iicw .eq. 1) then
         f_cw=fcw(1)
      elseif(iicw .eq. 2) then
         f_cw=.5*(fminorca + fmaxorca)
      elseif(iikpl .gt. 0) then
         f_cw=fkpl
      elseif(iirc .gt. 0) then
         f_cw=.5*(freq1 + freq2)
      else
         f_cw=fcw(1)
      endif
      if(f_cw .eq. 0.d0) then
         print *,'f_cw set to 50 since f_cw=0'
         f_cw=50.d0
      endif
      f_rat=f_cw/1000.
c
      if(iikrak .gt. 1) then
         meshfac=max(10,iikrak)
      else
         meshfac=20
      endif
c
      open(19,file=outroot(1:lout)//'_kr.env',form='formatted')
      write(19,100) svp_title
100   format(' '' ',a,' '' ')
      write(19,102) f_cw
102   format(f9.2,'    ! FREQUENCY')
cc    write(19,104) nlayb + nlayt + 1
      nmedia=1 + (jsurf-2) + (nlay-jobot-1)
      if(h(jsol(1,1)-1) .eq. 0.d0) nmedia=nmedia-1
      if(h(jsol(2,1)+1) .eq. 0.d0) nmedia=nmedia-1
      write(19,104) nmedia
104   format(i9,  '    ! NMEDIA')
      write(19,105) 
105   format(' ''NVW''       ! OPTIONS')
c
      do j=2,jsurf-1
         if(h(j) .ne. 0.d0) then
            nmesh=mesh(geo(1,1,j),h(j),f_cw,meshfac)
            write(19,110) nmesh,0.,zdep(j)
110         format(i5,2x,f5.2,f8.2,'    ! NMESH,SIGMA,LAST DEPTH') 
            call krak_wr(zdep,geo(1,1,j),j,1,fexp(1,j),f_rat)
            call krak_wr(zdep,geo(1,1,j),j,2,fexp(1,j),f_rat)
         endif
      enddo
c
c: SVP in ocean:
      hwat=zdep(jobot) - zdep(jsurf-1)
      nmesh=mesh(geo(1,1,jsurf),hwat,f_cw,meshfac)
      write(19,110) nmesh,0.,zdep(jobot)
      j=jsurf
      call krak_wr(zdep,geo(1,1,j),j,1,fexp(1,j),f_rat)
      do j=jsurf,jobot
         call krak_wr(zdep,geo(1,1,j),j,2,fexp(1,j),f_rat)
      enddo
c
c: Bottom layers:
      do j=jobot+1,nlay-1
         if(h(j) .ne. 0.d0) then
            nmesh=mesh(geo(1,1,j),h(j),f_cw,meshfac)
            write(19,110) nmesh,0.,zdep(j)
            call krak_wr(zdep,geo(1,1,j),j,1,fexp(1,j),f_rat)
            call krak_wr(zdep,geo(1,1,j),j,2,fexp(1,j),f_rat)
         endif
      enddo
c
c: Acoustic halfspace with no roughness:
      write(19,115)
115   format(' ''A'' 0.0    ! ACOUSTIC H-SPACE, NO ROUGHNESS')
      j=nlay
      call krak_wr(zdep,geo(1,1,j),j,1,fexp(1,j),f_rat)
      if(cphmax .lt. 1.e10 .and. cphmax .gt. 0.d0) then
         cphhigh=cphmax
      else
         cphhigh=w/dreal(kn(nmode)) + 10.d0
      endif
      write(19,120) cphlo,cphhigh,rmax
120   format(f9.2,2x,g13.6,'    ! CPHMIN,CPHMAX'/f7.1,'     ! RMAX',
     .   ' (NZS,ZS; NZR,ZREC ON NEXT 2 LINES:)')
      if(iitl .ne. 0 .or. iifft .ne. 0) then
         write(19,130) nzs,(zsrc(j),j=1,nzs)
         write(19,130) nrec,(zrec(j),j=1,nrec)
      elseif(iimf .ne. 0 .or. iisig .ne. 0) then
         write(19,130) 1,zmf(1)
         write(19,130) nzmf,(zmf(j),j=1,nzmf)
      else
         print *,'Note: No depths requested for ORCA.  Outputting '//
     .      'no depths for Kraken.'
         write(19,130) 1,10.
         write(19,130) 1,10.
      endif
130   format(i4,2x,2000(f8.2))
c
      close(19)
c
c: Write out a .FLP file also:
      open(19,file=outroot(1:lout)//'_kr.flp',form='formatted')
      write(19,150) 
150   format('/,              !TITLE'/
     .       ' ''RA''           ! OPT ''X/R'' ''C/A'' '/
     .       '9999            ! M (number of modes to include)'/
     .       '1  0.0          ! NPROF  RPROF(1:NPROF) (km)')
      if(nrng .gt. 4) then
         if(abs(range(nrng)-range(nrng-1) - (range(2)-range(1))) .lt.
     .      .001) then
            write(19,152) nrng,.001*range(1),.001*range(nrng),'/'
152         format(i4,2(1x,f8.3),1x,a1)
         else
            write(19,153) nrng,(.001*range(j),j=1,nrng)
         endif
      else
         write(19,153) nrng,(.001*range(j),j=1,nrng)
153      format(i4,1x,2000(f8.3,1x))
      endif
      write(19,130) nzs,(zsrc(j),j=1,nzs)
      write(19,130) nrec,(zrec(j),j=1,nrec)
      write(19,130) nrec,(0,j=1,nrec)
      close(19)
c
      return
      end
ccc
      subroutine krak_wr(zdep,geo,j,ii,fexp,f_rat)
c
      implicit none
      integer*4 j,ii,jj
      real*8 zdep(j),geo(2,5),fexp(2),f_rat,alpha_p,alpha_s
c
      alpha_p=geo(ii,4)*(f_rat**(fexp(1)-1))
      alpha_s=geo(ii,5)*(f_rat**(fexp(2)-1))
      write(19,120) zdep(j+ii-2),(geo(ii,jj),jj=1,3),
     .   .001*geo(ii,1)*alpha_p,.001*geo(ii,2)*alpha_s
120   format(f8.2,2x,f7.2,2x,f7.2,2x,f7.3,2(2x,f8.4))
c
      return
      end
ccc
      function mesh(geo,h,f_hz,meshfac)
c
      implicit none
      integer*4 mesh,meshfac
      real*8 geo(2,5),h,f_hz,lambda
c
      if(geo(1,2) .eq. 0.d0 .and. geo(2,2) .eq. 0.d0) then
         lambda=min(geo(1,1),geo(2,1))/f_hz
      else
         lambda=min(geo(1,2),geo(2,2))/f_hz
      endif
      mesh=max(1,nint(float(meshfac)*h/lambda))
c
      return
      end
ccc
      subroutine oases_out
c
c: Writes out an OASES .dat file corresponding to the current ORCA run.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 j,nw,nloas,ic2
      real*8 xnw,cphlox,cphlow,cphhigh,cph_fac,fmx,dr_orca,dr_oas,
     .   dc_inv,f_rat
c
      open(19,file=outroot(1:lout)//'_oas.dat',form='formatted')
      if(iicw .eq. 1) then
         write(19,100) svp_title,'N','J'
         write(19,102) fcw(1),fcw(nfcw),nfcw,0.0
         fmx=fcw(nfcw)
      elseif(iicw .eq. 2) then
         write(19,100) svp_title,'N','J','F','1'
         write(19,102) (fminorca+fmaxorca)/2.,(fminorca+fmaxorca)/2.,1,0.0
         fmx=fmaxorca
      endif
      f_rat=fmx/1000.
100   format(a/20(a1,1x))
102   format(f9.2,1x,f9.2,1x,i4,1x,i4,'    ! F1,F2,NF,COFF')
c
      nloas=0
      do j=1,nlay
         if(h(j) .ne. 0.d0) then
cc          if(isp(j) .eq. 0 .and. (geo(1,1,j) .eq. geo(2,1,j) .or.
cc   .         geo(1,2,j) .eq. 0.d0)) then
cc             nloas=nloas + 2
cc          else
               nloas=nloas + 1
cc          endif
         endif
      enddo
      write(19,104) nloas
104   format(i9,21x,'    ! NLAY')
c
c: Profile:
      do j=1,nlay
         if(h(j) .ne. 0.d0) then
            call oas_wr(zdep,geo(1,1,j),j,isp(j),fexp(1,j),f_rat)
            nloas=nloas + 1
         endif
      enddo
c
      write(19,115) zsrc(1)
      write(19,115) zrec(1),zrec(nrec),nrec,1
115   format(2(f9.3,1x),i4,1x,i4)
cc    cphhigh=1.5d0*w/dreal(kn(nmode))
      cphhigh=1.e8
      dr_orca=(range(nrng)-range(1))/max(1,nrng-1)
      cphlox=0.8d0*cphlo
      if(dr_orca .gt. 0.d0) then
         cphlow=min(cphlox,fmx*dr_orca)
      else
         cphlow=cphlox
      endif
      dc_inv=1.d0/cphlow - 1.d0/cphhigh
      if(iitl .ne. 0) then
         xnw=range(nrng)*fmx*dc_inv
      else
         xnw=10000.*fmx*dc_inv
      endif
c: Add 2 instead of 1 because adding 1 usually results in jitter at large R:
cxx   nw=2**(int(log10(xnw)/log10(2.)) + 2)
c: Set iioas=1 for default power of 2, iioas=2 for larger value:
      nw=2**(int(log10(xnw)/log10(2.)) + 1 + iioas)
c: Now fine tune cmin so that dr_oas=dr_orca or an integral fraction of it:
      if(dr_orca .gt. 0.d0) then
         dr_oas=(nw-1)/(nw*fmx*dc_inv)
         if(dr_oas .lt. .95*dr_orca) then
            dr_oas=dr_orca/(int(dr_orca/dr_oas)+1)
         else
            dr_oas=dr_orca
         endif
         cph_fac=(nw-1)/(nw*fmx*dr_oas) + 1.d0/cphhigh
         cphlow=1.d0/cph_fac
      endif
      ic2=nw
      if(cphlow .lt. cphlox) then
         dc_inv=1.d0/cphlow - 1.d0/cphhigh
         ic2=nint(nw*(1.d0/cphlox - 1.d0/cphhigh)/dc_inv)
      endif
cc    if(iicw .eq. 2) nw=-1
      write(19,120) cphlow,cphhigh
120   format(f9.2,2x,g13.6,6x,'    ! CPHMIN,CPHMAX')
      write(19,122) nw,1,ic2,0
122   format(4(i6,1x),'      ! NW(-1=auto),IC1,IC2,KERNEL PLOT ,',
     .   'INTERVAL')
c
      dr_orca=1.
      if(nrng .gt. 1) dr_orca=rkm(2) - rkm(1)
      if(iicw .eq. 1) then
cc         write(19,124) rkm(1),dr_orca,nrng,1
cc124      format(2(f9.4,1x),i5,1x,i2,9x,'! RMIN,DR,NR,DUMMY')
         write(19,124) rkm(1),rkm(nrng),nrng,1
124      format(2(f9.4,1x),i5,1x,i2,9x,'! RMIN, RMAX, DUMMY, DUMMY')
      else
       write(19,126) nfftbb,fminorca,fmaxorca,
     .      1./fsbb,rkm(1),dr_orca,nrng,1
126      format(i6,1x,f9.2,1x,f9.2,1x,f8.6,1x,
     .      2(f9.4,1x),i2,1x,i2,9x,'! NFFT,FMIN,FMAX,DT,',
     .      'RMIN,DR,NR,DUMMY')
      endif
c
      close(19)
c
      return
      end
ccc
      subroutine oas_wr(zdep,geo,j,isp,fexp,f_rat)
c
      implicit none
      integer*4 j,isp
      real*8 zdep(j),geo(2,5),fexp(2),f_rat,alpha_p,alpha_s
c
      alpha_p=geo(1,4)*(f_rat**(fexp(1)-1))
      alpha_s=geo(1,5)*(f_rat**(fexp(2)-1))
      if(isp .eq. 1) then
c: Isospeed fluid layer:
         write(19,120) zdep(max(1,j-1)),geo(1,1),geo(1,2),
     .      .001*geo(1,1)*alpha_p,.001*geo(1,2)*alpha_s,
     .      .5d0*(geo(1,3)+geo(2,3)),0.0
      else
         if(geo(1,1) .ne. geo(2,1) .and. geo(1,2) .eq. 0.d0) then
c: Airy layer in fluid:
            write(19,120) zdep(max(1,j-1)),geo(1,1),-geo(2,1),
     .         .001*geo(1,1)*alpha_p,.001*geo(1,2)*alpha_s,
     .         .5d0*(geo(1,3)+geo(2,3)),0.0
         else
c: Top of layer:
            write(19,120) zdep(max(1,j-1)),geo(1,1),geo(1,2),
     .         .001*geo(1,1)*alpha_p,.001*geo(1,2)*alpha_s,
     .         .5d0*(geo(1,3)+geo(2,3)),0.0
cc   .         geo(1,3),0.0
c: Bottom of layer:
cc          write(19,120) zdep(j),geo(2,1),geo(2,2),
cc   .         .001*geo(2,1)*geo(2,4),.001*geo(2,2)*geo(2,5),
cc   .         geo(2,3),0.0
         endif
      endif
120   format(f8.2,2(2x,f9.2),2(2x,f8.4),2x,f7.3,2x,f3.1)
c
      return
      end
ccc
      subroutine fepe_out
c
c: Writes out a FEPE fepe.in or ram.in file corresponding to the 
c: current ORCA run.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 j,ndr_fepe,ndz_fepe,npade,nlam_bot,nlam_false
      real*8 hwat,lambda,dr_fepe,dz_fepe,zmax_fepe,attn_bot,
     .   dr_orca,dz_orca,f_cw,f_rat,alpha_p1,alpha_p2
c
cc    nlam_bot=18
      nlam_bot=10
      nlam_false=3
      if(iicw .eq. 1) then
         f_cw=f_hz
      elseif(iicw .eq. 2) then
         f_cw=.5*(fminorca + fmaxorca)
      elseif(iikpl .gt. 0) then
         f_cw=fkpl
      elseif(iirc .gt. 0) then
         f_cw=.5*(freq1 + freq2)
      else
         f_cw=fcw(1)
      endif
      if(f_cw .eq. 0.d0) then
         print *,'f_cw set to 50 since f_cw=0'
         f_cw=50.d0
      endif
c
      if(iifepe .eq. 1) then
         open(19,file=outroot(1:lout)//'_fepe.in',form='formatted')
      else
         open(19,file=outroot(1:lout)//'_ram.in',form='formatted')
      endif
      write(19,100) svp_title
c100  format(' '' ',a,' '' ')
100   format(a)
      write(19,102) f_cw,zsr(mzsrc(1)),zsr(mzrec(1))
102   format(f9.2,2x,f9.2,2x,f9.2,11x,'    freq zs zr')
      lambda=cfmin/f_cw
      if(iifepe .eq. 1) then
         dr_fepe=lambda/12.
cc       dz_fepe=lambda/60.
         dz_fepe=dr_fepe/5.
      else
         dr_fepe=lambda/2.
         dz_fepe=dr_fepe/5.
      endif
      if(nrng .gt. 1) then
         dr_orca=(range(nrng) - range(1))/max(1,nrng-1) 
         if(dr_orca .eq. 0) dr_orca=100.
      else
         dr_orca=100.
      endif
      ndr_fepe=max(1,nint(dr_orca/dr_fepe))
      dr_fepe=dr_orca/ndr_fepe
      write(19,104) range(nrng)+dr_fepe/2,dr_fepe,ndr_fepe
104   format(f9.2,2x,f9.4,2x,i9,11x,'    rmax dr(!) ndr')
      if(nrec .gt. 1) then
         dz_orca=(zsr(mzrec(nrec)) - zsr(mzrec(1)))/(nrec-1)
      else
         dz_orca=5.
      endif
      ndz_fepe=max(1,nint(dz_orca/dz_fepe))
      dz_fepe=dz_orca/ndz_fepe
      zmax_fepe=zdep(nlay-1) + nlam_bot*lambda
      write(19,105) zmax_fepe,dz_fepe,ndz_fepe,zsr(nzsr)+dz_fepe/2
105   format(f9.2,2x,f9.4,2x,i9,2x,f9.2,'    zmax(!) dz(!) ndz zmplt')
      if(iifepe .eq. 1) then
         npade=2
         write(19,106) cfmin,npade,70.,0
106      format(f9.2,2x,i9,2x,f9.2,2x,i9,'    c0 npade(!) theta(!) ',
     .      'irefl')
         write(19,116) 5,2.*lambda,89.9,npade+1
116      format(i9,2x,f9.2,2x,f9.2,2x,i9,'    istrt rmin(5!) ',
     .      'thmax(3) mpade(5!)')
      else
         npade=5
         write(19,206) cfmin,npade,1,0.
206      format(f9.2,2x,i9,2x,i9,2x,f9.2,'    c0 npade(!) ns rs')
      endif
      hwat=zdep(jobot) - zdep(jsurf-1)
      write(19,107) 0.,hwat
107   format(f9.2,2x,f9.2,22x,'    bathymetry: r(m),z(m)'/'-1 -1')
      write(19,108) zdep(jsurf-1),geo(1,1,jsurf)
      do j=jsurf,jobot
         write(19,108) zdep(j),geo(2,1,j)
108      format(f9.2,2x,f9.2,22x,'    z cw')
      enddo
      write(19,109)
109   format('-1 -1')
c
c: cp in bottom layers:
      do j=jobot+1,nlay-1
         if(h(j) .ne. 0.d0) then
            write(19,120) zdep(j-1),geo(1,1,j)
            write(19,120) zdep(j),geo(2,1,j)
120         format(f9.2,2x,f9.2,22x,'    zb cb')
         endif
      enddo
      write(19,120) zdep(nlay-1),geo(1,1,nlay)
      write(19,109)
c
c: rho in bottom layers:
      do j=jobot+1,nlay-1
         if(h(j) .ne. 0.d0) then
            write(19,121) zdep(j-1),geo(1,3,j)
            write(19,121) zdep(j),geo(2,3,j)
121         format(f9.2,2x,f9.2,22x,'    zb rhob')
         endif
      enddo
      write(19,121) zdep(nlay-1),geo(1,3,nlay)
      write(19,109)
c
      f_rat=f_cw/1000.
c: attn in bottom layers:
      do j=jobot+1,nlay-1
         if(h(j) .ne. 0.d0) then
            alpha_p1=geo(1,4,j)*(f_rat**(fexp(1,j)-1))
            alpha_p2=geo(2,4,j)*(f_rat**(fexp(1,j)-1))
            write(19,122) zdep(j-1),.001*geo(1,1,j)*alpha_p1
            write(19,122) zdep(j),.001*geo(2,1,j)*alpha_p2
122         format(f9.2,2x,f9.2,22x,'    zb attnb')
         endif
      enddo
      j=nlay
      alpha_p1=geo(1,4,j)*(f_rat**(fexp(1,j)-1))
      attn_bot=.001*geo(1,1,j)*alpha_p1
      write(19,122) zdep(j-1),attn_bot
      write(19,123) zmax_fepe-nlam_false*lambda,attn_bot
      write(19,123) zmax_fepe,max(attn_bot,10.d0)
123   format(f9.2,2x,f9.2,22x,'    zb(!) attnb(!) (False Bottom)')
      write(19,109)
c
      close(19)
c
      return
      end
ccc
      subroutine modelab_out
c
c: Writes out a MODELAB .svp file corresponding to the current ORCA run.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 j
      real*8 f_cw,f_rat,alpha_p1,alpha_p2
c
      if(iicw .eq. 1) then
         f_cw=fcw(1)
      elseif(iicw .eq. 2) then
         f_cw=.5*(fminorca + fmaxorca)
      elseif(iikpl .gt. 0) then
         f_cw=fkpl
      elseif(iirc .gt. 0) then
         f_cw=.5*(freq1 + freq2)
      endif
      f_rat=f_cw/1000.
c
      open(19,file=outroot(1:lout)//'_mlab.svp',form='formatted')
      write(19,100) f_cw,geo(1,1,nlay),40,1.e-12
100   format(
     .   'INPUT FILE FOR MODELAB (program reads input after ',
     .   'starred lines)'/
     .   'f      = frequency'/
     .   'vphmax = max phase velocity'/
     .   'wgdf   = sampling factor (20 for single duct, ',
     .   '40 for multi-duct)'/
     .   'ktol   = error in k allowed (recommend 1.e-12)'/
     .   '*        f      vphmax  wgdf      ktol'/
     .   f10.4,2x,f10.3,2x,i4,2x,e8.1)
c
      write(19,110)
110   format('-----------------------------------------------',
     .   '------------------------'/
     .   '*       z           c  alpha(dB/m/kHz)          rho')
c: Ocean layers:
      do j=jsurf,jobot
         alpha_p1=geo(1,4,j)*(f_rat**(fexp(1,j)-1))
         alpha_p2=geo(2,4,j)*(f_rat**(fexp(1,j)-1))
         write(19,115) zdep(j-1),geo(1,1,j),alpha_p1,geo(1,3,j)
      enddo
c: Bottom of last ocean layer:
      j=jobot
      alpha_p2=geo(2,4,j)*(f_rat**(fexp(1,j)-1))
      write(19,115) zdep(j),geo(2,1,j),alpha_p2,geo(2,3,j)
c: Bottom layers:
      do j=jobot+1,nlay-1
         if(h(j) .ne. 0.d0) then
            alpha_p1=geo(1,4,j)*(f_rat**(fexp(1,j)-1))
            alpha_p2=geo(2,4,j)*(f_rat**(fexp(1,j)-1))
            write(19,115) zdep(j-1),geo(1,1,j),alpha_p1,geo(1,3,j)
            write(19,115) zdep(j),geo(2,1,j),alpha_p2,geo(2,3,j)
         endif
      enddo
c: Top of lower halfspace:
      j=nlay
      alpha_p1=geo(1,4,j)*(f_rat**(fexp(1,j)-1))
      write(19,115) zdep(j-1),geo(1,1,j),alpha_p1,geo(1,3,j)
115   format(f9.3,2x,f10.3,2x,f15.5,2x,e11.5)
c
      write(19,120)
120   format('--------------------------------------------------',
     .   '---------------------'/
     .   'Depths at which to compute mode functions'/
     .   'List depths at which to compute mode functions as positive ',
     .   'numbers.'/
     .   'In addition you may enter triplets of: -dz z1 z2'/
     .   '* z1 z2 ...  AND/OR  -dz z1 z2  (pos=indiv depth, neg=dz ',
     .   'from z1 to z2)')
      write(19,125) (zsr(j),j=1,nzsr) 
125   format(2000(f9.2,1x))
c
      close(19)
c
      return
      end
