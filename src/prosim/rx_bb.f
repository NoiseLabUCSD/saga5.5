      subroutine rx_bb(nsctr,jfhi,jflo,intrpl)
c
c: Does broadband calculations for modes on the real axis.
c
      implicit none
      include 'Parms_com'
      include 'i_o_opt_com'
      common /eig1_com/ eig_char(5,NM_MAX)
      complex*16 eig_char
      include 'i_o_1c_com'
      include 'i_o_2_com'
      include 'gen1_com'
      include 'gen3_com'
      include 'deltaf_com'
      common /ny_com/ narray(NVRMAX),act_sect
      integer*4 narray,act_sect
      common /jmhilo_com/ jmhi_old(RGMAX),jmlo_old(RGMAX)
      integer*4 jmhi_old, jmlo_old
c
      real*8 kmax,kmin,Dfn,logmode,pie
      real*4 f1,f2,fint_max,nfj,dfj
c pln: correction from version 1.1 done by EKW
      integer*4 ndv,nm1,nm2,jflo,jfhi,ndfmax,nm_fmax,
     .          ncx,iibad,ndf,j,j1,kdf,
     .          jmlo,jmhi,iick,isum,nsctr,intrpl(RGMAX)
      logical adpt
c
       pie=3.14159265358979d0
c
       act_sect=nsctr
c Set adpt=.true. for the adaptive Df
       adpt=.false.
c       adpt=.true.
c
c: Check input variables, sizes of arrays, and initialize variables:
       iick=0
       if(iidiag .eq. -10) iick=1
       ndv=6
       isum=nsctr+isubb
       call bb_init(isum)
       call rx_prep
       call rx_zmx_init(nrecusr)
c
       if(multcw) then
         ndfmax=2
         ndfmaxl=2
         adpt=.false.
       else
         if(isum .eq. 2) then
c: Make max frequency jump close to fint_max:
           fint_max=20.0
cpln           fint_max=15.0
cpln           fint_max=12.0
           if(adpt) fint_max=70.0
           fint_maxl=fint_max
c fbv           ndfmaxl=min(nfbb,nint(fint_maxl/df)+1)
cc         ndfmax=max(1,int((fmax-fmin)/fint_max+0.9999)) !# of subbands
         endif
       endif
c
c: Find frequencies limits in block:
       jfhi=nfbb
       if(isubb .ne. 1) then
c fbvfbv        nm_lim = nmode_init
        jfhi=jfhi_up
       endif
       f2=faxbb(jfhi)
c
c: Find (or recall) modes at higher frequency in block:
       f_hz=f2
       nctot=0
       nzoff=0
       call rx_freq_init(kmax,kmin)
c
       if(adpt) then
          if (isum .eq. 2) then
           call rx_mode_int(kmax,kmin,nm_lim,ndv,0,1,1)
cpln           call rx_mode_int(kmax,kmin,nm_lim,ndv,0,0,1)
           nfj=((nmode-2)*crmax)/300.0
           dfj=f2-nfj
           if(nmode .lt. 400) then
              logmode = 2.*log(float(nmode))
              Dfn=log10(dfj/logmode)
c fbv              Dfn=log10(f2/logmode)
           else
              logmode = log10(float(nmode))
              Dfn=log(dfj/logmode)
           endif
           fint_maxl = min(fint_maxl*nint(Dfn),NFCMAX)
          else
           if(nsctr .eq. 1) then
              fint_maxl = max(int(fint_maxl-(0.60*log10(f2))),30)
           endif
          endif
       endif
c     
       if(.not. multcw) then
c fbv         ndfmax=min(nfbb,nint(fint_maxl/df)+1)
         ndfmax=min(NFCMAX,nint(fint_maxl/df)+1)
c fbv         ndfmax=min(nfbb,nint(fint_max/df)+1)
c
         call var_lim(ndfmax,3,VLINE,48,
     .             'NDFMAX',8,iibad,0)
       endif
c
c fbv       if(ndfmax .gt. 200) then
c fbv        print *,'Increase array size in Parms_com: ',200,ndfmaxl
c fbv        stop
c fbv       endif
       jflo = max(jfhi - ndfmax+1, 1)
c
c: Check memory limits:

       if(multcw) then
         call mem_lim(nfbb,NFBBMAX,MLINE,LML,
     .             'NFBBMAX',26,'NFBBMAX',11,iibad,1)
       endif
c
       ndf=jfhi - jflo + 1
       ncx=0
       nmode=0
c If there are problems with memory uncomment the 
c following lines with 'c MEMORY'
c: Rcall modes at lower frequency in block:
       if(isubb .eq. 1) then
c pln: correction from version 1.1 done by EKW
          kim_min=1.d100
          call rx_mode_int(kmax,kmin,nm_lim,ndv,0,1,ndf)
cpln          call rx_mode_int(kmax,kmin,nm_lim,ndv,0,0,ndf)
c fbvfbv         nmode_init = nmode
          if(nmode .eq. 0) then
             write(6,*)'N. of modes = 0'
             iifail=1
             return
          end if
       else 
          call rx_mode_rcl(nsctr,nmode,nzsr,nctot,eig_char,
     .         dphiz_bb,ddphiz_bb,phi_bb,dphi_bb,kn,
     .         NVRMAX,NM_MAX,NFCMAX,ndf)
       endif
c: Check memory limits:

       if(isubb .eq. 1) then
         call mem_lim(nmode,NM_MAX,MLINE,LML,
     .             'NM_MAX',26,'NM_MAX',11,iibad,1)
       endif
c
       if (iiwrite.ge.2) then
          write(6,*)'Min range: ',rmin
          write(6,*)'Max range: ',rmax
          write(6,*)'Max no. of modes: ',nmode
       end if
       nm_fmax=nmode
c: Enter results in higher f of subband:
c       write(6,*)'inside rx_bb #1'
       call rx_bb_enter(1,nmode,ndf,nzsr,
     .                  kn,eig_char,knbb,vg_bb,dvg_bb,
     .                  NVRMAX,NM_MAX,NFCMAX)
c
c: Find modes at lower frequency in block:
       f1=faxbb(jflo)
       f_hz=f1
       call rx_freq_init(kmax,kmin)
       ncx=nctot
       nmode=0
c pln: correction from version 1.1 done by EKW
       kim_min=1.d100
       nm2=nm_fmax
       call rx_mode_int(kmax,kmin,nm2,ndv,0,1,1)
cpln       call rx_mode_int(kmax,kmin,nm2,ndv,0,0,1)
       if (iiwrite.ge.2)
     .    write(6,*)'Max no. of modes: ',nmode
       if(nmode .eq. 0) then
         write(6,*)'N. of modes = 0'
         iifail=1
         return
c         stop
       endif
c
c If there are problems with memory uncomment the 
c following lines with 'c MEMORY'
c: Save modes at lower frequency in block:
       call rx_mode_save(nsctr,nmode,nzsr,ncx,eig_char,
     .                   dphiz_bb,ddphiz_bb,phi_bb,dphi_bb,kn,
     .                   NVRMAX,NM_MAX,NFCMAX,1)

c: Enter results in lower f of subband:
c: write(6,*)'inside rx_bb #2'
       call rx_bb_enter(1,nmode,1,nzsr,
     .                  kn,eig_char,knbb,vg_bb,dvg_bb,
     .                  NVRMAX,NM_MAX,NFCMAX)
c
c: Perform interpolation in frequency:
       nm1=nmode

c: (Re)initialization of jmhi & jmlo
       if(isubb .gt. 1) then
          jmhi=jmhi_old(nsctr)
          jmlo=jmlo_old(nsctr)
       endif

       call rx_mcross(knbb,nm1,nm2,ndf,jmlo,jmhi,
     .                ndrx,f1,f2,ccr_lo,ccr_hi,pie,
     .                NM_MAX,NFCMAX)
c
       if(.not. multcw) then
c
        call rx_bb_mterp(nm1,nm2,jflo,jfhi,jmlo,jmhi,knbb,
     .                   vg_bb,dvg_bb,phi_bb,dphi_bb,dphiz_bb,
     .                   ddphiz_bb,ndv,NVRMAX,NM_MAX,NFCMAX)
c
cpln added 17-3-97
        if(iifail.gt.0) then
           if (iiwrite.ge.2) then
              write(6,*)
              write(6,*)'Brute force for frequency band:'
              write(6,*)jflo,jfhi,faxbb(jflo),faxbb(jfhi)
           end if
          iifail=0
          do j=jfhi,jflo,-1
            kdf=j-jflo+1
            f_hz=faxbb(j)
            call rx_freq_init(kmax,kmin)
            ncx=nctot
            nmode=0
            call rx_mode_int(kmax,kmin,nm_lim,ndv,0,1,kdf)
cpln            call rx_mode_int(kmax,kmin,nm_lim,ndv,0,0,kdf)
            if(iifail.gt.0) return
            if(j.eq.jfhi) then
               j1=nfbb*nmode/2 + 2
               nm_fmax=nmode
            end if
            nm1=nmode
            call rx_bb_enter(1,nmode,kdf,nzsr,
     .                       kn,eig_char,knbb,vg_bb,dvg_bb,
     .                       NVRMAX,NM_MAX,NFCMAX)
          end do
         end if
c
         if(iick .eq. 1) then
           call interp_check(phi_bb,jflo,jfhi,
     .                       NVRMAX,NM_MAX,NFCMAX)
        endif
       endif

c: Perform interpolation in range & evaluate the pressure:
       nmode=nm_fmax
c
       call rx_mode_field(jfhi,jflo,nsctr,intrpl,knbb,
     .                      phi_bb,nrec,nrecusr)
       if(nzoff .gt. 0) then
         if(iiwrite.ge.2)
     .    print *,'Info msg: # of zeros of mode function off by one ',
     .      nzoff,' times.  Can happen with multiple duct environments.'
c         write(2,*) 'Info msg: # of zeros of mode function off by one ',
c     .      nzoff,' times.  Can happen with multiple duct environments.'
       endif
c
c
       jmhi_old(nsctr)=jmhi
       jmlo_old(nsctr)=jmlo
       jmhi=0
       jmlo=0


       return
      end
ccc
      subroutine rx_bb_enter(jm1,nmode,jdf,nzsr,
     .   kn,eig_char,knbbx,vg_bb,Dvg_bb,
     .   dim1,dim2,dim3)
c
      implicit none
      include 'Parms_com'
      include 'freq2_com'
      integer*4 jm1,nmode,jdf,nzsr,
     .   jm,dim1,dim2,dim3
      real*8 vg_bb(dim2,dim3),Dvg_bb(dim2,dim3)
      complex*16 kn(0:dim2),eig_char(5,dim2),knbbx(dim2,dim3)
c
      do jm=jm1,nmode
         knbbx(jm,jdf)=kn(jm)
         vg_bb(jm,jdf)=dreal(eig_char(4,jm))
         Dvg_bb(jm,jdf)=dreal(eig_char(5,jm))
c: Updates mode counter
         mod_updt(jdf)=mod_updt(jdf) + 1
      enddo
c
      return
      end
ccc
      subroutine rx_bb_mterp(nm1,nm2,jflo,jfhi,jmlo,jmhi,
     . knbb,vg_bb,Dvg_bb,phi_bb,dphi_bb,dphiz_bb,ddphiz_bb,
     . ndv,dim1,dim2,dim3)
c
      implicit none
      include 'Parms_com'
      include 'i_o_opt_com'
      include 'i_o_2_com'
      include 'gen1_com'
      integer*4 nm1,nm2,jf1,jflo,jfhi,jmlo,jmhi,ndv,
     .   ndvm,jm1,nm1x,nm2x,nvgbad,nvgok,iibad,n2k,nm0,jdf1,
     .   dim1,dim2,dim3
      real*4 phi_bb(dim1,dim2,dim3)
      real*4 dphi_bb(dim1,dim2,dim3)
      real*4 dphiz_bb(dim1,dim2,dim3)
      real*4 ddphiz_bb(dim1,dim2,dim3)
      real*8 kmin
      real*8 vg_bb(dim2,dim3),Dvg_bb(dim2,dim3)
      complex*16 knbb(dim2,dim3)
c
      nvgbad=0
      nvgok=0
      n2k=0
c
      do jm1=1,jmlo-1
         call rx_bb_fint(jm1,jflo,jfhi,jflo,
     .      knbb,vg_bb,Dvg_bb,phi_bb,dphi_bb,dphiz_bb,
     .      ddphiz_bb,ndv,nvgok,nvgbad,n2k,dim1,dim2,dim3)
      enddo
      if(iifail .gt. 0) return
      if(jmlo .le. jmhi) then
c: Find modes from jmlo to jmhi, which can cross due to multiple ducts, 
c: using safe method:
         iibad=0
c: Min number of modes that should be found: jmhi (when nm2>jmhi) or
c: nm1 (when all modes are in ducts, nm1<=jmhi):
         nm0=min(nm1,jmhi)
c: Only compute derivatives needed for dispersion curves:
         ndvm=min(ndv,3)
         if(iidc .eq. 0) ndvm=2
         do jf1=jflo+1,jfhi-1
           jdf1=jf1-jflo+1
c: Use found mode at slightly lower freq as kmin since k increases with f:
           kmin=knbb(jmhi,jdf1-1)
           call rx_bb_mint(jmlo,jmhi,jf1,jflo,kmin,knbb,
     .       vg_bb,Dvg_bb,phi_bb,dphi_bb,dphiz_bb,ddphiz_bb,ndvm,
     .       dim1,dim2,dim3)
c added pln 21-5-97
cpln           if(jdf1-1.eq.1)knbb(jmhi,jdf1-1)=(0.d0,0.d0)
c: FIX (4-2-96):
            if(nmode .lt. nm0) then
cc          if(nmode .ne. jmhi) then
c: This should never happen:
               print *,'rx_bb_mint did not find all 1st try: ',
     .            jf1,jflo,jfhi,f_hz,jmlo,jmhi,nm1,nm2,nmode
c: Try one more time, setting kmin to a lower value:
               if(jmhi .lt. nm1) then
                  kmin=knbb(jmhi+1,1)
               else
                  kmin=0.d0
               endif
               call rx_bb_mint(nmode+1,jmhi,jf1,jflo,kmin,knbb,
     .          vg_bb,Dvg_bb,phi_bb,dphi_bb,dphiz_bb,ddphiz_bb,ndvm,
     .          dim1,dim2,dim3)
c added pln 21-5-97
cpln               if(jdf1-1.eq.1)knbb(jmhi,jdf1-1)=(0.d0,0.d0)
               if(nmode .lt. nm0) then
                  iibad=jf1
                  print *,'rx_bb_mint did not find all 1st try: ',
     .               jf1,jflo,jfhi,f_hz,jmlo,jmhi,nm1,nm2,nmode
               else
                  print *,'rx_bb_mint found 2nd try'
               endif
            endif
         enddo

         if(iibad .ne. 0) then
            jf1=iibad
            jdf1=jf1-jflo
            kmin=knbb(jmhi,jdf1)
            call rx_bb_mint(jmlo,jmhi,jf1,jflo,kmin,knbb,
     .       vg_bb,Dvg_bb,phi_bb,dphi_bb,dphiz_bb,ddphiz_bb,ndv,
     .       dim1,dim2,dim3)
c added pln 21-5-97
cpln            if(jdf1-1.eq.1)knbb(jmhi,jdf1-1)=(0.d0,0.d0)
cpln            stop
cpln added 261101 to prevent this loop
            iifail=1
            return
         endif
      endif
      do jm1=jmhi+1,nm1
cpln         write(6,*)'jm1,jflo,jfhi: ',jm1,jflo,jfhi
         call rx_bb_fint(jm1,jflo,jfhi,jflo,
     .    knbb,vg_bb,Dvg_bb,phi_bb,dphi_bb,dphiz_bb,
     .    ddphiz_bb,ndv,nvgok,nvgbad,n2k,dim1,dim2,dim3)
      enddo
      if(iifail .gt. 0) return
c
c: Do triangle from jflo to jfhi, where # modes increases from nm1 to nm2
      nm1x=nm1+1
      nm2x=min(nm1x+3,nm2)
      jf1=jflo+1
      kmin=0.d0
      do while(jf1 .lt. jfhi .and. nm1x .le. nm2x)
         call rx_bb_mint(nm1x,nm2,jf1,jflo,kmin,knbb,
     .    vg_bb,Dvg_bb,phi_bb,dphi_bb,dphiz_bb,ddphiz_bb,ndv,
     .    dim1,dim2,dim3)
c added pln 21-5-97
cpln         if(jdf1-1.eq.1)knbb(jmhi,jdf1-1)=(0.d0,0.d0)
         if(iifail.gt.0) return
cpln added 17/3/97
         if(nmode .ge. nm2x .and. jf1 .lt. jfhi-1) then
            do jm1=nm1x,nmode
               call rx_bb_fint(jm1,jf1,jfhi,jflo,
     .          knbb,vg_bb,Dvg_bb,phi_bb,dphi_bb,dphiz_bb,
     .          ddphiz_bb,ndv,nvgok,nvgbad,n2k,dim1,dim2,dim3)
            enddo
            if(iifail .gt. 0) return
            nm1x=nmode+1
            nm2x=min(nm1x+3,nm2)
         endif
         jf1=jf1 + 1
      enddo

      if (iiwrite.ge.2)
     .   print *,'Finished frequency block ',faxbb(jflo),faxbb(jfhi)
c
      return
      end
c
      subroutine interp_check(phi_bb,
     .                        jflo,jfhi,dim1,dim2,dim3)
c
      implicit none
      include 'Parms_com'
      include 'i_o_2_com'
      include 'gen1_com'
      integer*4 jflo,jfhi,jf,jm,jdf,
     .          dim1,dim2,dim3
      real*4 phi_bb(dim1,dim2,dim3),ph_err,phi_ex
      real*4 phiz(NVRMAX,NM_MAX,NFCMAX)
      real*8 kmax,kmin
c
c: TEMP: check how well interpolation worked:
      do jf=jflo+1,jfhi-1
         jdf=jf-jflo+1
         do jm=1,nmode
            phiz(2,jm,jdf)=phi_bb(2,jm,jdf)
         enddo
         nmode=0
         f_hz=faxbb(jf)
         call rx_freq_init(kmax,kmin)
c pln: correction from version 1.1 done by EKW
         kim_min=1.d100
cpln         call rx_mode_int(kmax,kmin,nm_lim,3,0,0,jdf)
         call rx_mode_int(kmax,kmin,nm_lim,3,0,1,jdf)
         do jm=1,nmode
            phi_ex=phiz(2,jm,jdf)
            ph_err=100.*abs(phi_ex - phi_bb(2,jm,jdf))/abs(phi_ex)
            if(ph_err .gt. 5. .and. abs(phi_ex) .gt. 1.e-4) then
               print *,'jf,jm,j0 = ',sngl(f_hz),jf,jm,2,ph_err,
     .            phi_ex,phi_bb(2,jm,jdf)
            endif
         enddo
      enddo
c
      return
      end
ccc
      subroutine rx_bb_mint(jm1,jm2,jf1,jflo,kminx,knbb,vg_bb,
     . Dvg_bb,phi_bb,dphi_bb,dphiz_bb,ddphiz_bb,ndv,
     . dim1,dim2,dim3)
c
      implicit none
      include 'Parms_com'
      common /eig1_com/ eig_char(5,NM_MAX)
      complex*16 eig_char
      common /out2_com/ kn(0:NM_MAX)
      complex*16 kn
      include 'i_o_2_com'
      include 'gen1_com'
      integer*4 jm1,jm2,jf1,jflo,ndv,jdf1,
     .   iim0,dim1,dim2,dim3
      real*4 phi_bb(dim1,dim2,dim3)
      real*4 dphi_bb(dim1,dim2,dim3)
      real*4 dphiz_bb(dim1,dim2,dim3)
      real*4 ddphiz_bb(dim1,dim2,dim3)
      real*8 kmin,kmax,kminx
      real*8 vg_bb(dim2,dim3),Dvg_bb(dim2,dim3)
      complex*16 knbb(dim2,dim3)
c
      f_hz=faxbb(jf1)
      call rx_freq_init(kmax,kmin)
      if(kminx .gt. kmin) kmin=kminx
      nmode=jm1-1
c
      jdf1=jf1 - jflo + 1
c pln: correction from version 1.1 done by EKW
      if(nmode .gt. 0) then
         kmax=dreal(knbb(nmode,jdf1))
         iim0=1
      else
         kim_min=1.d100
         iim0=0
      endif
c pln: correction from version 1.1 done by EKW
      call rx_mode_int(kmax,kmin,jm2,ndv,iim0,0,jdf1)
cpln added 17/3/97
      if(iifail.gt.0) return
      call rx_bb_enter(jm1,nmode,jdf1,nzsr,
     . kn,eig_char,knbb,vg_bb,Dvg_bb,
     . dim1,dim2,dim3)
c
      return
      end
ccc
      subroutine rx_bb_fint(jm1,jf1x,jf2x,jflo,
     .   knbb,vg_bb,Dvg_bb,phi_bb,dphi_bb,dphiz_bb,ddphiz_bb,
     .   ndv,nvgok,nvgbad,n2k,dim1,dim2,dim3)
c
      implicit none
      include 'Parms_com'
      include 'i_o_opt_com'
      common /eig1_com/ eig_char(5,NM_MAX)
      complex*16 eig_char
      include 'i_o_2_com'
      common /out2_com/ kn(0:NM_MAX)
      complex*16 kn
      include 'freq2_com'
      include 'gen1_com'
      integer*4 jf1x,jf2x,jflo,ndv,ii_ref,jlay_ref,
     .   jm1,jf1,jf2,jfmid,jdf1,jdf2,jdf,ncx,
     .   ndf,jf,nstack,stack(6,300),ntry,nvgbad,nvgok,
     .   nzero,jzsr,n2k,iim0,jdfmid,dim1,dim2,dim3
      real*4 phi1,phi2,Dphi1_f,Dphi2_f,
     .   f1x,f2x,twpiex,pmat(4,NVRMAX),delfx,att1,att2,delatt,attx
c fbv     .   f1x,f2x,twpiex,pmat(4,nzsr),delfx,att1,att2,delatt,attx
      real*4 phi_bb(dim1,dim2,dim3)
      real*4 dphi_bb(dim1,dim2,dim3)
      real*4 dphiz_bb(dim1,dim2,dim3)
      real*4 ddphiz_bb(dim1,dim2,dim3)
      real*4 phi_f(NVRMAX)
      real*8 vg_bb(dim2,dim3),Dvg_bb(dim2,dim3)
      real*8 f1,f2,k1,k2,k1p,k2p,
     .   k1pp,k2pp,kmin,kmax,klast,fm,fmsq,fmcu,kminx,kmid,ekmid(6),
     .   phmid,fac,dk_err,c1,c2,c3,c4,c5,c6,delf,twpie,
     .   vgmid,vgex,vg_err,vg,D2k_Df,D2k_err,kmidp,kmidpp
cc    real*8 err_w,err_k
      complex*16 knbb(dim2,dim3),knx
      data twpie/6.28318530717959/,twpiex/6.28318530718/
c
      nstack = 0
      jf1=jf1x
      jf2=jf2x
10    continue
c: Find mid frequency between f1 and f2:
      jfmid=(jf1 + jf2)/2
      f1=faxbb(jf1)
      f2=faxbb(jf2)
c: jdf's are frequency counters relative to current frequency block:
      jdf1=jf1 - jflo + 1
      jdf2=jf2 - jflo + 1
      jdf=jfmid - jflo + 1
      jdfmid = jfmid - jflo + 1
      ndf=jf2 - jf1 + 1
      fac=dfloat(ndf-1)
      fm=(jdf - jdf1)/fac
c: Initialize frequency dependent mode eigenvalue variables:
      f_hz=faxbb(jfmid)
      call rx_freq_init(kmax,kminx)
c
      k1=dreal(knbb(jm1,jdf1))
      k2=dreal(knbb(jm1,jdf2))
      k1p=twpie/vg_bb(jm1,jdf1)
      k2p=twpie/vg_bb(jm1,jdf2)
      k1pp=Dvg_bb(jm1,jdf1)
      k2pp=Dvg_bb(jm1,jdf2)
c      write(75,*)jm1,jf1,jf2
c      write(77,*)jdf1,jdf2
c      write(75,*)k1,k2
c      write(77,*)k1p,k2p
c      write(77,*)k1pp,k2pp
c      write(77,*)
      call poly5_fit(f1,f2,k1,k2,k1p,k2p,k1pp,k2pp,
     .   c1,c2,c3,c4,c5,c6,delf)
      call poly5_eval(c1,c2,c3,c4,c5,c6,delf,fm,kmid,kmidp,kmidpp,
     .   fmsq,fmcu)
      nctot=nctot + 1
      call rx_ek_calc(kmid,ekmid,phmid,ndv,iifail)
      if(iifail .gt. 0) return
      dk_err=ekmid(1)/ekmid(2)
c
      vgmid=twpie/kmidp
c: This way of getting vg may not be kosher due to exponents:
cc    vgex=-ekmid(2)/ekmid(3)
cc    err_k=g22(2,1) + migam2(1)*g12(2,1) + migam2(2)*g12(1,1)
cc    err_w=g22(3,1) + migam2(1)*g12(3,1) + migam2(3)*g12(1,1)
cc    vgex=-err_k/err_w
c
      if(abs(dk_err) .le. dk_max0) then
c pln         call rx_bb_norm_Dvg(kmid,ekmid,dk_err,vgex,D2k_Df)
         call rx_bb_Dvg(g12(1,2,1),g22(1,2,1),xkh,xkhsq,xksq(1,nlay),
     .      w,wsq,twpie,vgex,D2k_Df)

         vg_err=dabs(vgex-vgmid)/vgex
         D2k_err=dabs((D2k_Df-kmidpp)/D2k_Df)
cc       print *,'vg = ',vgex,vg,dabs(vgex-vg)/vg
         if(vg_err .lt. 1.e-4 .and. D2k_err .ge. .1d0) n2k=n2k+1
         if(vg_err .lt. 1.d-4 .and. D2k_err .lt. .1d0) then
cc    print *,'D2k_err = ',jfmid,jm1,D2k_df,kmidpp,D2k_err
c: Fit good: interpolate from jf1 to jf2:
            nvgok=nvgok + 1
            f1x=sngl(f1)
            f2x=sngl(f2)
c: Do cubic fit of mode function at each depth:
            do jzsr=1,nzsr
               phi1=phi_bb(jzsr,jm1,jdf1)
               phi2=phi_bb(jzsr,jm1,jdf2)
               Dphi1_f=twpiex*dphiz_bb(jzsr,jm1,jdf1)
               Dphi2_f=twpiex*dphiz_bb(jzsr,jm1,jdf2)
               call cub_fit_r4(f1x,f2x,phi1,phi2,Dphi1_f,Dphi2_f,
     .            pmat(1,jzsr),pmat(2,jzsr),pmat(3,jzsr),pmat(4,jzsr),
     .            delfx)
            enddo
c: Interpolate attenuation linearly:
            att1=dimag(knbb(jm1,jdf1))
            att2=dimag(knbb(jm1,jdf2))
            delatt=att2 - att1
            do jf=jf1+1,jf2-1
               fm=(jf - jf1)/fac
               call poly5_eval(c1,c2,c3,c4,c5,c6,delf,fm,kmid,
     .            kmidp,kmidpp,fmsq,fmcu)
               attx=att1 + delatt*fm
               knx=dcmplx(kmid,dble(attx))
               jdf=jf-jflo + 1
               knbb(jm1,jdf)=knx
               vg_bb(jm1,jdf)=twpie/kmidp
               do jzsr=1,nzsr
                  phi_f(jzsr)=pmat(1,jzsr) + pmat(2,jzsr)*fm + 
     .               pmat(3,jzsr)*fmsq + pmat(4,jzsr)*fmcu
                  phi_bb(jzsr,jm1,jdf)=phi_f(jzsr)
               enddo
c: Updates mode counter
               mod_updt(jdf)=mod_updt(jdf) + 1
c               write(77,*)jm1,jf,kmid
            enddo
c
            if(nstack .gt. 0) then
               jf1=stack(1,nstack)
               jf2=stack(2,nstack)
               nstack=nstack-1
               goto 10
            else
               goto 99
            endif
         else
            nvgbad=nvgbad + 1
         endif
      else
c: Compute error in group velocity for next comparison:
         vgex=-ekmid(2)/ekmid(3)
         vg_err=dabs(vgex-vgmid)/vgex
      endif
c
c: Fit bad: find mode exactly and enter into knbb,vg_bb,Dvg_bb,phi_bb:
c: Allow more than 100* error???
      if(abs(dk_err) .le. 100.d0*dk_max0 .and. vg_err .lt. 1.d-2) then
c: Guess close enough to get eigenvalue by polishing root:
         ntry=0
35       if(abs(dk_err) .gt. dk_max) then
            ntry=ntry+1
            if(ntry .gt. 2) then
               if(iidiag .gt. 0) print *,'polish ntry>2: ',ntry,jm1,
     .            jf1,jf2
            endif
c: If trouble finding mode, give up and use mode_int:
            if(ntry .gt. 4) goto 45
            kmid=kmid - dk_err
            nctot=nctot + 1
            call rx_ek_calc(kmid,ekmid,phmid,ndv,iifail)
            if(iifail .gt. 0) return
            dk_err=ekmid(1)/ekmid(2)
            goto 35
         endif
c
         nmode=jm1
         call rx_enter(nmode,ncalc,ncx,nctot,kn,kmid,
     .      ekmid,ekn,phmid,mode_phz,w,iidiag)
         call rx_norm(jlay_ref,ii_ref,ndv,vg,D2k_Df)
         if(iifail .gt. 0) return
         call rx_mf_lay(jlay_ref,ii_ref,vg,nzero)
         if(nzero .eq. nmode) then
          call rx_mode_fun_Dw(nmode,
     .         jlay_ref,ii_ref,vg,jdf)
c ffbbvv          call rx_mode_fun_Dw(nmode,phi,dphi,Dphiz_w,Ddphiz_w,
c ffbbvv     .         jlay_ref,ii_ref,vg)
c: Success using kmid guess.  Jump to rx_bb_enter:
            goto 50
c         else
c: Found wrong mode, so use mode_int below:
c           print *,'Informative Message: nzero~=nmode in rx_bb_fint ',
c    .         nzero,nmode,f_hz
         endif
      endif
c
45    continue
      nmode=jm1-1
      kmin=max(kminx,kmid - 100.d0*dk_max0)
      if(nmode .gt. 0) then
         klast=dreal(knbb(nmode,jdfmid))
         if(klast .lt. kmin) then
c: This is nothing to worry about, just a bad guess for kmid:
cc          print *,'kmin > klast: ',kmid,kmin,kminx,k1,k2
            kmin=max(k1,kminx)
         endif
         kmax=klast
         iim0=1
      else
         klast=kmax + 2.d0*dk_max0
         iim0=0
      endif
      ntry=0
c: Find next mode only:
30    continue
      ncx=nctot
      ntry=ntry + 1
cc    print *,'mode_int 30: ',ntry,kmax,kmin,knbbx(jfmid,nmode),
cc   .   kminx,kmid,k1
c pln: correction from version 1.1 done by EKW
      call rx_mode_int(kmax,kmin,jm1,ndv,iim0,0,jdf)
 9999 if(nmode .lt. jm1) then
c: No mode found over that interval, so increase interval:
         if((ntry .gt. 1).or.(iifail.gt.0)) then
            iifail=1
            write(6,*)'k1: ',k1
            write(6,*)'k2: ',k2
            print *,'trouble'
            nstack=0
            return
         else
            kmax=kmin
            kmin=max(k1,kminx)
            iim0=0
         endif
cc    print *,'trying kminx = ',kminx,kmax
         goto 30
      endif
      if(abs(klast - dreal(kn(nmode))) .lt. 2.d0*dk_max0) then
         if(iidiag .gt. 0) print *,'Info msg: dup mode??',
     .      klast,kn(nmode)
         nmode=jm1-1
         if(klast .lt. dreal(kn(nmode))) then
            kmax=klast - 2.d0*dk_max0
            iim0=0
         else
c added by PG and PLN 14/3-97--IS THIS CORRECT WE ASK
c pln            kmax=dreal(kn(nmode))-0.001*dk_max0
            kmax=dreal(kn(nmode))
            iim0=1
         endif
         kmin=max(k1,kminx)
         ntry=0
         goto 30
      endif
c
cc    print *,'success by mode_int: ',k1,kn(nmode),k2
c
50    continue
c      write(6,*)'inside rx_bb_fint'
      call rx_bb_enter(jm1,jm1,jdf,nzsr,
     . kn,eig_char,knbb,vg_bb,Dvg_bb,
     . dim1,dim2,dim3)
c
c: Push right interval onto stack if it exists:
      if(jf2 .gt. jfmid+1) then
         nstack=nstack+1
         stack(1,nstack)=jfmid
         stack(2,nstack)=jf2
      endif
      if(jfmid .eq. jf1+1) then
c: Left interval empty, check stack:
         if(nstack .gt. 0) then
            jf1=stack(1,nstack)
            jf2=stack(2,nstack)
            nstack=nstack-1
            goto 10
         endif
c: Fall through to mode loop when no intervals left on stack.
      else
c: Prepare to do left interval:
         jf2=jfmid
         goto 10
      endif
99    continue
c
      return
      end
ccc
      subroutine rx_bb_Dvg(g12,g22,xkh,xkhsq,Ksq,w,wsq,twpie,
     .   vg,D2k_Df)
c
      implicit none
      real*8 g12(6),g22(6),xkh,xkhsq,Ksq,w,wsq,twpie,vg,D2k_Df,e_k(6),
     .   migam2(6),gamsq,Dvg_w
c
      gamsq=Ksq - xkhsq
      migam2(1)=sqrt(-gamsq)
      e_k(1)=g22(1) + migam2(1)*g12(1)
      migam2(2)=xkh/migam2(1)
      migam2(3)=-Ksq/(migam2(1)*w)
      e_k(2)=g22(2) + migam2(1)*g12(2) + migam2(2)*g12(1)
      e_k(3)=g22(3) + migam2(1)*g12(3) + migam2(3)*g12(1)
      migam2(4)=-migam2(2)*migam2(3)/migam2(1)
      migam2(5)=Ksq/(migam2(1)*gamsq)
      migam2(6)=Ksq*(Ksq/gamsq - 1.d0)/(wsq*migam2(1))
      e_k(4)=g22(4) + migam2(1)*g12(4) + migam2(4)*g12(1) + 
     .   migam2(2)*g12(3) + migam2(3)*g12(2)
      e_k(5)=g22(5) + migam2(1)*g12(5) + migam2(5)*g12(1) + 
     .   2.d0*migam2(2)*g12(2)
      e_k(6)=g22(6) + migam2(1)*g12(6) + migam2(6)*g12(1) + 
     .   2.d0*migam2(3)*g12(3)
c
      vg=-e_k(2)/e_k(3)
c: Total derivative of group velocity:
      Dvg_w=(vg*vg*e_k(6) + 2.d0*vg*e_k(4) + e_k(5))/e_k(2)
      D2k_Df=(-(twpie/vg)**2)*Dvg_w
c
      return
      end
ccc
      subroutine rx_mcross(knbb,nm1,nm2,jdf2,jmlo,jmhi,
     .  ndrx,f1,f2,ccr_lo,ccr_hi,pie,dim2,dim3)
c
c: Finds mode interval over which modes might cross between frequencies
c: jf1 and jf2.
c
      implicit none
      integer*4 nm1,nm2,jmlo,jmhi,ndrx,jm,nm_min
      integer*4 jdf2,dim2,dim3
      real*8 ccr_lo,ccr_hi,pie,kcr_lo,kcr_hi
      real*4 f1,f2
      complex*16 knbb(dim2,dim3)
c
      if(ndrx .le. 1) then
c: Easy case of zero or one duct (one-to-one correspondence):
         jmlo=1
         jmhi=0
      else
         kcr_hi=2.d0*pie*f1/ccr_lo
c: EKW FIX 6-2-97: Do not let jmhi exceed # modes at high or low freqs:
         nm_min=min(nm1,nm2)
c: Find first mode at f1 that could be propagating at ss ccr_lo:
c: EKW FIX 5-9-97 (in case all modes are well trapped):
         jmlo=nm_min+1
         do jm=1,nm_min
            if(dreal(knbb(jm,1)) .le. kcr_hi) then
               jmlo=jm
               goto 10
            endif
         enddo
10       continue
         kcr_lo=2.d0*pie*f2/ccr_hi
c: Find highest order mode at f2 that could be trapped in duct with
c: highest sound speed between ducts of ccr_hi:
         do jm=nm_min,1,-1
            if(dreal(knbb(jm,jdf2)) .ge. kcr_lo) then
               jmhi=jm
               goto 20
            endif
         enddo
20       continue
c: To find modes more safely (avoid possibility of missing two at end of
c: interval), increase jmhi by one if possible:
         jmhi=min(jmhi+1,nm_min)
      endif
c
      return
      end
ccc
      subroutine cub_fit_r4(k1,k2,y1,y2,y1p,y2p,c1,c2,c3,c4,delk)
c
c: This subroutine fits a cubic polynomial y(x)=c1 + c2*x + c3*x**2 +
c: c4*x**3, where x=(k-k1)/(k2-k1), to the points (k1,y1) and (k2,y2)
c: and the derivatives y1p and y2p at those points.
c: Note that k is normalized from (k1,k2) to (0,1), so the fit values
c: must be computed accordingly: k=k1 + kx*delk.
c
      implicit none
      real*4 k1,k2,y1,y2,y1p,y2p,c1,c2,c3,c4,b1,b2,delk,y1px,y2px
c
      delk=k2 - k1
      y1px=y1p*delk
      y2px=y2p*delk
      c1=y1
      c2=y1px
      b1=y2 - y1 - y1px
      b2=y2px - y1px
      c3=3.d0*b1 - b2
      c4=b2 - 2.d0*b1
c
      return
      end
cc
