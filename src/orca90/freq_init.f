      subroutine freq_init
c
c: Initializes variables associated with a new frequency.
c
      use parms_com
      use i_o_com
      use gen_com
c
      integer*4 j,ii,jsr,inc,ii1,ii2,jlay
      real*8 zfac
c
      w=twpie*f_hz
      wsq=w*w
c: Compute minimum Re(k) to search for eigenvalues:
      if(cphmax .gt. 0.) then
         kremin=w/cphmax
      else
         kremin=-1.d100
      endif
      kcrmin=w/crmax
      khspmin=w/chspmax
c
      if(nthorpe .gt. 0) call thorpe_attn
c
c: Set Airy halfspace thicknesses to 3 wavelengths:
      if(ihf(1) .eq. 1) h(1)=3.d0*geo(2,1,1)/f_hz
      if(ihf(nlay) .eq. 1) h(nlay)=3.d0*geo(1,1,nlay)/f_hz
c
      do j=1,nlay
         call k_attn(w,geo(1,1,j),geo(1,4,j),fexp(1,j),xk(1,j),
     .      xksq(1,j))
         call k_attn(w,geo(2,1,j),geo(2,4,j),fexp(1,j),xk(2,j),
     .      xksq(2,j))
         call eta_calc(xksq(1,j),xksq(2,j),h(j),eta(j),etasq(j),
     .      isp(j),j,nlay)
         if(iisol(j) .eq. 1) then
            call k_attn2(w,geo(1,2,j),geo(1,5,j),fexp(2,j),xb(1,j),
     .         xbsq(1,j),xbsqinv(1,j))
            call k_attn2(w,geo(2,2,j),geo(2,5,j),fexp(2,j),xb(2,j),
     .         xbsq(2,j),xbsqinv(2,j))
            call eta_calc(xbsq(1,j),xbsq(2,j),h(j),etb(j),etbsq(j),
     .         iss(j),j,nlay)
         endif
      enddo
c
      kw0=w/cfmin
      xkref=xk(isvmin,nsvmin)
c: Set branch point ratio for reference depth:
      xkrat_ref(1)=xkref/kw0
      kw=dreal(xkref)
      cref=w/kw0
      errdkms=dmin1(errdkms,(1.d-5*kw)**2)
      errdk100=-100.d0*dsqrt(errdkms)
c
      do ii=1,2 
c: xkbp(iihs,iips) holds the branch points for halfspace iihs (1=lower,
c: 2=upper) and wave type iips (1=p-wave,2=s-wave):
         ii1=ii
         ii2=3 - ii
         xkbp(ii,1)=xk(ii1,jsol(ii,2)+jflu(ii,3))
         xkbp(ii,2)=xb(ii1,jsol(ii,2)+jflu(ii,3))
         xkrat(ii,1)=xkbp(ii,1)/kw0
         xkrat(ii,2)=xkbp(ii,2)/kw0
         if(iidiag .ge. 2) then
            print *,'BP = ',ii,xkbp(ii,1),xkbp(ii,2)
            print *,'k/kw = ',ii,xkrat(ii,1),xkrat(ii,2)
         endif
         inc=jflu(ii,3)
         if(allf(ii) .eq. 0) then
            do jsr=1,nzsr
               if(kksh(jsr) .eq. 1) then
                  jlay=jsr2j(jsr)
                  zfac=(zsr(jsr)-zdep(jlay-1))/h(jlay)
                  ksm2_sr(jsr)=1.d0/(xbsq(1,jlay) + 
     .               zfac*(xbsq(2,jlay)-xbsq(1,jlay)))
               else
                  ksm2_sr(jsr)=(0.d0,0.d0)
               endif
            enddo
         endif
      enddo
c
      return
      end
ccc
      subroutine k_attn(w,c0,alpha,fexp,xk,xksq)
c
c: Computes complex wavenumber for medium with sound speed c0 and
c: attenuation alpha (in dB/(m-kHz)).
      implicit none
      complex*16 xk,xksq
      real*8 w,c0,alpha,fexp,db2nep,w0,alpha_w
c: db2nep = 20*log10(e), conversion from dB to nepers:
      data db2nep/8.68588963806504/
c: w0 = 2*pi*f0 = 2*pi*1000:
      data w0/6.283185307179586d+03/
c: fac=2*pi*1000*20*log(e):
cc    data fac/5.457505415367364d+04/
       
c
      if(c0 .ne. 0.) then
c: Convert alpha (db/m-kHz) to alphap (db/m):
cc       xk=(w/c0)*dcmplx(1.d0,alpha*c0/fac)
c: New 3-26-99. Incorporate power law for attenuation vs freq such that
c: when fexp=1 or w=w0, we get alpha_w=alpha:
         alpha_w=alpha * (w/w0)**fexp
         xk=dcmplx(w/c0,alpha_w/db2nep)
         xksq=xk*xk
      endif
c: alpha(dB/m/kHz) = (f/1000)^(fexp-1) * alpha(dB/m @ 1kHz)
c
      return
      end
ccc
      subroutine k_attn2(w,c0,alpha,fexp,xk,xksq,xksqinv)
c
c: Computes complex wavenumber for medium with sound speed c0 and
c: attenuation alpha (in dB/(m-kHz)).
      implicit none
      complex*16 xk,xksq,xksqinv
      real*8 w,c0,alpha,fexp,db2nep,w0,alpha_w
c: db2nep = 20*log10(e), conversion from dB to nepers:
      data db2nep/8.68588963806504/
c: w0 = 2*pi*f0 = 2*pi*1000:
      data w0/6.283185307179586d+03/
c: fac=2*pi*1000*20*log(e):
cc    data fac/5.457505415367364d+04/
c
      if(c0 .ne. 0.) then
c: Convert alpha (db/m-kHz) to alphap (db/m):
cc       xk=(w/c0)*dcmplx(1.d0,alpha*c0/fac)
c: New 3-26-99. Incorporate power law for attenuation vs freq such that
c: when fexp=1 or w=w0, we get alpha_w=alpha:
         alpha_w=alpha * (w/w0)**fexp
         xk=dcmplx(w/c0,alpha_w/db2nep)
         xksq=xk*xk
         xksqinv=1.d0/xksq
      else
         xksqinv=(0.d0,0.d0)
      endif
c: alpha(dB/m/kHz) = (f/1000)^(fexp-1) * alpha(dB/m @ 1kHz)
c
      return
      end
ccc
      subroutine eta_calc(xk1sq,xk2sq,h,eta,etasq,iso,jlay,nlay)
c
      implicit none
      integer*4 iso,jlay,nlay
      complex*16 xk1sq,xk2sq,eta,etasq,dk,eip23
      real*8 h,third,sqrt3,pie23
      data third/0.333333333333333d0/,sqrt3/1.73205080756888/
     .   eip23/(-0.5d0,0.86602540378444d0)/,
     .   pie23/2.09439510239320d0/
c
      if(iso .eq. 1 .or. h .le. 0.d0) then
         eta=(0.d0,0.d0)
         etasq=(0.d0,0.d0)
         return
      endif
      if(jlay .ne. 1 .and. jlay .ne. nlay) then
         dk=xk2sq - xk1sq
         if(dreal(dk) .gt. 0.d0) then
            eta=(dk/h)**third
         else
c: For dk/h in left half plane, make sure eta lies close to neg
c: real axis:
            eta=-((-dk/h)**third)
         endif
      else
c: For Airy halfspaces, make sure line in -eta direction in xi-plane 
c: points between +-pi/3 direction:
         dk=xk2sq - xk1sq
         eta=-((-dk/h)**third)
      endif
      etasq=eta*eta
c
      return
      end
ccc
      subroutine freq_chng
c
c: Changes variables associated with a new frequency after freq_init has 
c: already been called.
c
      use parms_com
      use i_o_com
      use gen_com
c
      integer*4 j,ii,jx1,jx,jsr,inc
      real*8 w_old,w_rat,w_ratsq,w_23,w_43,rat23
c
c: Since attenuation with non-unity frequency exponents (fexp) possible
c: always call freq_init:
      call freq_init
      return
c
c: If Thorpe layers exist, might as well call freq_init:
      if(nthorpe .gt. 0) then
         call freq_init
         return
      endif
c
      w_old=w
      w=twpie*f_hz
      w_rat=w/w_old
      w_ratsq=w_rat*w_rat
      rat23=2.d0/3.d0
      w_23=w_rat**rat23
      w_43=w_23*w_23
c: Compute minimum Re(k) to search for eigenvalues:
      if(cphmax .gt. 0.) then
         kremin=w/cphmax
      else
         kremin=-1.d100
      endif
      kcrmin=w/crmax
c
      do j=1,nlay
         call xk_chng(xk(1,j),xksq(1,j),w_rat,w_ratsq,eta(j),
     .      etasq(j),w_23,w_43,isp(j),h(j),ihf(j))
         if(iisol(j) .eq. 1) then
            call xk_chng2(xb(1,j),xbsq(1,j),w_rat,w_ratsq,etb(j),
     .         etbsq(j),w_23,w_43,iss(j),xbsqinv(1,j),ihf(j))
         endif
      enddo
c
      xkref=xk(isvmin,nsvmin)
      kw=dreal(xkref)
      kw0=w/cfmin
      errdkms=dmin1(errdkms,(1.d-5*kw)**2)
      errdk100=-100.d0*dsqrt(errdkms)
c
      do ii=1,2 
c: xkbp(iihs,iips) holds the branch points for halfspace iihs (1=lower,
c: 2=upper) and wave type iips (1=p-wave,2=s-wave):
         xkbp(ii,1)=xk(ii,jsol(ii,2)+jflu(ii,3))
         xkbp(ii,2)=xb(ii,jsol(ii,2)+jflu(ii,3))
         xkrat(ii,1)=xkbp(ii,1)/kw0
         xkrat(ii,2)=xkbp(ii,2)/kw0
         inc=jflu(ii,3)
         if(allf(ii) .eq. 0) then
            do j=jsol(ii,1),jsol(ii,2)+inc,inc
               jx1=jzmx(j)
               do jx=jx1,jx1+nzmx(j)-1
                  jsr=jsrmx(jx)
                  ksm2_sr(jsr)=ksm2_sr(jsr)/w_ratsq
               enddo
            enddo
         endif
      enddo
c
      return
      end
ccc
      subroutine xk_chng(xk,xksq,w_rat,w_ratsq,eta,etasq,w_23,
     .   w_43,isp,h,ihf)
c
      implicit none
      integer*4 isp,ihf
      complex*16 xk(2),xksq(2),eta,etasq
      real*8 w_rat,w_ratsq,w_23,w_43,h
c
      xk(1)=xk(1)*w_rat
      xksq(1)=xksq(1)*w_ratsq
      if(isp .eq. 0) then
         xk(2)=xk(2)*w_rat
         xksq(2)=xksq(2)*w_ratsq
         if(ihf .eq. 0) then
            eta=eta*w_23
            etasq=etasq*w_43
         else
            eta=eta*w_rat
            etasq=etasq*w_ratsq
c: Layer thickness h inversely proportional to frequency:
            h=h/w_rat
         endif
      else
         xk(2)=xk(1)
         xksq(2)=xksq(1)
c: Layer thickness h inversely proportional to frequency:
         if(ihf .ne. 0) h=h/w_rat
      endif
c
      return
      end
ccc
      subroutine xk_chng2(xk,xksq,w_rat,w_ratsq,eta,etasq,w_23,
     .   w_43,isp,xbsqinv,ihf)
c
      implicit none
      integer*4 isp,ihf
      complex*16 xk(2),xksq(2),eta,etasq,xbsqinv(2)
      real*8 w_rat,w_ratsq,w_23,w_43
c
      xk(1)=xk(1)*w_rat
      xksq(1)=xksq(1)*w_ratsq
      xbsqinv(1)=xbsqinv(1)/w_ratsq
      if(isp .eq. 0) then
         xk(2)=xk(2)*w_rat
         xksq(2)=xksq(2)*w_ratsq
         xbsqinv(2)=xbsqinv(2)/w_ratsq
         if(ihf .eq. 0) then
            eta=eta*w_23
            etasq=etasq*w_43
         else
            eta=eta*w_rat
            etasq=etasq*w_ratsq
         endif
      else
         xk(2)=xk(1)
         xksq(2)=xksq(1)
         xbsqinv(2)=xbsqinv(1)
      endif
c
      return
      end
ccc
      subroutine thorpe_attn
c
c: Set attenuation in Thorpe attenuation layers:
c
      use parms_com
      use i_o_com
      use gen_com
c
      real*8 fsq,attn_thorpe,a_nu_db_kyd
      integer*4 j,jlay
c
      fsq=(1.d-3*f_hz)**2
      do j=1,nthorpe
         jlay=jthorpe(j)
	 a_nu_db_kyd=.1*fsq/(1.+fsq)+40.*fsq/(4100.+fsq)+2.75e-4*fsq
	 attn_thorpe=a_nu_db_kyd*(1./914.4)/(.001*f_hz)
         geo(1,4,jlay)=attn_thorpe
         geo(2,4,jlay)=attn_thorpe
      enddo
c
      return
      end
