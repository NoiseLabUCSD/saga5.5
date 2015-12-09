      subroutine bb_fmarch(jmk,jmo,k,r1r2,kp0,Lp0,dL_dwp0,dL_dkp0,
     .   iisg,dfinc,dfincx,nm2,jfbb2,iixi_cut,kdone)
c
c: Marches mode to new frequency.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jmk,jmo,iisg,ntry,nm2,jfbb2,iixi_cut,kdone
      complex*16 k,r1r2(3,4),k0,rr0(3,4),kp0,Lp0,dL_dwp0,dL_dkp0,
     .   kp,Lp,dL_dwp,dL_dkp,lnr1r2,branch,delk,dW_dk
      real*8 f_hz0,delf,dfinc,dfincx,df_fac,Lmagsq,Lerrok,magsq,
     .   frac_done,f_final
      real*8 xi_real,xi_rat
c
      frac_done=0.
      f_hz0=f_hz
      f_final=f_hz + iisg*df
c: Save input values of kp, Lp, etc. in case merged modes found:
      kp=kp0
      Lp=Lp0
      dL_dwp=dL_dwp0
      dL_dkp=dL_dkp0
c
10    ntry=0
20    continue
      ntry=ntry + 1
      delf=iisg*dfinc*df
      f_hz=f_hz0 + delf
      if(f_hz .lt. fmin) then
         f_hz=fmin
         delf=f_hz - f_hz0
         dfinc=delf/(iisg*df)
      elseif(f_hz .gt. fmax) then
         f_hz=fmax
         delf=f_hz - f_hz0
         dfinc=delf/(iisg*df)
      endif
c: Obtain guess for eigenvalue at new frequency using derivatives with
c: respect to k and w:
cxx   print *,'delk dir = ',atan2(dimag(delk),dreal(delk))*180./pie
      call freq_chng
      delk=-(Lp + dL_dwp*delf*twpie)/dL_dkp
      k=kp + delk
      call r1r2_calc(k,r1r2,3,1,jjfail)
      if(jjfail .gt. 0) return
      nctot=nctot + 1
      lnr1r2=r1r2(1,4)
      if(ntry .gt. 5) then
         print *,'jm,f,ntry = ',jmo,f_hz,ntry,dfinc
         print *,'Re(k): ',dreal(kp0)/kw0,dreal(kp)/kw0,dreal(k)/kw0
         print *,'Im(k): ',dimag(kp0)*8686.,dimag(kp)*8686.,
     .      dimag(k)*8686.
         print *,'Lp0,Lp,L = ',Lp0,Lp,lnr1r2,Lmagsq,Lerrok
         print *,'delk = ',dreal(delk)/kw0,dimag(delk)*8686
         print *,'iish,iish_ref,iishn = ',iish,iish_ref,iishn(jmk)
c        call r1r2_calc(kp,r1r2,3,1)
c        print *,'kp again: ',Lp,r1r2(1,4),dL_dkp,r1r2(2,4),
c    .      dL_dwp,r1r2(3,4)
cc       iidiag=1
         if(ntry .gt. 8) then
            print *,'Stopping ...'
            stop
         endif
      endif
c
      Lmagsq=magsq(lnr1r2)
      Lerrok=magsq(r1r2(2,4))*errdkms
      dfincx=dfinc
c: FIX 4-12-93: Check magnitude squared of ln(R1*R2) also:
      if(Lmagsq .le. Lerrok .and. Lmagsq .lt. .01d0) then
c: Guess perfect:
         if(abs(dimag(lnr1r2)) .gt. .157) then
            print *,'arg(r1r2) > .157 !!',nmode,f_hz,lnr1r2
         endif
         if(abs(dreal(lnr1r2)) .gt. .1) then
            print *,'mag(r1r2) > .1 !!',nmode,f_hz,lnr1r2
         endif
         if(Lerrok/Lmagsq .gt. 6.d0) dfincx=2.*dfinc
cc       if(Lerrok/Lmagsq .gt. 6.d0) dfinc=2.*dfinc
      elseif(abs(dreal(lnr1r2)) .lt. .1 .and.
     .   abs(dimag(lnr1r2)) .lt. .3927) then
c: Guess ok, but go back if trying to jump more than df since linear 
c: interpolation won't be so good here:
         if(dfinc .gt. 1.01d0) then
            df_fac=(Lerrok/Lmagsq)**(0.25)
ctt   print *,'Guess OK, but dfinc > 1: ',dfinc,df_fac,dfinc*df_fac
            dfinc=dmax1(1.d0,dmin1(dfinc-1.d0,
     .         dfloat(nint(dfinc*df_fac))))
ctt   print *,'dfinc set to ',dfinc
            call xkh_backup
            goto 20
         endif
c: Use eig_final to zero in on eigenvalue:
         call eig_final(k,r1r2,errdkms,0,kw0,iish,iifail,jjfail)
         if(iifail .eq. 1) then
            print *,'iifail=1 in eig_final from bb_fmarch'
            return
         endif
      else
c: Guess bad:
         dfinc=.5d0*dfinc
         if(dfinc .gt. 1.d0) dfinc=dfloat(nint(dfinc))
         call xkh_backup
         goto 20
      endif
c
      if(dimag(k) .lt. 0.d0 .and. dimag(k)*dimag(k) .gt. errdkms) 
     .   then
cxx      if(iidiag .ne. 0) print *,'MODE FOUND IN LOWER HP: ',
cxx  .      nmode,k,k/kw
cxx      iilhp=1
cxx      goto 888
         if(iidiag .ge. 1) print *,'neg im(k): ',k,k/kw0,f_hz,iish
      endif
c: Check for negative sheets, but don't change sheet to pos if found:
      if(isp(nlay) .eq. 1) then
         if(iish(1,1) + iish(2,1) .lt. 2) then
            kdone=1
            print *,'Mode crossed Pekeris cut: ',jmk,
     .         dreal(k)/kw0,dimag(k)*8685.9
            return
         endif
         call sheet_look(0,branch)
         if(branch .ne. (0.,0.)) then
            if(dfinc .gt. 1.001d0) then
               dfinc=1.d0
               call xkh_backup
               goto 20
            endif
            if(dreal(k) .lt. dreal(branch) .or. dimag(k) .lt. 
     .         dimag(branch)) then
               call fix_path(k0,rr0,1)
               if(iidone .eq. 1 .or. iifail .eq. 1) then
                  print *,'bb fix_path trouble',iidone,iifail
                  return
               endif
c: EKW FIX: make phcut=2*pi so that dph_left is close to zero in dph_calc:
               phcut=twpie
               call eig_find0(k0,rr0,k,r1r2,nmode+10)
               if(iidone .eq. 1 .or. iifail .eq. 1) then
                  print *,'bb eig_find0 trouble',iidone,iifail
                  return
               endif
            else
               if(iidiag .ge. 1) print *,'OK(?) MODE ON -1 SHEET ',
     .            'FOUND: ',nmode,sngl(f_hz),iish
            endif
         endif
cc    elseif(dfinc .lt. .999d0) then
cc    elseif(iixi_cut .eq. 0) then
      else
         if(-iisg*dimag(xi_hsp(1,1)) .lt. 0.d0) then
            kdone=1
            print *,'fmarch crossed gradient BL: ',jmk,xi_hsp(1,1),
     .         dreal(k)/kw0,dimag(k)*8685.9
            return
         elseif(iisol(nlay) .eq. 1) then
            if(-iisg*dimag(xi_hsp(1,2)) .lt. 0.d0) then
               kdone=1
               print *,'fmarch crossed gradient BL: ',jmk,xi_hsp(1,1),
     .            dreal(k)/kw0,dimag(k)*8685.9
               return
            endif
         endif
         xi_real=dreal(xi_hsp(1,1))
         xi_rat=dabs(xi_real/dimag(xi_hsp(1,1)))
         if(nsvmin .ne. nlay) then
            xi_real=dreal(xi_hsp(1,1))
            xi_rat=dabs(xi_real/dimag(xi_hsp(1,1)))
            if(xi_rat .gt. 20. .and. xi_real .lt. -1.) then
ccc         dW_dk=-2.*r1r2(1,1)*gamiref*r1r2(2,3)
ccc         call mode_fun(k,r1r2(1,1),r1r2(1,2),dW_dk,phi,dphi,
ccc  .         psi,dpsi,exp_gbs,nm2)
c: Compute mode amplitudes at duct reference depths:
ccc         call phimag_check(phi,dphi,phim_bb(jfbb2),k,nm2)
ccc         if(iidone .eq. 1) then
      print *,'BLM CROSS: ',jmo,sngl(f_hz),sngl(dfinc),xi_rat,xi_real
ccc            iidone=0
               f_hz=f_final
               call freq_chng
               call bb_blm_exit2(k,r1r2)
               iixi_cut=1
               dfinc=1.
               dfincx=1.
            endif
         endif
      endif
      nsave=nsave + max(0,nint(dfinc-1.d0))
cpln      write(6,*)'nsave,max0: ',nsave,max(0,nint(dfinc-1.d0)),dfinc
cpln      pause
c     
c: Check if entire freq bin has been crossed:
      if(dfinc .lt. .999d0) then
c: If only jumped by a fraction of df, check if a whole df done:
         frac_done=frac_done + dfinc
c: Allow dfinc to increase again:
         dfinc=dfincx
         if(frac_done .lt. .999) then
c: If total fraction less than 1, update f_hz0 and jump more in f:
            if(iift .eq. 1) call bb_write(jmo,jfbb2,faxbb,k,kw0,iish)
            f_hz0=f_hz
            kp=k
            Lp=r1r2(1,4)
            dL_dkp=r1r2(2,4)
            dL_dwp=r1r2(3,4)
            dfinc=min(dfinc,1.-frac_done)
            goto 10
         endif
      endif
c
      return
      end
ccc
      subroutine bb_blm_exit(k,r1r2)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      real*8 I0_re,del_xi_zero
      complex*16 k,r1r2(3,4),k0,rr0(3,4),beta,delk,lnr1r2,I_k,del_I,
     .   xi_lower
c
      k0=k
      beta=(etasq(nlay)/k0)*dcmplx(0.d0,-1.d0)
      delk=(.01d0*pie/Htot)*beta/cdabs(beta)
      k0=k0 + delk
      xi_lower=xi_hsp(1,1) - 
     .   dcmplx(0.d0,min(2.d0,-.05d0*dreal(xi_hsp(1,1))))
      k0=cdsqrt(xksq(1,nlay) + xi_lower*etasq(nlay))

      del_xi_zero=pie/dsqrt(dabs(dreal(xi_hsp(1,1))))
      xi_lower=xi_hsp(1,1) - dcmplx(0.d0,4.d0*del_xi_zero)
      k0=cdsqrt(xksq(1,nlay) + xi_lower*etasq(nlay))
      call r1r2_calc(k0,rr0,3,1,jjfail)
      if(jjfail .gt. 0) return
c
10    continue
      lnr1r2=rr0(1,4)
      I0_re=dreal(lnr1r2)
      if(dabs(I0_re) .gt. .05d0) then
         I_k=rr0(2,4)
         beta=(-I_k*etasq(nlay)/k0)
         del_I=dcmplx(-I0_re,-I0_re*dimag(beta)/dreal(beta))
         delk=del_I/I_k
         k0=k0 + delk
c
         call r1r2_calc(k0,rr0,3,1,jjfail)
         if(jjfail .gt. 0) return
         print *,'I = ',dreal(k0)/kw0,dimag(k0)*8685.9,lnr1r2,
     .      rr0(1,4),xi_hsp(1,1)
         goto 10
      endif
c
c: Cal traj_mag if not close enough to mode:
      iimst=-1
      phcut=2.d0*pie
      call eig_find0(k0,rr0,k,r1r2,0)
cc    dW_dk=-2.*r1r2(1,1)*gamiref*r1r2(2,3)
cc    call mode_fun(k,r1r2(1,1),r1r2(1,2),dW_dk,phi,dphi,
cc   .   psi,dpsi,exp_gbs,nm2)
c
      return
      end
ccc
      subroutine bb_blm_exit2(k,r1r2)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      real*8 del_xi_zero,ln_re1,ln_re2,ln_re0
      complex*16 k,r1r2(3,4),k0,rr0(3,4),xi1,xi2,xi0
      integer*4 ntry
c
      xi1=xi_hsp(1,1)
      call i_xi_calc(xi1,k0,rr0,ln_re1)
c
      ntry=0
      if(ln_re1 .lt. 0.d0) then
10       continue
         xi2=dcmplx(2.d0*dreal(xi1),dimag(xi1))
         ntry=ntry + 1
         if(ntry .gt. 10) print *,'ntry trouble in exit'
         call i_xi_calc(xi2,k0,rr0,ln_re2)
         if(ln_re2 .lt. 0.d0) then
            xi1=xi2
            goto 10
         endif
      else
         xi2=xi1
15       continue
         xi1=dcmplx(0.50*dreal(xi2),dimag(xi2))
         ntry=ntry + 1
         if(ntry .gt. 10) print *,'ntry trouble in exit'
         call i_xi_calc(xi1,k0,rr0,ln_re1)
         if(ln_re1 .gt. 0.d0) then
            xi2=xi1
            goto 15
         endif
      endif
c
c: We now have |R1R2|=1 bracketed on xi axis:
20    xi0=0.5d0*(xi1 + xi2)
      call i_xi_calc(xi0,k0,rr0,ln_re0)
      ntry=ntry + 1
      if(ntry .gt. 20) print *,'ntry trouble in exit'
      if(dabs(ln_re0) .gt. 0.2d0) then
         if(ln_re0 .lt. 0.d0) then
            xi1=xi0
         else
            xi2=xi0
         endif
         goto 20
      endif
      print *,'ntry final = ',ntry
c
c: Call traj_mag if not close enough to mode:
      iimst=-1
      phcut=2.d0*pie
cc    iimst=0
      call eig_find0(k0,rr0,k,r1r2,0)
c
      return
      end
ccc
      subroutine i_xi_calc(xi0,k0,rr0,i0)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      complex*16 xi0,k0,rr0(3,4)
      real*8 i0,del_xi_zero
c
      del_xi_zero=pie/dsqrt(dabs(dreal(xi0)))
      xi0=dcmplx(dreal(xi0),dimag(xi0)-4.d0*del_xi_zero)
      k0=cdsqrt(xksq(1,nlay) + xi0*etasq(nlay))
      call r1r2_calc(k0,rr0,3,1,jjfail)
      if(jjfail .gt. 0) return
      i0=dreal(rr0(1,4))
c
      return
      end
