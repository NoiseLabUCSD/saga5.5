       subroutine bb_modes
c
c: Computes the broadband field from fmin to fmax in steps of df
c: using mode theory.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jm,nm_fmax,nm_fmin,nmbbtot,iibad,nmiss,ndied,nrose
      real*4 ncperm
c
      nctot=0
      nclast=0
      nmerge=0
      ncperm=0.e0
      nctot=0
      nclast=0
      nmode=0
      nm_fmax=0
      nsave=0
c
c: Check input variables, sizes of arrays, and initialize variables:
      call bb_init
c
c
      nmode=0
      call zmx_init
      f_hz=fmax
      call mode_find(1)
      if(jjfail.gt.0) return
      kim_bb(nfbb)=kim_min
      phim_bb(nfbb)=phi_mag_max
      nm_fmax=nmode
      print *,'fmax,nmode = ',fmax,nmode
c
c: Check if knbb and eig_bb are large enough:
      iibad=0
      call mem_lim(nmode*nfbb,NM_NF_MAX,MLINE,LML,'nmode*nfbb',10,
     .   'NM_NF_MAX',9,iibad,1)
      if(iibad .eq. 1) stop
c
c: Open mode characteristics file if desired:
      if(iiout .ne. 0) call bb_out_init(nm_fmax)
c
      call bb_ftinit
c
      do jm=1,nm_fmax
         kn_indx(jm)=jm
      enddo
      nmbbtot=nm_fmax
      call bb_mloop(0,nm_fmax,nmbbtot,-1,fmax,nfbb,ndied)
      if(isp(nlay) .eq. 0) then
         call bb_blmodes
      endif
c
c: Find modes at lowest frequency:
      f_hz=fmin
      call mode_find(0)
      if(jjfail.gt.0)return
      nm_fmin=nmode
c
      call bb_miss(knbb,nm_fmax,kn(1),nm_fmin,kn(nm_fmin+3),nmiss)
      if(nmiss .gt. 0) then
         call bb_mloop(nm_fmax,nmiss,nmbbtot,1,fmin,1,nrose)
         nm_fmax=nm_fmax + nmiss
      else
         nrose=0
      endif
      if(ndied .gt. nrose) then
         print *,'ndied>nrose: ',ndied,nrose
      endif
c
      call bb_done(nm_fmax)
c
      ncperm=float(nctot)/float(max(1,nmbbtot))
      print *,'nmbbtot,nctot,nsave= ',nmbbtot,nctot,nsave,ncperm
      if(nmerge .gt. 0) then
         print *,'# modes found to merge: ',nmerge
         print *,'# modes skipped due to merging = ',nskip
      endif
      write(lusvp,120) nmbbtot,nctot,ncperm,nsave
120   format('TOTAL # MODES = ',i8,'; # R1R2 CALCS = ',i8,
     .   '; #CALCS/MODE = ',f5.2,'; NSAVE = ',i6)
c
      return
      end
ccc
      subroutine bb_interp(phi,dphi,psi,dpsi,exp_gbs,phix,dphix,
     .   psix,dpsix,expx_gbs,nzsr,nm1,nm2,nm3,iisg,nfinc,knbb,tf,
     .   nfbb,jfbb,nmode,jmo,faxbb,cfmin,iish,eig_bb,iift,kim_bb,
     .   phim_bb,iish_bb)
c
      implicit none
      integer*4 nzsr,nm1,nm2,nm3,iisg,nfinc,nfbb,jfbb,nmode,jmo,
     .   iish(2,2)
      integer*4 jf1,jf2,jf,jz,jchar,iift,iish_bb(nfbb,nmode),jfx
      real*4 faxbb(nfbb),phim_bb(nfbb)
      real*8 kim_bb(nfbb),cfmin,kw0x,twpie
      complex*8 phi(nzsr,nm3),dphi(nzsr,nm3),psi(nzsr,nm3),
     .   dpsi(nzsr,nm3),exp_gbs(nzsr,nm3),phix(nzsr),dphix(nzsr),
     .   psix(nzsr),dpsix(nzsr),expx_gbs(nzsr),tf
      complex*16 knbb(nfbb,nmode),eig_bb(5,nfbb,nmode),fac
      data twpie/6.28318530717959/
c
c: Interpolate eigenvalues and mode characteristics from jf1 to jf2:
      jf1=jfbb
      jf2=jfbb + iisg*nfinc
      fac=(knbb(jf2,jmo) - knbb(jf1,jmo))/nfinc
      do jf=1,nfinc-1
         jfx=jf1+iisg*jf
         knbb(jfx,jmo)=knbb(jf1,jmo) + jf*fac
         iish_bb(jfx,jmo)=iish_bb(jf1,jmo)
      enddo
      do jchar=1,5
         fac=(eig_bb(jchar,jf2,jmo) - eig_bb(jchar,jf1,jmo))/nfinc
         do jf=1,nfinc-1
            jfx=jf1+iisg*jf
            eig_bb(jchar,jfx,jmo)=eig_bb(jchar,jf1,jmo) + jf*fac
         enddo
      enddo
c
c: Use nm3 column as temporary storage for interpolation factors:
      do jz=1,nzsr
         phi(jz,nm3)=(phi(jz,nm2)-phi(jz,nm1))/nfinc
         dphi(jz,nm3)=(dphi(jz,nm2)-dphi(jz,nm1))/nfinc
         psi(jz,nm3)=(psi(jz,nm2)-psi(jz,nm1))/nfinc
         dpsi(jz,nm3)=(dpsi(jz,nm2)-dpsi(jz,nm1))/nfinc
         exp_gbs(jz,nm3)=(exp_gbs(jz,nm2)-exp_gbs(jz,nm1))/nfinc
      enddo
      do jf=1,nfinc-1
         do jz=1,nzsr
            phix(jz)=phi(jz,nm1)+jf*phi(jz,nm3)
            dphix(jz)=dphi(jz,nm1)+jf*dphi(jz,nm3)
            psix(jz)=psi(jz,nm1)+jf*psi(jz,nm3)
            dpsix(jz)=dpsi(jz,nm1)+jf*dpsi(jz,nm3)
            expx_gbs(jz)=exp_gbs(jz,nm1)+jf*exp_gbs(jz,nm3)
         enddo
c: Compute field at all receivers at interpolated freq between f_hz0 & f_hz:
         jfx=jf1+iisg*jf
         call bb_field(knbb(jfx,jmo),phix,dphix,dpsix,expx_gbs,1,tf,
     .      jfx,jmo)
         if(iift .eq. 1) then
            kw0x=twpie*faxbb(jfx)/cfmin
            call bb_write(jmo,jfx,faxbb,knbb(jfx,jmo),kw0x,iish)
         endif
         kim_bb(jfx)=dmin1(kim_bb(jf1),kim_bb(jf2))
         phim_bb(jfx)=min(phim_bb(jf1),phim_bb(jf2))
      enddo
c
      return
      end
ccc
      subroutine bb_enter(k,lnrr,vg,R1,knbb,eig_bb,jfbb,jmo,nfbb,
     .   nmode,iish,iish_ref,iish_bb,nmbb)
c
c: Fills knbb and eig_bb with eigenvalue and mode characteristics.
c
      implicit none
      integer jfbb,nfbb,jmo,nmode,iish(2,2),iish_ref(2),
     .   iish_bb(nfbb,nmode),nmbb(nfbb)
      complex*16 k,lnrr(3),vg,R1,knbb(nfbb,nmode),eig_bb(5,nfbb,nmode)
c
      knbb(jfbb,jmo)=k
      eig_bb(1,jfbb,jmo)=lnrr(1)
      eig_bb(2,jfbb,jmo)=lnrr(2)
      eig_bb(3,jfbb,jmo)=lnrr(3)
      eig_bb(4,jfbb,jmo)=vg
      eig_bb(5,jfbb,jmo)=R1
      nmbb(jfbb)=jmo
      call iish_code(iish,iish_ref,iish_bb(jfbb,jmo),1)
c
      return
      end
ccc
      subroutine bb_write(jmo,jfbb,faxbb,k,kw0,iish)
c
      implicit none
      integer*4 jmo,jfbb,iish(2,2)
      real*4 faxbb(jfbb)
      real*8 kw0,twpie
      complex*16 k
      data twpie/6.28318530717959/
c
ccc   write(14,100) jmo,jfbb,faxbb(jfbb),dreal(k)/kw0,dimag(k)/kw0,
ccc  .   iish(1,1),iish(2,1),iish(1,2),iish(2,2)
c: Output IM(k) in terms of dB/km:
      write(14,100) jmo,jfbb,faxbb(jfbb),dreal(k)/kw0,8685.9*dimag(k),
     .   iish(1,1),iish(2,1),iish(1,2),iish(2,2)
100   format(i3,1x,i4,1x,f8.3,1x,e14.8,1x,e14.8,4(1x,i2))
c
      return
      end
ccc
      subroutine bb_align(tf,tstart,wbb,nfbb,nrec,iifull)
c
      implicit none
      integer*4 nfbb,nrec,iifull,jzr,jf
      complex*8 tf(nfbb,nrec),eiph
      real*8 tstart,phase,wbb(nfbb)
c
      do jf=1,nfbb
         phase=-wbb(jf)*tstart
         eiph=cmplx(cos(phase),sin(phase))
         do jzr=1,nrec
            tf(jf,jzr)=tf(jf,jzr)*eiph
         enddo
      enddo
c: If highest frequency is Nyquist, make sure imaginary part is zero:
      if(iifull .eq. 1) then
         do jzr=1,nrec
            tf(nfbb,jzr)=cmplx(real(tf(nfbb,jzr)),0.e0)
         enddo
      endif
c
      return
      end
ccc
      subroutine bb_out(nfbb,nmode,nzsr,jmo,knbb,phibb,dpsibb,nh_off)
c
      implicit none
      integer*4 nfbb,nmode,nzsr,jmo,nh_off,jf,jsr
      complex*16 knbb(nfbb,nmode)
      complex*8 phibb(nfbb,nzsr),dpsibb(nfbb,nzsr)
c
      write(33,rec=nh_off+jmo) (knbb(jf,jmo),jf=1,nfbb),
     .   ((phibb(jf,jsr),jf=1,nfbb),jsr=1,nzsr),
     .   ((dpsibb(jf,jsr),jf=1,nfbb),jsr=1,nzsr)
c
      return
      end
ccc
      subroutine bb_merge(k,knbb,nfbb,nmode,jfbb2,jmo,errdk2,iimrg,
     .   jm_mrg,k_mrg,nmerge,faxbb,iidiag)
c
c: Checks if this mode jmo has merged with any of previous modes found.
c
      implicit none
      integer*4 nfbb,nmode,jfbb2,jmo,jm_mrg,iimrg,nmerge,iidiag
      complex*16 knbb(nfbb,nmode),k,k_mrg
      real*8 magsq,kdif,errdk2
      real*4 faxbb(nfbb)
c
      do jm_mrg=jmo-1,1,-1
         kdif=magsq(knbb(jfbb2,jm_mrg) - k)
         if(kdif .le. errdk2) then
            if(iidiag .ge. 2) then
               print *,'Merged modes found: ',faxbb(jfbb2),jm_mrg,
     .            jmo,knbb(jfbb2,jm_mrg),k
               print *,'k = ',knbb(jfbb2+1,jm_mrg),knbb(jfbb2+1,jmo)
            endif
            iimrg=iimrg + 1
c: Mode found at previous (higher) frequency that was valid for sure:
            k_mrg=knbb(jfbb2+1,jm_mrg)
            if(iimrg .eq. 1) nmerge=nmerge + 1
            return
         endif
      enddo
      iimrg=0
c
      return
      end
c
      subroutine bb_fft_out
c
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 irec,jj0,jsrc,jrec,j,mx,nrd
      real*8 vgmax,tbuf
      real*4 rspace,freqs,dt1,sd
      complex*8 zzero
      character*16 trfext
      logical bintrf
      data zzero/(0.,0.)/
c
c      bintrf=.true.
      bintrf=.false.
c
c     rspace = range increment for receiver position
      if (nsrc.eq.1) then
         rspace=rkm(1)
      else
         rspace=(rkm(nsrc)-rkm(1))/float(nsrc-1)
      end if
      mx=nf1+nfbb-1
      dt1=1./fsbb
      freqs=(fmax-fmin)/2.+fmin
      sd=zsr(mzsrc(1))
      outroot(1:lout)='cprosim_out'
      outfile=outroot(1:lout)//'.dat'
      svp_title='TRF OUTPUT FROM SAGA'
c      write(6,*)'1st rec depth: ',zrec(1)
c      write(6,*)'last rec depth: ',zrec(nrec)
c      write(6,*)'Delta time: ',dt1
c      write(6,*)'Freq components: ',mx
c      write(6,*)'Center frequency: ',freqs
c      write(6,*)'Source depth: ',sd
c      write(6,*)'Range space: ',rspace
c      write(6,*)'Logical bintrf ',bintrf
c      write(6,*)'Before call to trfhead'
c      pause
      if(tilth) then
        call trfhead(outroot,svp_title,zrec(1),zrec(1),
     &    rkm(1),rspace,nfftbb,nf1,mx,dt1,freqs,sd,
     &    bintrf,1,nzs,nsrc)
      else
c        write(*,*)' calling trfhead nrec=',nrec
        call trfhead(outroot,svp_title,zrec(1),zrec(nrec),
     &    rkm(1),rspace,nfftbb,nf1,mx,dt1,freqs,sd,
     &    bintrf,nrec,nzs,nsrc)
      end if
c      write(6,*)'After call to trfhead'
c
      if(iifft .ne. 0) then
         tbuf=0
         vgmax= cfmin
c         do jsrc=1,nsrc
c            jj0=(jsrc-1)*nfbb*nrec
c            call bb_align(tf(jj0+1),nfbb,nrec,iifull)
c         enddo
c
         if (bintrf) then
            if(tilth) then
               do j=1,nfbb
                  do jsrc=1,nsrc
                     jj0=(jsrc-1)*nfbb*nrec+(jsrc-1)*nfbb
                     write(luttrf)real(tf(j+jj0)),
     .                    -imag(tf(j+jj0))
                  end do
               end do
            else
               do j=1,nfbb
                  do jsrc=1,nsrc
                     jj0=(jsrc-1)*nfbb*nrec
c     write(6,*)'No of nsrc: ',nsrc
                     do jrec=1,nrec
c     write(6,*)'No of nrec: ',nrec
c     write(6,*)'range: ',range(jsrc*jrec)
                        write(luttrf)real(tf(j+jj0)),-aimag(tf(j+jj0))
                        jj0=jj0 + nfbb
                     end do
                  end do
               end do
            end if
         else
            if(tilth) then
               do j=1,nfbb
                  do jsrc=1,nsrc
                     jj0=(jsrc-1)*nfbb*nrec+(jsrc-1)*nfbb
                     write(luttrf,*)real(tf(j+jj0)),
     .                    -imag(tf(j+jj0))
                  end do
               end do
            else
               do j=1,nfbb
                  do jsrc=1,nsrc
                     jj0=(jsrc-1)*nfbb*nrec
                     do jrec=1,nrec
                        write(luttrf,*)real(tf(j+jj0)),-aimag(tf(j+jj0))
                        jj0=jj0 + nfbb
                     end do
                  end do
               end do
            end if
         end if
c     
         close(luttrf)
      endif
c
      return
      end
c
      subroutine trfhead(filenm,title,rd,rdlow,r0,rspace,
     &  nx,lx,mx,dt1,freqs,sd,bintrf,ir,is,
     &  nplots)
c
      implicit none
c
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
c      include 'i_o_com'
c
      logical bintrf
      integer*4 i,idummy,inttyp,isrow,icdr,msuft,nx,lx,mx,
     &          nplots,ir,is,nout,icnt
      real*4 dummy,adummy(10),omegim,r0,rd,rspace,rdlow,dt1,
     &       sd,freqs
      character*6 prognm
      character*8 fileid
      character*64 filenm
      character*64 title
      character signn
c >>> Dummy real*4 to force single prec trf file.
c
      omegim=0.
      if (bintrf) then
         if(is.gt.1) stop'in bb_fft_out'
c         open(luttrf,file=filenm(1:lout)//'.dat',
c     &        status='unknown',form='unformatted')
         fileid='PULSETRF'
         write(luttrf) fileid
         prognm='ORCA'
         write(luttrf) prognm
c         nout=1
c         write(luttrf) nout
         write(luttrf) int(1.)
c         icnt=1
c         write(luttrf) icnt
         write(luttrf) int(1.)
         write(luttrf) title(1:64)
         signn='+'
         write(luttrf) signn
c     CENTER FREQUENCY
         dummy=freqs
         write(luttrf) dummy
c     SOURCE DEPTH
         dummy=sd
         write(luttrf) dummy
c     UPPER MOST RECEIVER DEPTH
         adummy(1)=rd
c     LOWER MOST RECEIVER DEPTH
         adummy(2)=rdlow
c     IR=NO OF RECEIVERS BETWEEN R0 AND RDLOW
         write(luttrf) adummy(1),adummy(2),ir
C
C     MAY BE ADDED IN ORCA
C     IF (IR.LT.0) THEN
c     do L = 1,abs(ir)
c     dummy=rdc(L)
c     WRITE(LUTTRF) dummy
c     end do
c     RDC(L) CONTAINS ALL THE RECEIVER DEPTHS
C     write(luttrf) (rdc(l),l=1,abs(ir))
C     END IF
C
C     R0= THE FIRST RANGE TO PLOT IN RANGE STACKED PLOT
         adummy(1)=r0
C     RSPACE= THE RANGE INCREMENT TO PLOT IN RANGE STACKED PLOT
         adummy(2)=rspace
c     WRITE R0, RSPACE AND THE NO OF PLOTS ASSOCIATED WITH IR
         WRITE(LUTTRF) adummy(1),adummy(2),nplots
C     DT=TIME SAMPLING
         dummy=dt1
C     NX=NO OF TIME SAMPLES (DENOTED NT IN OASES MANUAL) MAYBE?
C     LX=INDEX OF FIRST FREQUENCY COMPONENT (INT(FR1*DT))
C     MX=INDEX OF LAST FREQUENCY COMPONENT (INT(FR2*DT))
         write(luttrf) nx,lx,mx,dummy
         icdr=0
         write(luttrf) icdr
         dummy=omegim
         write(luttrf) dummy
C     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         write(luttrf) msuft
c         write(6,*) 'trfhead: msuft=',msuft
         isrow=1
         write(luttrf) isrow
         inttyp=1
         write(luttrf) inttyp
         idummy=0
         do 300 i=1,2
            write(luttrf) idummy
 300     continue
         dummy=0
         do 400 i=1,5
            write(luttrf) dummy
 400     continue
      else
c         write(6,*)'No of char in filenm: ',lout
c         write(6,'(a)')'Filename: ',filenm(1:lout)
c         pause
c         open(luttrf,file=filenm(1:lout)//'.dat',
c     &        status='unknown',form='formatted')
c         open(luttrf,file='pek.asc',
c     &        status='unknown',form='formatted')
         fileid='PULSETRF'
         prognm='ORCA'
         write(luttrf,'(1x,a)') fileid
         write(luttrf,'(1x,a)') prognm
c     WRITE(LUTTRF,*) NOUT
         write(luttrf,*) int(1.)
         icnt=0
         write(luttrf,*) int(1.)
         write(luttrf,'(1x,a)') title(1:64)
         signn='+'
         write(luttrf,'(1x,a)') signn
         write(luttrf,*) freqs
         write(luttrf,*) sd
c     IR=1
         write(luttrf,*) rd,rdlow,ir
C
C     MAY BE ADDED IN ORCA
C     IF (IR.LT.0) THEN
C     WRITE(LUTTRF,*) (RDC(L),L=1,ABS(IR))
C     END IF
C     NPLOTS=1
C
         write(luttrf,*) r0,rspace,nplots
         write(luttrf,*) nx,lx,mx,dt1
         icdr=0
         write(luttrf,*) icdr
         write(luttrf,*) omegim
c     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         write(luttrf,*) msuft
         isrow=1
         write(luttrf,*) isrow
         inttyp=1
         write(luttrf,*) inttyp
         idummy=0
         do 301 i=1,2
            write(luttrf,*) idummy
 301     continue
         dummy=0e0
         do 401 i=1,5
            write(luttrf,*) dummy
 401     continue
      end if
      return
      end
ccc
      subroutine bb_fft_out_original
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 irec,jj0,jsrc,jrec,j,lfft
      real*8 vgref,tbuf,tstart
      complex*8 zzero
      data zzero/(0.,0.)/
c
      if(iifft .ne. 0) then
c: Open output FFT file:
         call openfftout2(10,outroot,lout,fftfile,lfft,nfftbb,
     .      fsbb,fmin,fmax,xhbb,NRECL)
cc       open(10,file=outroot(1:lout)//'_fft',access='direct',
cc   .      recl=NRECL*(22+nfftbb))
         irec=0
c: Set up FFT header:
c: Flag that FFT is a simulation, not data:
         xhbb(1)=-1.
         xhbb(3)=1
         xhbb(10)=0.
         xhbb(11)=nrec
         xhbb(12)=1
c: Make ref time for impulse resp 10% of the way through total time window:
         tbuf=.10*Tw
c: Output transfer functions to the FFT file:
c: Make reference sound speed, used to shift impulse responses as a 
c: function of range, the minimum speed in the profile.
         vgref=cfmin
         do jsrc=1,nsrc
            tstart=rng_sr(jsrc)/vgref - tbuf
            jj0=(jsrc-1)*nfbb*nrec
            call bb_align(tf(jj0+1),tstart,wbb,nfbb,nrec,iifull)
            xhbb(14)=tstart
            do jrec=1,nrec
               xhbb(4)=jsrc
c: Not compatible with new FFT format:
cc             xhbb(15)=zsr(mzrec(jrec))
cc             xhbb(17)=cp_sr(mzrec(jrec))
               xhbb(18)=rng_sr((jrec-1)*nsrc + jsrc)
               xhbb(19)=zsr(mzrec(jrec))
               irec=irec + 1
c: Note that CONJUGATE of transfer function is output to FFT file:
               write(10,rec=irec) (xhbb(j),j=1,20),
     .            (conjg(tf(j)),j=jj0+1,jj0+nfbb)
               jj0=jj0 + nfbb
            enddo
         enddo
         close(10)
      endif
c
      return
      end
ccc
      subroutine bb_mfout(jmo,jfbb,nfout,phibbx,phi_re,phi_im)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
c
      integer*4 jmo,jfbb,nfout,jrec,jz,jf,jfx,jfile
      complex*8 phibbx(nfbb,nzsr)
      real*4 phi_re(nfout,nrec),phi_im(nfout,nrec)
      character*10 mrsuf,misuf
      data mrsuf/'phire_m000'/,misuf/'phiim_m000'/
c
      do jrec=1,nrec
         jz=mzrec(jrec)
         do jf=jfbb,nfbb
            jfx=jf-jfbb+1
            phi_re(jfx,jrec)=real(phibbx(jf,jz))
            phi_im(jfx,jrec)=aimag(phibbx(jf,jz))
         enddo
      enddo
c
      write(mrsuf(8:10),'(i3.3)') jmo
      write(misuf(8:10),'(i3.3)') jmo
      jfile=0
      if(mod(jmo,10) .eq. 1) jfile=2
      call out_writex(outroot,lout,SUFX//mrsuf,11,phi_re,zrec,
     .   faxbb(jfbb),nrec,nfout,dlab,flab,z4,z4,z4,z4,jfile,
     .   'Re(phi) vs z,f',' ',' ',' ','f7.2','f7.2','f8.5',ncall)
      call out_writex(outroot,lout,SUFX//misuf,11,phi_im,zrec,
     .   faxbb(jfbb),nrec,nfout,dlab,flab,z4,z4,z4,z4,jfile,
     .   'Im(phi) vs z,f',' ',' ',' ','f7.2','f7.2','f8.5',ncall)
c
      return
      end
ccc
      subroutine bb_ftinit
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
c: Open file for the frequency trajectories of the modes:
      if(iift .eq. 1) then
         open(14,file=outroot(1:lout)//'_ftraj',status='unknown',
     .        form='formatted')
         write(14,100) -1,faxbb(1),faxbb(nfbb),fsbb,nfftbb,kw0,0,0,0
100      format(i3,2x,3(f8.3,2x),i5,2x,e10.4,3(1x,i2))
c: Write branch point(s) to frequency trajectory file:
         call bb_write(0,nfbb,faxbb,xkbp(1,1),kw0,iish)
         call bb_write(0,nfbb,faxbb,xkbp(1,2),kw0,iish)
         if(geo(2,1,1) .gt. 500.) then
            call bb_write(0,nfbb,faxbb,xkbp(2,1),kw0,iish)
            call bb_write(0,nfbb,faxbb,xkbp(2,2),kw0,iish)
         endif
      endif
c
      return
      end
ccc
      subroutine bb_done(nm_fmax)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
      integer*4 nm_fmax,jsr,jf
c
c: Close files:
      if(iimt .ne. 0) close(21)
      if(iift .ne. 0) close(14)
      if(iiout .ne. 0) then
         close(33)
         open(33,file=outroot(1:lout)//'_bbeig',status='unknown',
     .      access='direct',
     .      recl=NRECL*lheadw)
         write(33,rec=1) fsbb,nfftbb,fmin,fmax,nfbb,nzsr,nm_fmax,
     .      rmin,rmax,sngl(cfmin),
     .      (zsr(jsr),jsr=1,nzsr),(nmbb(jf),jf=1,nfbb),iirx,
     .      (sngl(rho_sr(jsr)),jsr=1,nzsr)
cc       print *,'nmbb at end= ',(nmbb(jf),jf=1,nfbb)
cc       print *,'lheadw,lrecw = ',lheadw,lrecw,nh_off
         close(33)
      endif
c
c: Output FFT file:
c pln 020500       call bb_fft_out
c
c: Only output kni if not already done in disp_curv:
c pln 020500     if(iidc .ne. 2 .and. iidc .ne. 3) then
c         call out_writex(outroot,lout,SUFX//'kni',4,r4mat2,faxbb,xmode,
c     .        nfbb,nm_fmax,flab,mnlab,z4,z4,z4,z4,2,'Im(kn) vs Mode No',
c     .        ' ',' ','dB/km','f5.0','f6.1','f7.1',ncall)
c      endif
c
      return
      end
ccc
      subroutine bb_miss(knbbx,nm_fmax,knx,nm_fmin,kn_trk,nmiss)
c
c: Compares modes at fmin that were tracked from fmax with those that
c: were found at fmin using mode_find.  Finds any missing modes.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 nm_fmax,nm_fmin,nmiss,
     .   nm_trk,jm,jm_trk,nm1,jmf
      complex*16 knbbx(nfbb,nm_fmax),knx(nm_fmin),kn_trk(nm_fmax),
     .   k_found
      real*8 magsq
c
      nm_trk=0
      nm1=nm_fmax+1
      do jm=1,nm_fmax
c: Copy all non-zero eigenvalues to new list in knbb(:,nm1)
         if(knbbx(1,jm) .ne. (0.d0,0.d0)) then
            nm_trk=nm_trk + 1
            kn_trk(nm_trk)=knbbx(1,jm)
         endif
      enddo
c: Sort eigenvalues in knbb(:,nm1) for comparison with knx(1:nm_fmin:
      call hpsort_indx_c16(nm_trk,kn_trk,kn_indx)
      print *,'Sorted kn_trk: ',nm_trk,(kn_trk(jm),jm=1,nm_trk)
      print *,'kn at fmin: ',nm_fmin,(knx(jm),jm=1,nm_fmin)
c
      nmiss=0
      if(nm_trk .eq. nm_fmin) then
         print *,'BB tracking found same # modes as CW at fmin: ',
     .      nm_trk,nm_fmin
      elseif(nm_trk .gt. nm_fmin) then
         print *,'BB tracking found more modes than CW at fmin: ',
     .      nm_trk,nm_fmin
      else
         print *,'BB tracking found fewer modes than CW at fmin: ',
     .      nm_trk,nm_fmin
      endif
c
      jm_trk=1
      kn_trk(nm_trk+1)=(1.d100,0.d0)
      do jmf=1,nm_fmin
         k_found=knx(jmf)
10       continue
         if(magsq(k_found-kn_trk(jm_trk)) .gt. errdk2) then
            if(dreal(k_found) .gt. dreal(kn_trk(jm_trk))) then
               nmiss=nmiss + 1
               kn_indx(nmiss)=jmf
            else
               jm_trk=jm_trk + 1
               goto 10
            endif
         else
            jm_trk=jm_trk + 1
         endif
      enddo
      print *,'# missing modes: ',nmiss
      print *,'Missing modes: ',(knx(kn_indx(jm)),jm=1,nmiss)
c
      return
      end
      subroutine cub_fit_new(k1,k2,y1,y2,y1p,y2p,c1,c2,c3,c4,delk)
c
c: This subroutine fits a cubic polynomial y(x)=c1 + c2*x + c3*x**2 + 
c: c4*x**3, where x=(k-k1)/(k2-k1), to the points (k1,y1) and (k2,y2) 
c: and the derivatives y1p and y2p at those points.
c
      implicit none
      real*8 k1,k2,y1,y2,y1p,y2p,c1,c2,c3,c4,b1,b2,delk,y1px,y2px
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
c
      return
      end 
ccc
      subroutine cub_root_new(c1,c2,c3,c4,xhit,kcub)
c
c: Finds the root, xhit, between 0 and 1 of the
c: polynomial f(x)=c1 + c2*x + c3*x**2 + c4*x**4.
c
      implicit none
      integer*4 kcub
      real*8 c1,c2,c3,c4,xhit,a1,a2,a3,a1sq,q,r,
     .   dd,rdd,rq,arg,th,pietth,piefth,one_third
      data pietth/2.09439510239320/,piefth/4.18879020478639/,
     .   one_third/0.33333333333333/
c
      kcub=1
      a1=c3/c4
      a2=c2/c4
      a3=c1/c4
      a1sq=a1*a1
      q=(3.d0*a2 - a1sq)/9.
      r=(9.d0*a1*a2 - 27.d0*a3 - 2.d0*a1*a1sq)/54.d0
      dd=q**3 + r**2
      if(dd .ge. 0.d0) then
         rdd=sqrt(dd)
         xhit=dsign(1.d0,r+rdd)*abs(r+rdd)**(one_third) +
     .      dsign(1.d0,r-rdd)*abs(r-rdd)**(one_third) - a1*one_third
      else
         rq=sqrt(-q)
         arg=r/rq**3
         if((arg .lt. -1.d0) .or. (arg .gt. 1.d0)) then
            xhit=0.5d0
            kcub=0
            return
         endif
         th=acos(arg)/3.d0
         xhit=2.d0*rq*cos(th) - a1*one_third
         if((xhit .lt. 0.d0) .or. (xhit .gt. 1.d0)) then
            xhit=2.d0*rq*cos(th + pietth) - a1*one_third
            if((xhit .lt. 0.d0) .or. (xhit .gt. 1.d0)) then 
               xhit=2.d0*rq*cos(th + piefth) - a1*one_third
            endif
         endif
      endif
      if((xhit .le. 0.d0) .or. (xhit .ge. 1.d0)) then
         xhit=.5d0
         kcub=0
      endif
c     if(kcub .eq. 0) print *,'kcub=0 in cubroot' 
c
      return
      end 
ccc
