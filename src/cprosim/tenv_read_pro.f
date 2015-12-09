c: *******************************
c: *   AUTHOR (1996):            *
c: *     P R O S I M  GROUP      *
c: *     S A C L A N T C E N     *
c: *******************************
      subroutine tenv_read_pro(psect,isubx)
c     
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
c      include 'i_o_com'
      include 'sector_env_com'
      include 'debug_com'
c:    Local variables
      integer*4 i,j,j1,psect,isubx,iierr,nline
      character*64 eline
      common /info/ienv_info,lu_env
      common/iterpar/iter,iforwpop,ipop
      integer iter, iforwpop, ipop
      integer lu_env,ienv_info,ii
      real*8 hhb
      common /tiltparm/tiltv,tilth,dtiltv,dtilth
      logical tiltv,tilth
      real dtiltv,dtilth
c
      data eline/'INVALID INPUT IN TENV FILE: '/
c
      nline=0
      j1=1
      geot(2,1,j1)= 343.0
      geot(2,2,j1)= 0.0
      geot(2,3,j1)= 1.21e-3
      geot(2,4,j1)= 0.0
      geot(2,5,j1)= 0.0
      ht(j1)=1.d+20
      ktt(j1)= 1
      do j=1,5
         geot(1,j,j1)= geot(2,j,j1)
      enddo
c
      call svp_check_val_lay(2,0.d0,geot(1,1,j1),'upper h-space',13,
     .      eline,nline,iierr)
c
      nsvp=R_ND0(psect)
c
      if(nrec .gt. 1) then
         if(tiltv) then
            do i=1,nrec
               dtiltvp(i)=dtiltv/(nrec-1)*(i-1)
            end do
         else
            do i=1,nrec
               dtiltvp(i)=0.0e0
            end do
         end if
c      else
c         tiltv=.false.
      end if
c
      if(nsrc .gt. 1) then
         if(tilth) then
            dtilthp=dtilth/float(nsrc-1)
            do i=2,nsrc
               zrec(i)=zrec(1)+dtilthp*(i-1)
            end do
c            write(6,*)
c            write(6,*)'Horiz tilt: ',dtilth,zrec(1),zrec(nsrc)
         else
            do i=1,nsrc
               dtilthp=0.0e0
            end do
         end if
c      else
c         tilth=.false.
      end if
c
      if(iidebug .eq. 1)
     . write(6,*)'No of layers in water: ',nsvp

      do j=1,nsvp
         zsvp(j)=R_Z0(j,psect)
         csvp(j)=R_C0(j,psect)
         if(iidebug .eq. 1)
     .        write(6,*)zsvp(j),csvp(j)
      end do
c
      do j=2,nsvp
         call check_val_r8(zsvp(j),zsvp(j-1),1.d+10,eline,27,nline,
     .      'zsvp(j)',7,iierr)
         call check_val_r8(csvp(j),1.d-10,1.d+10,eline,27,nline,
     .      'csvp(j)',7,iierr)
      enddo
c
c This is ONLY for IVWK
cpln      zsvp(nsvp)=r_h0(psect)
cpln      csvp(nsvp)=1495.-0.04*zsvp(nsvp)
cpln      R_C0(nsvp,psect)=csvp(nsvp)
c
      rho_svp=1.0
      alpha_svp=0.0
c
      if (R_ND1(psect).eq.0) then
         nlayb=0
      else
         nlayb=R_ND1(psect)-1
      endif
      if(iidebug .eq. 1)
     . write(6,*)'No of layers in bottom: ',nlayb
c
      do j=1,nlayb
         ktb(j)= 1
         hb(j)=R_Z1(j+1,psect)-R_Z1(j,psect)
cpln     SAGA RETURNS LAYER THICKNESS AND NOT ABSOLUTE
cpln     SEDIMENT DEPTH RE WATER-SEDIMENT INTERFACE
cpln         hb(j)=R_Z1(j+1,psect)
         geob(1,1,j)= R_C1(j,psect) ! P sound speed


cpln Original
         geob(2,1,j)= R_C1(j+1,psect) !P sound speed
cpln ONLY ISO SPEED LAYERING IN BOTTOM:
cpln         geob(2,1,j)= R_C1(j,psect) !P sound speed


c Only for IVWK workshop
cpln         geob(2,1,j)= R_C1(j,psect) !P sound speed
c
cpln         write(6,*)hb(j),geob(1,1,j),geob(2,1,j)
cpln         write(6,*)hb(j),R_C1(j,psect),R_C1(j+1,psect)
cpln ONLY ASCOT01
cpln         if(nlayb.gt.2) then
cpln            geob(2,1,j)= R_C1(j+1,psect) !P sound speed
cpln         else
cpln            geob(2,1,j)= R_C1(j,psect) !P sound speed
cpln         end if
         geob(1,2,j)= 0.0       ! shear speed
         geob(2,2,j)= 0.0       ! shear speed
         geob(1,3,j)= R_R1(psect)
         geob(2,3,j)= R_R1(psect)
C ONLY FOR ITWK
cpln         geob(1,4,j)= R_BETA(j,psect)
cpln         geob(2,4,j)= R_BETA(j,psect)
         geob(1,4,j)= R_BETA(1,psect)
         geob(2,4,j)= R_BETA(1,psect)
         geob(1,5,j)= 0.0       ! shear att
         geob(2,5,j)= 0.0       ! shear att
         if(iidebug .eq. 1) then
            write(6,*)hb(j),geob(1,1,j),geob(2,1,j)
            write(6,*)geob(1,2,j),geob(2,2,j)
            write(6,*)geob(1,3,j),geob(2,3,j)
            write(6,*)geob(1,4,j),geob(2,4,j)
            write(6,*)geob(1,5,j),geob(2,5,j)
         end if
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
         call svp_check_val_lay(1,hb(j),geob(1,1,j),'bottom layer',12,
     .      eline,nline,iierr)
         call svp_check_val_lay(2,hb(j),geob(1,1,j),'bottom layer',12,
     .      eline,nline,iierr)
         call zero_sh(geob(1,2,j),geob(1,5,j),j,'top   ','bottom')
         call zero_sh(geob(2,2,j),geob(2,5,j),j,'bottom','bottom')
      end do
cpln      pause
c
c     read subbottom parameters from temp file LUSVP
      j1=nlayb+1
      
      hb(j1)=1.d+20
      ktb(j1)=1
      geob(1,1,j1)= R_C2(psect) ! P sound speed
      geob(2,1,j1)= R_C2(psect) !P sound speed
      geob(1,2,j1)= 0.0         ! shear speed
      geob(2,2,j1)= 0.0         ! shear speed
      geob(1,3,j1)= R_R2(psect)
      geob(2,3,j1)= R_R2(psect)
C ONLY ITWK
cpln      geob(1,4,j1)= R_BETA(j1,psect)
cpln      geob(2,4,j1)= R_BETA(j1,psect)
      geob(1,4,j1)= R_BETA(2,psect)
      geob(2,4,j1)= R_BETA(2,psect)
      geob(1,5,j1)= 0.0         ! shear att
      geob(2,5,j1)= 0.0         ! shear att

cpln      do ii=1,j1
cpln         write(6,*)'R_BETA: ',R_BETA(ii,psect)
cpln      end do
cpln      pause
c
c: c_hsp < 0 means use previous layer as halfspace:
c: Set thickness of halfspaces to large numbers (for use in zmx_init):
      call svp_check_val_lay(1,0.d0,geob(1,1,j1),'lower h-space',
     .     13,eline,nline,iierr)
      call zero_sh(geob(1,2,j1),geob(1,5,j1),j1,'top   ','bottom')
c
      if(iidebug .eq. 1) then
         write(6,*)'Sub-bottom'
         write(6,*)geob(1,1,j1),geob(2,1,j1)
         write(6,*)geob(1,2,j1),geob(2,2,j1)
         write(6,*)geob(1,3,j1),geob(2,3,j1)
         write(6,*)geob(1,4,j1),geob(2,4,j1)
         write(6,*)geob(1,5,j1),geob(2,5,j1)
      endif
c
      if(iiwrite .gt. 0) then
c     
C WRITE OUT ENVIRONMENT TO INPUT FILE FOR STAND ALONE PROSIM
         write(lualon,*)'OUTPUT FROM SAGA'
         write(lualon,'(I5)')0
         write(lualon,'(2I5)')0,0
         write(lualon,'(3I5)')1,0,0
         write(lualon,'(2I5)')0,0
         write(lualon,'(F33.26,2I5,F5.1,I5)')zsvp(nsvp),
     .                 0,0,10.0,0
         do j=1,nsvp
            write(lualon,'(2F33.26)')zsvp(j),csvp(j)
         end do
         hhb=0
         if (R_ND1(psect).gt.0) then
            do j1=1,nlayb
               hhb=hhb+hb(j1)
            end do
            write(lualon,'(3F33.26)')hhb,geob(1,3,nlayb),
     .            geob(1,4,nlayb)
            hhb=0
            write(lualon,'(2F33.26)')hhb,geob(1,1,1)
            do j1=1,nlayb
               hhb=hhb+hb(j1)
               write(lualon,'(2F33.26)')hhb,geob(2,1,j1)
            end do
         else
            write(lualon,'(3F33.26)')hhb,geob(1,3,nlayb+1),
     .            geob(1,4,nlayb+1)
         end if
         j1=nlayb+1
         write(lualon,'(3F33.26)')geob(1,3,j1),geob(1,4,j1),
     .        geob(1,1,j1)
         write(lualon,'(2I5)')0,0
         write(lualon,*)'PULSE'
         write(lualon,'(1I5,3F33.26)')nfftbb,fmin,fmax,-fsbb
         write(lualon,'(2F33.26,i5)')rkm(1)*1.D-3,rkm(nsrc)*1.D-3,nsrc
         write(lualon,'(2F33.26,i5)')zrec(1),zrec(nrec),
     .        nrec
         write(lualon,'(F33.26)')zsrc(1)
         close(lualon)

      end if
c
      nlayt=0
c
      return
      end
c
      subroutine zero_sh(cs,as,j,ch_tb1,ch_tb2)
c
      implicit none
      include 'Parms_com'
      real*8 cs,as
      integer*4 j
      character*6 ch_tb1,ch_tb2
c
      if(cs .eq. 0.d0 .and. as .ne. 0.d0) then
         write(lusvp,125) ch_tb1,ch_tb2,j
125      format('/*** SHEAR ATTENUATION SET TO 0 for 0 shear speed at ',
     .         a6,' of ',a6,' layer # ',i3,'***'/)
         as=0.d0
      endif
c
      return
      end
c
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
