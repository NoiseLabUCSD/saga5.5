c: *******************************
c: *   AUTHOR (1996):            *
c: *     P R O S I M  GROUP      *
c: *     S A C L A N T C E N     *
c: *******************************
      subroutine tenv_read_tilt(psect,intrpl,isubx)
c     
      implicit none
      include 'Parms_com'
      include 'i_o_svp_com'
      include 'i_o_opt_com'
      include 'i_o_1b_com'
      include 'sector_env_com'
      common /out_com3/ dmmy1,dmmy2,NSCTOR,dmmy3
      integer*4 dmmy1,dmmy2,NSCTOR,dmmy3,isum,ii
      include 'depth_com'
c TEMP PLN
      common /out_com2b/ faxbb(NFBBMAX),
     .   mode_phz(3,0:NM_MAX),iigeom,
     .   nfftbb,nfft0
      integer*4 nfftbb,iigeom,nfft0,isctr
      real*4 faxbb,mode_phz,seclngtot

c:    Local variables
      real*8 zsx1(2*NVRMAX),zsx2(2*NVRMAX),tan_a1,tan_ax,zra
      real*8 a1,ax,hhb
      integer*4 i,j,j1,psect,isubx,lu_env,ienv_info
      integer*4 index,intrpl(RGMAX),displcm
      common /info/ienv_info,lu_env
      integer*4 iter,iforwpop,ipop
      common/iterpar/iter,iforwpop,ipop
c
      common /tiltparm/tilt,dtilt
      logical tilt
      real dtilt
c
      j1=1
      geot(2,1,j1)= 343.0
c      geot(2,1,j1)= 1.0
      geot(2,2,j1)= 0.0
c      geot(2,3,j1)= 1.0e-5
c      geot(2,3,j1)= 1.0e-8
      geot(2,3,j1)= 1.21e-3
      geot(2,4,j1)= 0.0
      geot(2,5,j1)= 0.0
      do j=1,5
         geot(1,j,j1)= geot(2,j,j1)
      enddo
      nsvp=R_ND0(psect)

      do j=1,nsvp
         zsvp(j)=R_Z0(j,psect)
         csvp(j)=R_C0(j,psect)
      end do
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
      do j=1,nlayb
         ktb(j)= 1.0
         hb(j)=R_Z1(j+1,psect)-R_Z1(j,psect)
cpln         hb(j)=R_Z1(j+1,psect)
         geob(1,1,j)= R_C1(j,psect) ! P sound speed
cpln Original
        geob(2,1,j)= R_C1(j+1,psect) !P sound speed
cpln ONLY ISO SPEED LAYERING IN BOTTOM
cpln         geob(2,1,j)= R_C1(j,psect) !P sound speed
cpln ONLY for IVWK
cpln         geob(2,1,j)= R_C1(j,psect) !P sound speed
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
cONLY for IVWK
cpln         geob(1,4,j)= R_BETA(j,psect)
cpln         geob(2,4,j)= R_BETA(j,psect)
         geob(1,4,j)= R_BETA(1,psect)
         geob(2,4,j)= R_BETA(1,psect)
         geob(1,5,j)= 0.0       ! shear att
         geob(2,5,j)= 0.0       ! shear att
      end do
c
c     read subbottom parameters from temp file LUSVP
      j1=nlayb+1
      
      geob(1,1,j1)= R_C2(psect) ! P sound speed
      geob(2,1,j1)= R_C2(psect) !P sound speed
      geob(1,2,j1)= 0.0         ! shear speed
      geob(2,2,j1)= 0.0         ! shear speed
      geob(1,3,j1)= R_R2(psect)
      geob(2,3,j1)= R_R2(psect)
cONLY for IVWK
cpln      geob(1,4,j1)= R_BETA(j1,psect)
cpln      geob(2,4,j1)= R_BETA(j1,psect)
      geob(1,4,j1)= R_BETA(2,psect)
      geob(2,4,j1)= R_BETA(2,psect)
      geob(1,5,j1)= 0.0         ! shear att
      geob(2,5,j1)= 0.0         ! shear att
c
cpln      do ii=1,j1
cpln         write(6,*)'R_BETA: ',R_BETA(ii,psect)
cpln      end do
cpln      pause
c
      rmin=rkm(1)
c  added for SAGA test 22/11/99
      rmin=1.e-20
c
      if(nrecusr .gt. 1) then
         if(tilt) then
            do i=1,nrecusr
               dtiltp(i)=dtilt/(nrecusr-1)*(i-1)
            end do
         else
            do i=1,nrecusr
               dtiltp(i)=0.0e0
            end do
         end if
      else
         tilt=.false.
      end if
c
c
      if((iiwrite.gt.0).and.(isubx.eq.1)) then

C WRITE OUT ENVIRONMENT TO FILE LU_ENV
         write(lu_env,*)'Population: ',ipop,
     &                  '  forward model: ',iforwpop

         write(lu_env,*)'***** SECTOR NO',psect
         write(lu_env,*)'First layer (top)'
         write(lu_env,*)(geot(1,j,1),j=1,5)
         write(lu_env,*)(geot(2,j,1),j=1,5)
c     
         write(lu_env,*)
         write(lu_env,*)'Water column'
         do j=1,nsvp
            write(lu_env,*)zsvp(j),csvp(j)
         end do
c     
         if (R_ND1(psect).gt.0) then
            write(lu_env,*)
            write(lu_env,*)'Sediment'
            do j1=1,nlayb
               write(lu_env,*)hb(j1),(geob(1,j,j1),j=1,5)
               write(lu_env,*)hb(j1),(geob(2,j,j1),j=1,5)
            end do
         end if
c
         write(lu_env,*)
         write(lu_env,*)'Sub-bottom'
         write(lu_env,*)(geob(1,j,j1),j=1,5)
         write(lu_env,*)(geob(2,j,j1),j=1,5)
        
         write(lu_env,*)
         write(lu_env,*)' First ranges:',  rkm(1)
         write(lu_env,*)' Source depth', zsrc(1)
         write(lu_env,*)' First receiver depth',  zrecusr(1)
         write(lu_env,*)
         write(lu_env,*)' Tilt of array:',  dtilt
         write(lu_env,*)
         write(lu_env,*)

C WRITE OUT ENVIRONMENT TO INPUT FILE FOR STAND ALONE PROSIM
         if(psect .eq. 1) then
            write(lualon,*)'OUTPUT FROM SAGA'
            write(lualon,'(I5)')0
            write(lualon,'(2I5)')500,500
            write(lualon,'(3I5)')1,500,500
            write(lualon,'(2I5)')0,0
         end if

         seclngtot=0
         do isctr=1,psect
            seclngtot=seclngtot+secleng(psect)
         end do

         write(lualon,'(F15.8,2I5,F10.3,I5)')zsvp(nsvp),
     .        0,0,seclngtot/1000.0,0
           
         do j=1,nsvp
            write(lualon,'(2F15.8)')zsvp(j),csvp(j)
         end do
         if (R_ND1(psect).gt.0) then
            hhb=0
            do j1=1,nlayb
               hhb=hhb+hb(j1)
            end do
            write(lualon,'(2F15.8)')hhb,geob(1,3,nlayb)
cpln     .            geob(1,4,nlayb)
            hhb=0
            write(lualon,'(3F15.8)')hhb,geob(1,1,1),geob(1,4,1)
            do j1=1,nlayb-1
               hhb=hhb+hb(j1)
               write(lualon,'(3F15.8)')hhb,geob(1,1,j1+1),
     .               geob(1,4,j1+1)
            end do
            j1=nlayb
            hhb=hhb+hb(j1)
            write(lualon,'(3F15.8)')hhb,geob(1,1,j1),
     .               geob(1,4,j1)
            j1=j1+1
            write(lualon,'(3F15.8)')geob(1,3,j1),geob(1,4,j1),
     .            geob(1,1,j1)
            write(lualon,'(2I5)')0,0
            if(psect .eq. nsctor) then
               write(lualon,*)'PULSE'
               write(lualon,'(2I5)')0,0
               write(lualon,'(1I5,3F10.1)')nfftbb,fmin,fmax,-fsbb
               write(lualon,*)rkm(1)*1.E-3,rkm(nsrc)*1.E-3,nsrc
               write(lualon,*)zrecusr(1),zrecusr(nrecusr),
     .              nrecusr
               write(lualon,*)zsrc(1)
            end if
         end if


C
      endif
      nlayt=0
c
      isum=isubx+psect
      if (NSCTOR .gt. 1) then
         if(psect .eq. 1) then
            index = 1
            nrecold = 0
         else
            index = intrpl(psect-1) + 1
         endif
         if(index .eq. 0) index=1
         
         displcm = 0
         if (intrpl(psect) .ne. 0) then
            tan_a1 = (secz(psect)-secz(psect+1))/secleng(psect)
cpln            a1 = datan(tan_a1)
            do j=index, intrpl(psect)
               if(psect .eq. 1) then
                  if(rkm(j) .le. secleng(psect))
     &                 zra = secz(psect) - tan_a1*rkm(j)
               else
                  zra = secz(psect) - tan_a1*(rkm(j)-
     &                 (rstart+secleng(psect)))
               end if
               do i=1, nrecusr
                  ax=(zrecusr(i)*tan_a1)/zra
cpln                  tan_ax=dtan(ax)
cpln zsx1 is virtual receivers to left of sector
                  zsx1(i+displcm)=zrecusr(i)+
     &                 (ax*(rkm(j)-rstart))
cpln zsx1 is virtual receivers to right of sector
                  zsx2(i+displcm)=zsx1(i+displcm)-
     &                 (ax*secleng(psect))
               enddo
               displcm=displcm+nrecusr
            enddo
            do j=1, nrecold
               zrec(j) = zsx2old(j)
            enddo
cpln            write(6,*)'displcm,index,nrec,nrecusr: ',displcm,
cpln     .                 index,nrec,nrecusr,intrpl(psect)
cpln            write(6,*)
cpln            pause
            do j=1,displcm
               zrec(nrecold+j)=zsx1(j)
               zsx2old(j)=zsx2(j)
            enddo
            nrec=nrecold+displcm
            nrecold = displcm
         elseif(isum .eq. 2) then
            zrec(1)=zrecusr(1)
            zrec(2)=zrecusr(nrecusr)
            call uni_space(-nzs,zsrc,1.e0)
            call uni_space(-nrecusr,zrec,1.e0)
cpln            write(6,*)'RD environment'
         end if
      else
         if(isum .eq. 2) then
            zrec(1)=zrecusr(1)
cpln As in binigauss:
cpln            zrec(2)=zrecusr(nrecusr)
            zrec(2)=zrecusr(nrec)
            call uni_space(-nzs,zsrc,1.e0)
cpln As in binigauss:
            call uni_space(-nrec,zrec,1.e0)
cpln            call uni_space(-nrecusr,zrec,1.e0)
cpln            write(6,*)'No RD environment'
         end if
      end if
         
      if((iiwrite.gt.0).and.(isubx.eq.1)) then
cpln      if(isubx.eq.1) then
         write(lu_env,*)
         write(lu_env,*)'ZRECUSR'
         write(lu_env,*)(zrecusr(j),j=1,nrecusr)
         write(lu_env,*)
         write(lu_env,*)'ZREC'
         write(lu_env,*)(zrec(j),j=1,nrec)
      end if

      return
      end

