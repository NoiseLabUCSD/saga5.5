       subroutine bb_brute
c
c: Computes the broadband field from fmin to fmax in steps of df
c: using brute force calls of the usual CW routines
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
      include 'sector_env_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
      integer*4 jm,jfbb,nmbbtot,nm_fmax,nctot0,ncjf,jzjm,j
      real*8 watscale 
      real*4 ncperm
c
c      write(*,*)'hello'
      nm_put=0
      nmode=0
      nctot0=0
      nctot=0
      ncperm=0
      nm_tot=0
      nm_fmax=0
c
      call bb_init
      call zmx_init
      f_hz=faxbb(nfbb)
      call freq_init
      call bb_ftinit
c
      nmbbtot=0
      if(i_call_porter .eq. 1) then
         if(iiwrite. gt. 0)
     .        write(6,*)'Full-blown calculation'
         do jfbb=nfbb,1,-1
            nctot0=nctot
            f_hz=faxbb(jfbb)
            call freq_chng
            call mode_find(0)
            if(jjfail.gt. 0) return
            if(jfbb .eq. nfbb) then
               if(iiAih(1) .eq. -1 .and. iiAih(2) .eq. -1) then
                  nm_fmax=nm_tot
               else
                  nm_fmax=max(1.1*nm_tot,float(nm_tot) + 12)
               endif
               if(iiout .ne. 0) call bb_out_init(nm_fmax)
            endif
            do jm=1,nm_tot
               call bb_field(kn(jm),phi,dphi,dpsi,exp_gbs,
     .                            jm,tf,jfbb,jm)
               call bb_enter(kn(jm),eig_char(1,jm),
     .              eig_char(4,jm),eig_char(5,jm),knbb,eig_bb,
     .              jfbb,jm,nfbb,nm_tot,iish,iish_ref,iish_bb,nmbb)
               if(iift .eq. 1) call bb_write(jm,jfbb,faxbb,kn(jm),
     .                                       kw0,iish)
            enddo
            nmbb(jfbb)=min(nm_fmax,nm_tot)
            if(iiout .ne. 0) then
               if (nm_tot .gt. nm_fmax) print *,'nm_tot>nm_fmax',
     .              f_hz,nm_tot,nm_fmax
               write(33,rec=nh_off+jfbb) (kn(jm),jm=1,nmbb(jfbb)),
     .              (phi(jzjm),jzjm=1,nzsr*nmbb(jfbb))
            endif
            nmbbtot=nmbbtot + nm_tot
            ncjf=nctot - nctot0
            if(iiwrite .gt. 0)
     .           print *,'Done f,nm_tot,#R1R2/mode = ',
     .           sngl(f_hz),nm_tot,float(ncjf)/max(1,nm_tot)
            if(i_geom.eq.0) then
               do jzjm=1,nzsr*nmbb(jfbb)
                  phi_reuse(jfbb,jzjm)=phi(jzjm)
               end do
            end if
         enddo
      else
         if(iiwrite .gt. 0)
     .        write(6,*)'Reuse eigenvalues and mode functions'
c: pln
c: scale mode functions when searching for water depth
         if(r_h0(1).ne.iwat0) then
            watscale=r_h0(1)/iwat0
            do j=1,nzsrgeom
               zsrgeom(j)=zsrgeom(j)*watscale
            end do
            if(iiwrite.gt.0) then
               write(6,*)'New water depth:  ',r_h0(1)
               write(6,*)'Ref. water depth: ',iwat0
               write(6,*)
            end if
         end if
         do jfbb=nfbb,1,-1
            do jzjm=1,nzsrgeom*nmbb(jfbb)
               phi(jzjm)=phi_reuse(jfbb,jzjm)
            end do
            do jm=1,nmbb(jfbb)
               call bb_field_reuse(knbb((jm-1)*nfbb+jfbb),
     .              phi,exp_gbs,jm,tf,jfbb)
            end do
            nmbbtot=nmbbtot + nm_tot
            ncjf=nctot - nctot0
            if(iiwrite .gt. 0)
     .           print *,'Done f,nm_tot,#R1R2/mode = ',
     .           faxbb(jfbb),nm_tot,float(ncjf)/max(1,nm_tot)
         end do
c: pln
c: scale mode functions back to original for water depth
         if(r_h0(1).ne.iwat0) then
            watscale=iwat0/r_h0(1)
            do j=1,nzsrgeom
               zsrgeom(j)=zsrgeom(j)*watscale
            end do
         end if
      end if
c
      call bb_done(nm_fmax)
c
c: Output FFT file:
      
c no time for this Peter gerstoft      
c      call bb_fft_out
c
      ncperm=float(nctot)/float(max(1,nmbbtot))
      write(lusvp,120) nmbbtot,nctot,ncperm
120   format('CW LOOP # MODES = ',i8,'; # R1R2 CALCS = ',i8,
     .   '; #CALCS/MODE = ',f5.2)
c
      return
      end
