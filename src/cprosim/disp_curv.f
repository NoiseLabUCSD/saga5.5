      subroutine disp_curv(iidc,knbb,eig_bb,nfbb,nmode,faxbb,xmode,
     .   r4mat1,r4mat2,outroot,lout,ncall,SUFX)
c
      implicit none
      integer*4 iidc,nfbb,nmode,lout,ncall,jm,jf
      complex*16 knbb(nfbb,nmode),eig_bb(5,nfbb,nmode)
      real*4 xmode(nmode),r4mat1(nmode,nfbb),r4mat2(nmode,nfbb),
     .   faxbb(nfbb),twpie
      character*16 outroot
      character*1 SUFX
      include 'lab_com'
      data twpie/6.283185307/
c
      do jm=1,nmode
         xmode(jm)=jm
      enddo
c
      if(iidc .eq. 1 .or. iidc .eq. 3) then
         do jf=1,nfbb
            do jm=1,nmode
               r4mat1(jm,jf)=real(eig_bb(4,jf,jm))
            enddo
         enddo
         write(6,*)'CALLING OUT_WRITEX FOR VG AND VPH'
         call out_writex(outroot,lout,SUFX//'vg',3,r4mat1,faxbb,xmode,
     .      nfbb,nmode,flab,mnlab,z4,z4,z4,z4,2,'Vg vs Mode No',
     .      ' ',' ','m/s','f7.2','f5.0','f8.2',ncall)
      endif
c
      if(iidc .eq. 2 .or. iidc .eq. 3) then
         do jf=1,nfbb
            do jm=1,nmode
               if(real(knbb(jf,jm)) .ne. 0.d0) then
                  r4mat1(jm,jf)=twpie*faxbb(jf)/real(knbb(jf,jm))
                  r4mat2(jm,jf)=dimag(knbb(jf,jm))
               else
                  r4mat1(jm,jf)=0.e0
                  r4mat2(jm,jf)=0.e0
               endif
            enddo
         enddo
c
         call out_writex(outroot,lout,SUFX//'vph',4,r4mat1,faxbb,xmode,
     .      nfbb,nmode,flab,mnlab,z4,z4,z4,z4,2,'Vph vs Mode No',
     .      ' ',' ','m/s','f7.2','f5.0','f8.2',ncall)
         call out_writex(outroot,lout,SUFX//'kni',4,r4mat2,faxbb,xmode,
     .      nfbb,nmode,flab,mnlab,z4,z4,z4,z4,2,'kni vs Mode No',
     .      ' ',' ','m/s','f5.0','f7.2','f8.2',ncall)
      endif
c
      return
      end
ccc
      subroutine disp_fill(nmode,nfcw,jfcw,kn,eig_char,knbb,eig_bb,
     .   nm_cw_max,fcw)
c
c: Fills eig_bb(4,jfcw,1:nmode) with the group velocities in 
c: eig_char(4,1:nmode) for the dispersion curve output option for CW
c: mode computations.
c
      implicit none
      include 'Parms_com'
      integer*4 nmode,nfcw,jfcw,nm_cw_max,jf,jm,nmm,iibad
      real*4 fcw(nfcw),fm
      complex*16 kn(0:nmode),eig_char(5,nmode),knbb(nfcw,nmode),
     .   eig_bb(5,nfcw,nmode)
c
c: If first CW frequency, initialize the arrays:
      if(jfcw .eq. 1) then
         fm=fcw(1)
         do jf=2,nfcw
            fm=max(fm,fcw(jf))
         enddo
c: Take a conservative guess at the maximum # of modes at max CW frequency:
         nmm=nint((fm/fcw(jfcw))*(nmode+2)) + 10
         call mem_lim(nfcw*nmm,NM_NF_MAX,MLINE,LML,'nfcw*nmode_max',
     .      16,'NM_NF_MAX',9,iibad,1)
         do jf=1,nfcw
cpln            do jm=1,nmm
            do jm=1,nmode
               eig_bb(4,jf,jm)=dcmplx(-999.d0,-999.d0)
               knbb(jf,jm)=dcmplx(1.d100,0.d0)
            enddo
         enddo
         nm_cw_max=nmode
      endif
c
      nm_cw_max=max(nm_cw_max,nmode)
      do jm=1,nmode
         eig_bb(4,jfcw,jm)=eig_char(4,jm)
         knbb(jfcw,jm)=kn(jm)
      enddo
c
      return
      end
