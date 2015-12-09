      subroutine mfun_fill(phi_,mfun,mfunph,nmm,ir,jrun,nrunx,
     .   jset0,iirix,iimpx,jfcw,iidphi)
c
c: Outputs Re/Im and/or Mag/Phase of mode function phi_.
c
      use parms_com
      use i_o_com
      use gen_com
      implicit none
      integer*4 nmm,ir,jrun,nrunx,jset0,iirix,iimpx,jfcw,
     .   jzmf,jm,jsr,iidphi,nc
      complex*8 phi_(nzsr,nmode)
      real*4 mfun(nzmf,nmm),mfunph(nzmf,nmm),piedeg
      character*8 re,im
c
      piedeg=acos(-1.)/180.
      do jm=1,nmode
         xmode(jm)=jm
      enddo
c
      if(iidphi .eq. 0) then
         re(1:7)='Re(phi)'
         im(1:7)='Im(phi)'
         nc=7
      elseif(iidphi .eq. 1) then
         re(1:8)='Re(dphi)'
         im(1:8)='Im(dphi)'
         nc=8
      elseif(iidphi .eq. 2) then
         re(1:7)='Re(psi)'
         im(1:7)='Im(psi)'
         nc=7
      endif
c
cc      print *,'mfun_fill: ',iirix,iimpx,nzmf,nmode,nmm,NaN     
      if(iirix .ne. 0) then
         do jzmf=1,nzmf
            jsr=mzmf(jzmf)
            do jm=1,nmode
               mfun(jzmf,jm)=real(phi_(jsr,jm))
               mfunph(jzmf,jm)=aimag(phi_(jsr,jm))
            enddo
            do jm=nmode+1,nmm
               mfun(jzmf,jm)=NaN
               mfunph(jzmf,jm)=NaN
            enddo
         enddo
c      
cpg         if(iirix .eq. 1 .or. iirix .eq. 3) then
cpg            call hdf_write_gen(3+ir,outroot,lout,zmf,nzmf,xmode,nmm,
cpg     .         fcw,nfcw,var_ax,nrunx,1.,1,mfun,1,nzmf,1,nmm,
cpg     .         jfcw,jfcw,jrun,jrun,1,1,'Depth - m',9,'Mode',4,
cpgcpg     .         'Frequency - Hz',14,'Parameter',9,' ',1,re,
cpg     .         nc,jset0,1,2)
cpg         endif
cpg         if(iirix .eq. 2 .or. iirix .eq. 3) then
cpgcpg            call hdf_write_gen(3+ir,outroot,lout,zmf,nzmf,xmode,nmm,
cpg     .         fcw,nfcw,var_ax,nrunx,1.,1,mfunph,1,nzmf,1,nmm,
cpgcpg     .         jfcw,jfcw,jrun,jrun,1,1,'Depth - m',9,'Mode',4,
cpg     .         'Frequency - Hz',14,'Parameter',9,' ',1,im,
cpg     .         nc,jset0+1,1,2)
cpg         endif
      endif
      if(iimpx .ne. 0) then
         do jzmf=1,nzmf
            jsr=mzmf(jzmf)
            do jm=1,nmode
               mfun(jzmf,jm)=abs(phi_(jsr,jm))
               mfunph(jzmf,jm)=atan2(aimag(phi_(jsr,jm)),
     .            real(phi_(jsr,jm)))/piedeg
            enddo
            do jm=nmode+1,nmm
               mfun(jzmf,jm)=NaN
               mfunph(jzmf,jm)=NaN
            enddo
         enddo
cpg         if(iimpx .eq. 1 .or. iimpx .eq. 3) then
cpg            call hdf_write_gen(3+ir,outroot,lout,zmf,nzmf,xmode,nmm,
cpg     .         fcw,nfcw,var_ax,nrunx,1.,1,mfun,1,nzmf,1,nmm,
cpg     .         jfcw,jfcw,jrun,jrun,1,1,'Depth - m',9,'Mode',4,
cpg     .         'Frequency - Hz',14,'Parameter',9,' ',1,'mag(phi)',
cpg     .         8,jset0+2,1,2)
cpg         endif
cpg         if(iimpx .eq. 2 .or. iimpx .eq. 3) then
cpg            call hdf_write_gen(3+ir,outroot,lout,zmf,nzmf,xmode,nmm,
cpgcpg     .         fcw,nfcw,var_ax,nrunx,1.,1,mfunph,1,nzmf,1,nmm,
cpg     .         jfcw,jfcw,jrun,jrun,1,1,'Depth - m',9,'Mode',4,
cpg     .         'Frequency - Hz',14,'Parameter',9,' ',1,'arg(phi)',
cpg     .         8,jset0+3,1,2)
cpg         endif
      endif
c
      return
      end
