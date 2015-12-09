      subroutine eig_list(jrun,nrunx,jf)
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jrun,nrunx,jf,jm
      real*8 cphase,db_km,vg,zt_bot,zt_top,s_km
cnn   real*8 null_sp
c
      if(jrun .eq. 1 .and. jf .eq. 1) then
         open(15,file=outroot(1:lout)//'_eig',form='formatted')
      endif
      write(15,100) f_hz
100   format('Frequency = ',f9.2/
     .'MODE   V_PH(m/s)  V_G(m/s) TIME(s/km)     dB/km  ZT_TOP  ',
     .   'ZT_BOT')
cnn  .   'ZT_BOT NULL SP')
c
      do jm=1,nmode
         cphase=w/real(kn(jm))
         db_km=8685.9d0*dimag(kn(jm))
         vg=dreal(eig_char(4,jm))
         call z_turn(cphase,1,zt_bot)
         call z_turn(cphase,2,zt_top)
         zt_top=max(-999.99d0,zt_top)
         zt_bot=min(9999.99d0,zt_bot)
         s_km=1000.d0/vg
         if(cphase .lt. 1.d8 .and. vg .lt. 1.d7 .and. vg .gt. -1.d6 
     .      .and. s_km .lt. 1.d3 .and. s_km .gt. -1.d2 .and. 
     .      db_km .lt. 1.d5) then
            write(15,110) jm,cphase,vg,s_km,db_km,zt_top,zt_bot
110         format(i4,1x,f11.2,1x,f9.3,1x,f10.6,1x,f9.3,2(1x,f7.2))
         else
            write(15,112) jm,cphase,vg,s_km,db_km,zt_top,zt_bot
112         format(i4,1x,e11.5,1x,e9.3,1x,e10.3,1x,e9.3,2(1x,f7.2))
         endif
cnn      if(jm .eq. 1) then
cnn         write(15,110) jm,cphase,vg,s_km,db_km,zt_top,zt_bot
cnn      else
cnn         null_sp=dreal(.001*2.*pie/(kn(jm-1) - kn(jm)))
cnn         write(15,120) jm,cphase,vg,1000.d0/vg,db_km,
cnn  .         max(-999.99d0,zt_top),min(9999.99d0,zt_bot),
cnn  .         min(999.999d0,null_sp)
cnn120         format(i4,1x,f11.2,1x,f9.3,1x,f10.7,1x,g7.2,2(1x,f7.2),
cnn  .         1x,f7.3)
cnn      endif
      enddo
      if(jrun .eq. nrunx .and. jf .eq. nfcw) then
         close(15)
      endif
c
      return
      end
ccc
      subroutine z_turn(cphase,ii,zturn)
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 ii,ii2,j,j1,inc
      real*8 cphase,zturn
c
      nsvmin=jduct(1,kduct)
      isvmin=jduct(2,kduct)
      ii2=3 - ii
      inc=jflu(ii,3)
      j1=nsvmin
      if(ii .ne. isvmin) j1=j1+inc
      do j=j1,jsol(ii,2),inc
         if(cphase .le. geo(ii,1,j)) then
            zturn=zdep(j-1)
            return
         elseif(cphase .le. geo(ii2,1,j)) then
            if(ii .eq. 1) then
               call z_interp(zdep(j-1),zdep(j),geo(1,1,j),geo(2,1,j),
     .            geo(1,4,j),geo(2,4,j),isp(j),cphase,zturn)
            else
               call z_interp(zdep(j),zdep(j-1),geo(2,1,j),geo(1,1,j),
     .            geo(2,4,j),geo(1,4,j),isp(j),cphase,zturn)
            endif
            return
         endif
      enddo
      if(ii .eq. 1) then
         zturn=zdep(nlay-1)
      else
         zturn=zdep(1)
      endif
c
      return
      end
ccc
      subroutine z_interp(z1,z2,c1,c2,alp1,alp2,isp,c,z)
c
c: Interpolates 1/c**2 linear profile, including attenuation.
c
      implicit none
      integer*4 isp
      real*8 z1,z2,c1,c2,alp1,alp2,c,z,alp,fac
      complex*16 xk1sq,xk2sq,xksq
c: fac=2*pi*1000*20*log(e):
      data fac/5.457505415367364d+04/
c
c: Assume attenuation is approximately linear:
      if(isp .eq. 1) then
         z=z1
         return
      endif
      alp=alp1 + (c-c1)*(alp2-alp1)/(c2-c1)
      xk1sq=dcmplx(1.d0/c1,alp1/fac)**2
      xk2sq=dcmplx(1.d0/c2,alp2/fac)**2
      xksq=dcmplx(1.d0/c,alp/fac)**2
      z=z1 + dreal((xksq-xk1sq)*(z2-z1)/(xk2sq-xk1sq))
c
      return
      end
