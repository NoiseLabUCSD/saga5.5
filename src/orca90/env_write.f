      subroutine env_write(env_mat)
c
c: Writes out file containing the environmental parameters vs. depth.
c
      use parms_com
      use i_o_com
      use gen_com
c
c: Interpolates geoacoustic profile at n_env points, fills
c: env_mat(n_env,0:5) with profiles of cp,cs,rho,kp,ks, and writes
c: out the profile to an ASCII file if desired.
c
      integer*4 j,jlay,nvalx,jj,jj1,jpt,j_env
      real*8 delz,facz,env_mat(n_env+20,0:5)
c
      delz=abs(zdep(nlay-1)-zdep(1))/max(1,n_env-1)
      open(20,file=outroot(1:lout)//'_prof',form='formatted')
      j_env=0
cc    print *,'zdep = ',nlay,(zdep(j),j=1,nlay)
      do jlay=2,nlay
         if(jlay .eq. jlfake(1) .or. jlay .eq. jlfake(2)) goto 99
         nvalx=1
         jpt=j_env + 1
         if(jlay .eq. 1) then
            j_env=j_env + 1
            env_mat(j_env,0)=zdep(1)
         elseif(jlay .eq. nlay) then
            j_env=j_env + 1
            env_mat(j_env,0)=zdep(nlay-1)
         else
            nvalx=floor(h(jlay)/delz)
            facz=h(jlay)/max(1,nvalx)
            jj1=1
c: If mismatch above, include first point at top of layer:
            if(mm(jlay-1) .eq. 1) jj1=jj1 - 1
            do jj=jj1,nvalx
               j_env=j_env + 1
               env_mat(j_env,0)=zdep(jlay-1) + jj*facz
            enddo
         endif
cc    print *,'jlay,j_env = ',jlay,j_env,nlay,n_env
c: Interpolate cp and kp:
         call c_interp(zdep(jlay-1),zdep(jlay),geo(1,1,jlay),
     .      geo(2,1,jlay),geo(1,4,jlay),geo(2,4,jlay),geo(1,3,jlay),
     .      geo(2,3,jlay),env_mat(jpt,0),env_mat(jpt,1),env_mat(jpt,4),
     .      env_mat(jpt,3),j_env-jpt+1)
c: Interpolate cs and ks:
         call c_interp(zdep(jlay-1),zdep(jlay),geo(1,2,jlay),
     .      geo(2,2,jlay),geo(1,5,jlay),geo(2,5,jlay),geo(1,3,jlay),
     .      geo(2,3,jlay),env_mat(jpt,0),env_mat(jpt,2),env_mat(jpt,5),
     .      env_mat(jpt,3),j_env-jpt+1)
         do jj=jpt,j_env
            write(20,100) (env_mat(jj,j),j=0,5)
         enddo
99       continue
      enddo
100   format(f10.3,1x,2(f10.4,1x),g10.4,1x,2(f10.6,1x))
c
      return
      end
ccc
      subroutine c_interp(z1,z2,c1,c2,alp1,alp2,rho1,rho2,z,c,alp,rho,
     .   nvalx)
c
c: Interpolates 1/c**2 linear profile, including attenuation.
c
      implicit none
      integer*4 nvalx,j
      real*8 z1,z2,c1,c2,alp1,alp2,rho1,rho2,fac,drho_dz
      real*8 z(nvalx),c(nvalx),alp(nvalx),rho(nvalx)
      complex*16 xk1sq,xk2sq,xk,dkdz
c: fac=2*pi*1000*20*log(e):
      data fac/5.457505415367364d+04/
c
      if(c1 .eq. 0.d0 .or. c2 .eq. 0.d0) then
         do j=1,nvalx
            c(j)=0.e0
            alp(j)=0.e0
         enddo
         return
      endif
c
      xk1sq=dcmplx(1.d0/c1,alp1/fac)**2
      xk2sq=dcmplx(1.d0/c2,alp2/fac)**2
      dkdz=(xk2sq - xk1sq)/max(1.d-100,z2 - z1)
      drho_dz=(rho2 - rho1)/max(1.d-100,z2 - z1)
      do j=1,nvalx
         xk=cdsqrt(xk1sq + (z(j)-z1)*dkdz)
         c(j)=1.d0/dreal(xk)
         alp(j)=fac*dimag(xk)
         rho(j)=rho1 + (z(j)-z1)*drho_dz
      enddo
c
      return
      end
