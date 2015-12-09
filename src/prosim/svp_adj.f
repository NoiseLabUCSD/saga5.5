c: *******************************
c: *   AUTHOR (1996):            *
c: *         E M G  GROUP        *
c: *     S A C L A N T C E N     *
c: *******************************
      subroutine svp_adj(intrpl)
c
      implicit none
      include 'Parms_com'
      include 'i_o_opt_com'
      common /ny_com/ narray(NVRMAX),act_sect
      integer*4 narray,act_sect
      common /out_com3/ nzsr, nlay,NSCTOR,dmmy1
      integer*4 nzsr, nlay,NSCTOR,dmmy1
      include 'depth_com'
c: Local variables:
      integer*4 displc,r,intrpl(RGMAX),nsrec,ISECT
      real*8 sec_int, sect_nd, rnge
c
      sec_int = 0.0
      sect_nd = 0.0
      displc = 0.0
      nsrec = nsrc
c      nsrec = -nsrc
c     EVALUATE THE NUMBER OF ARRAY IN EACH SECTOR
      intrpl(NSCTOR) = 0
      narray(NSCTOR) = 0
      DO ISECT = 1, NSCTOR-1
c       sect_nd = secleng(ISECT)
       sect_nd = secleng(ISECT) + sect_nd
       intrpl(ISECT) = 0
       narray(ISECT) = 0
c       DO r = 1+displc, nsrec
       DO r = 1, nsrec
         rnge = rkm(r)
         if((sec_int.lt.rnge).and.(sect_nd.ge.rnge))then 
            intrpl(ISECT) = intrpl(ISECT) + 1
            narray(ISECT) = narray(ISECT) + 1
         endif
       END DO
c       displc = narray(ISECT)
       sec_int = sect_nd
      END DO 
c
      if(NSCTOR .gt. 1) displc = narray(ISECT-1)
      DO r = 1+displc, nsrec
        rnge = rkm(r)
        if(sec_int.lt.rnge)then 
           intrpl(ISECT) = intrpl(ISECT) + 1
           narray(ISECT) = narray(ISECT) + 1
        endif
      END DO
c
c     SUM THE NUMBER OF ARRAYS
      DO ISECT = 2, NSCTOR
        intrpl(ISECT)=intrpl(ISECT)+intrpl(ISECT-1)
      END DO 
c      DO ISECT = NSCTOR,2,-1
c        narray(ISECT)=narray(ISECT)+narray(ISECT-1)
c      END DO     
c
      return
c
      end
