c:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************
c: *******************************
c: *     REVISION (1996):        *
c: *         E M G  GROUP        *
c: *     S A C L A N T C E N     *
c: *******************************                                          
      subroutine bb_init(isum)
c
c: Checks and initializes variables for broadband mode calculations.
c
      implicit none
      include 'Parms_com'
      include 'i_o_opt_com'
      common /file_com/ svp_file,opt_file,outroot,outfile
      include 'i_o_2_com'
      include 'gen3_com'
c
      integer*4 n,j,iibad,isum
      character*64 svp_file,opt_file,outroot,outfile
c
c: Initialize source and receiver geometry arrays:
      call sr_geom(isum)
c
      if(isum .eq. 2) then
         if(.not.multcw) then
         if(iiwrite.ge.2) then
            write(6,*)'Delta frequency= ',df
            write(6,*)'Maximum time window= ',1./df
c     fmc
            print *,' from subroutine bb_init.f :'
            print *,' fs,   nfft, Tw ',  fsbb, nfftbb, Tw
            print *,' fmin, fmax, df :', fmin, fmax, df
            print *,' nf1,  nf2,  nfreq :', nf1, nf2, nfbb
c     fmc
            iibad=0
            call mem_lim(nfbb,NFBBMAX,MLINE,LML,
     .           'nfbb',4,'NFBBMAX',7,iibad,0)
            call mem_lim(nfbb*nrec*nsrc,NTFMAX,
     .           MLINE,LML,'nfbb*nrec*nsrc',14,
     .           'NTFMAX',6,iibad,0)
            if(iibad .eq. 1) stop
c     
         endif
         endif
      endif
      return
      end
