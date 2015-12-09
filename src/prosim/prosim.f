         program prosim
c
c: Normal mode model for acoustic propagation in range dependent ocean
c: environments with fluid and/or solid layered structures above and below.  
c
c:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************
c: References for the code as of 3/22/95 are:
c:    E. K. Westwood, C. T. Tindle, and N. R. Chapman, "A normal mode
c:       model for multilayered acoustoelastic ocean environments based on
c:       an analytic reflection coefficient method," J. Acoust. Soc. Am.,
c:       95, No. 5, Pt. 2, 2908 (1994).
c:    E. K. Westwood, "An efficient broadband normal-mode model for
c:       acoustoelastic ocean environments," J. Acoust. Soc. Am., 96,
c:       No. 5, Pt. 2, 3352 (1994).
c
c: *******************************
c: *     REVISION (1996):        *
c: *         E M G  GROUP        *
c: *     S A C L A N T C E N     *
c: *******************************

      implicit none
      include 'Parms_com'
      include 'i_o_svp_com'
      include 'i_o_opt_com'
      include 'i_o_1b_com'
      include 'i_o_2_com'
      include 'deltaf_com'
      include 'depth_com'

                      
      integer*4 isectr,jflo,jfhi,intrp(RGMAX)

      lout= 10
      outroot(1:lout)= 'prosim_out'
      Outfile=outroot(1:lout)//'.dat'
      loutf=lout + 4
c      open(2,file=outfile(1:loutf),form='formatted',
c     &     status= 'UNKNOWN')
      print *,'Running PROSIM ...'
c      call flush(6)
c      ncall=0

c: Read input file:
        PRINT *, ' prosim, calling svp_read '
        call svp_read
        close(10)

c: Initialize flags
        isubb = 1
        ishft = 1
        ncountr = 0
        jfhi_up=2
c: Set iiwrite=1 so that output files and messages are sent:
        iiwrite=1
        iifft =0
        iifail =0
        nmode = 0

        do while(jfhi_up .gt. 1) !START SUBBAND CALCULATIONS

          print *, '# OF SUBBAND:', isubb

c: Recall & check input options:
          call opt_read
          if(isubb .eq. 1) then
            call svp_adj(intrp)
          endif

c: Check input parameters:
          isectr=1
          rstart=0.0
          rend=0.0

c   Initialize nmodemin until # of freq. within subbands
          nmodemin=100000

          do isectr=1,nsctor     !START SECTOR CALCULATIONS
            call tenv_read_tilt(isectr,intrp,isubb)
            call svp_check
            call svp_check2

            print *, '# OF SECTOR: ', isectr
            print *, 'LENGHT OF SECT: ', secleng(isectr)
            call rx_bb(isectr,jfhi,jflo,intrp)

            rstart = rend
            rend=secleng(isectr)+rstart

          enddo                  !FINISH SECTOR CALCULATIONS

          isubb=isubb+1
          jfhi_up=jflo

        enddo                    !FINISH SUBBAND CALCULATIONS

c      if(iifft .ne. 0) then
        call bb_fft_out()
c        call bart_cor()
c      endif
c
      close(77)
      close(78)
      close(88)
      close(59)
      close(51)
      print *,'PROSIM finished. Listing on file ',outfile(1:loutf)

      stop
      end
