      subroutine vbl_init
c
c: Initializes variables in scairy_com.
c
      use scairy_com
      use gen_com
      character*8 chnan
c
      ei13=( 0.5D0, 0.8660254037844387D0)
      ei23=(-0.5D0, 0.8660254037844387D0)
      ei16=( 0.8660254037844387D0, 0.5D0)
      ei56=(-0.8660254037844387D0, 0.5D0)
      eim13=( 0.5D0, -0.8660254037844387D0)
      eim23=(-0.5D0, -0.8660254037844387D0)
      eim16=( 0.8660254037844387D0, -0.5D0)
      eim56=(-0.8660254037844387D0, -0.5D0)
      eye=(0.D0,1.D0)
cc    zpt=(-0.9D0,2.8D0)
      pie23=2.094395102393195D0
      pie_x=3.1415926535897932D0
      pie_inv=0.31830988618379D0
      det_bi=(0.31830988618379,0.)
      det_pos=(0.13783222385545,0.07957747154595)
      det_neg=(0.13783222385545,-0.07957747154595)
      sqrt3=1.73205080756888
c
c: Create real*4 with NaN in it to initialize non-full HDF arrays:
      chnan='7FC00000'
      read(chnan,'(z8)') NaN

      return
      end
