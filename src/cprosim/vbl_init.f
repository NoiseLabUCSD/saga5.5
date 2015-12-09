      subroutine vbl_init
c
c: Initializes variables in scairy_com and lab_com.
c
      implicit none
      include 'scairy_com'
      include 'lab_com'
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
cc
      flab='Frequency  - Hz'
      rlab='Range  - km'
      mnlab='Mode Number'
      tlab='Time  - s'
      krlab='RE[k]'
      kilab='IM[k]'
      thlab='Grazing Angle  - deg'
      mlab='|R| - dB'
      phlab='Phase of R  - deg'
      dlab='Depth  - m'
      pllab='Propagation Loss - dB'
      mrlab='RE[Mode Amplitude]'
      milab='IM[Mode Amplitude]'
      malab='Mode Amplitude'
      mplab='Mode Phase'
      dblab='dB'
      mtlab='m'
      kmlab='km'
      z4=0. 
c
      return
      end
