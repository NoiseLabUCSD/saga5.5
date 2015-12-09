      Subroutine forwardmodel(iopt,mopt)
      integer   mopt,i,iopt(mopt)
      DO i=1,mopt
         iopt(i)=0.
      ENDDO
      iopt(1)=2
      iopt(30)=3 
      end
      subroutine map2forw
      end
