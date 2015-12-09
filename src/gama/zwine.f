      function zwine(sg,zz,g0,g1,g2,phi1,f1,f4)
c
c: this function returns the depth at which the ray characterized by
c: phi1 (a function of the snell invariant) turns.
c
c: due to roundoff errors, radical can go negative
      implicit integer*4(i-n)
      rad=sqrt(max(0.,(g0/2.)**2 - f4*phi1**2)) 
      bb=-1.*(g2*phi1**2 + g0/2.)
      sg=1.
      zwine=(bb + rad)/f1
      if((zwine .lt. 0.) .or. (zwine .gt. zz)) then
         sg=-1.
         zwine=(bb - rad)/f1
         if((zwine .lt. 0.) .or. (zwine .gt. zz)) then
            print *,'illegal zwine calculated: ',g0,g1,g2,zz,zwine,
     .      (bb + rad)/f1
         endif
      endif
c
      return
      end 
