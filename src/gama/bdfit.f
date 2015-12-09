      subroutine bdfit(kbdf1,kbdf2)
c
c: this subroutine computes the interpolation coefficients for the
c: eigenray characteristics that vary with frequency as a result 
c: of having beam displacement.
c
      implicit integer*4(i-n)
      include 'common/bdcom'
      include 'common/bdchar'
      include 'common/pii'
c
c: this bd eigenray was found for freqs bdf(kbdf1) through bdf(kbdf2):
      kk=kbdf2+1
      do 10 k=kbdf1,kk
         nparm(k)=6
c: make eig parameters constant for f<bdf(kbdf1) and f>bdf(kbdf2) by
c: setting j=k except for k=kbdf2+1, where j=kbdf2:
         j=k
         if(k .eq. kk) j=kbdf2
c: setting k=k-1 except for k=kbdf1, where km1=k
         km1=k-1
         if(k .eq. kbdf1) km1=k
         if((jc(km1) .eq. 0) .and. (jc(j) .eq. 0)) then
c: normal to normal:
            call lnfit(k,tbd(km1),tbd(j),parma(k,1),parmb(k,1))
            call lnfit(k,trpbd(km1),trpbd(j),parma(k,2),parmb(k,2))
            gtf1=gbd(km1)*trmbd(km1)*expo(bdf(km1),ebd(km1))
            gtf2=gbd(j)*trmbd(j)*expo(bdf(j),ebd(j))
            call lnfit(k,gtf1,gtf2,parma(k,3),parmb(k,3))
            nparm(k)=3
         elseif((jc(km1) .eq. 1) .and. (jc(j) .eq. 1)) then
c: shadow to shadow:
            call lnfit(k,tbd(km1),tbd(j),parma(k,1),parmb(k,1))
            call lnfit(k,trpbd(km1),trpbd(j),parma(k,2),parmb(k,2))
            gtf1=expo(bdf(km1),ebd(km1))
            gtf2=expo(bdf(j),ebd(j))
            call lnfit(k,gtf1,gtf2,parma(k,3),parmb(k,3))
            xairy1=bdom23(km1)*xaibd(km1)
            xairy2=bdom23(j)*xaibd(j)
            call lnfit(k,xairy1,xairy2,parma(k,4),parmb(k,4))
            s11=s1fbd(km1)*bdom16(km1)
            s12=s1fbd(j)*bdom16(j)
            call lnfit(k,s11,s12,parma(k,5),parmb(k,5))
            s21=s2fbd(km1)/bdom16(km1)
            s22=s2fbd(j)/bdom16(j)
            call lnfit(k,s21,s22,parma(k,6),parmb(k,6))
         elseif((jc(km1) .eq. 2) .and. (jc(j) .eq. 2)) then
c: doublet to doublet:
            tav1=.5*(tbd(km1)+t2bd(km1))
            tav2=.5*(tbd(j)+t2bd(j))
            call lnfit(k,tav1,tav2,parma(k,1),parmb(k,1))
            phav1=.5*(trpbd(km1)+trp2bd(km1))+piequ
            phav2=.5*(trpbd(j)+trp2bd(j))+piequ
            call lnfit(k,phav1,phav2,parma(k,2),parmb(k,2))
            gtf1=expo(bdf(km1),.5*(ebd(km1)+e2bd(km1)))
            gtf2=expo(bdf(j),.5*(ebd(j)+e2bd(j)))
            call lnfit(k,gtf1,gtf2,parma(k,3),parmb(k,3))
            sg=sign(1.,t2bd(j) - tbd(j))
            xairy1=(.75*abs(bdom(km1)*(t2bd(km1)-tbd(km1)) + 
     .         (trp2bd(km1)-trpbd(km1))))**(.66667) 
            xairy2=(.75*abs(bdom(j)*(t2bd(j)-tbd(j)) + (trp2bd(j)-
     .         trpbd(j))))**(.66667) 
            xaip1=xairy1**(.25)
            xaip2=xairy2**(.25)
            call lnfit(k,xairy1,xairy2,parma(k,4),parmb(k,4))
            gstr12=g2bd(km1)*trm2bd(km1)
            gstr1=gbd(km1)*trmbd(km1)
            gstr22=g2bd(j)*trm2bd(j)
            gstr2=gbd(j)*trmbd(j)
            s11=sqpie*(gstr12 + gstr1)*xaip1
            s21=sg*sqpie*(gstr12 - gstr1)/xaip1
            s12=sqpie*(gstr22 + gstr2)*xaip2
            s22=sg*sqpie*(gstr22 - gstr2)/xaip2
            call lnfit(k,s11,s12,parma(k,5),parmb(k,5))
            call lnfit(k,s21,s22,parma(k,6),parmb(k,6))
         elseif((jc(km1) .eq. 2) .and. (jc(j) .eq. 0)) then
c: doublet changing to 1 normal ray (lateral wave probably too weak):
            nparm(k)=3
c: find stronger of two doublet rays and fit with single normal ray:
            if(g2bd(km1) .gt. gbd(km1)) then
               call lnfit(k,t2bd(km1),tbd(j),parma(k,1),parmb(k,1))
               call lnfit(k,trp2bd(km1),trpbd(j),parma(k,2),parmb(k,2))
               gtf1=g2bd(km1)*trm2bd(km1)*expo(bdf(km1),e2bd(km1))
               gtf2=gbd(j)*trmbd(j)*expo(bdf(j),ebd(j))
               call lnfit(k,gtf1,gtf2,parma(k,3),parmb(k,3))
            else
               call lnfit(k,tbd(km1),tbd(j),parma(k,1),parmb(k,1))
               call lnfit(k,trpbd(km1),trpbd(j),parma(k,2),parmb(k,2))
               gtf1=gbd(km1)*trmbd(km1)*expo(bdf(km1),ebd(km1))
               gtf2=gbd(j)*trmbd(j)*expo(bdf(j),ebd(j))
               call lnfit(k,gtf1,gtf2,parma(k,3),parmb(k,3))
            endif
         elseif((jc(km1) .eq. 0) .and. (jc(j) .eq. 2)) then
c: 1 normal ray changing to doublet (lateral wave probably too weak):
            nparm(k)=3
c: find stronger of two doublet rays and fit with single normal ray:
            if(g2bd(j) .gt. gbd(j)) then
               call lnfit(k,tbd(km1),t2bd(j),parma(k,1),parmb(k,1))
               call lnfit(k,trpbd(km1),trp2bd(j),parma(k,2),parmb(k,2))
               gtf1=gbd(km1)*trmbd(km1)*expo(bdf(km1),ebd(km1))
               gtf2=g2bd(j)*trm2bd(j)*expo(bdf(j),e2bd(j))
               call lnfit(k,gtf1,gtf2,parma(k,3),parmb(k,3))
            else
               call lnfit(k,tbd(km1),tbd(j),parma(k,1),parmb(k,1))
               call lnfit(k,trpbd(km1),trpbd(j),parma(k,2),parmb(k,2))
               gtf1=gbd(km1)*trmbd(km1)*expo(bdf(km1),ebd(km1))
               gtf2=gbd(j)*trmbd(j)*expo(bdf(j),ebd(j))
               call lnfit(k,gtf1,gtf2,parma(k,3),parmb(k,3))
            endif
         elseif((jc(km1) .eq. 1) .and. (jc(j) .eq. 2)) then
c: shadow changing to doublet:
            tav2=.5*(tbd(j)+t2bd(j))
            call lnfit(k,tbd(km1),tav2,parma(k,1),parmb(k,1))
            phav2=.5*(trpbd(j)+trp2bd(j))+piequ
            call lnfit(k,trpbd(km1),phav2,parma(k,2),parmb(k,2))
            gtf1=expo(bdf(km1),ebd(km1))
            gtf2=expo(bdf(j),.5*(ebd(j)+e2bd(j)))
            call lnfit(k,gtf1,gtf2,parma(k,3),parmb(k,3))
            xairy1=bdom23(km1)*xaibd(km1)
            xairy2=(.75*abs(bdom(j)*(t2bd(j)-tbd(j)) + (trp2bd(j)-
     .         trpbd(j))))**(.66667) 
            call lnfit(k,xairy1,xairy2,parma(k,4),parmb(k,4))
c
            sg=sign(1.,t2bd(j) - tbd(j))
            xaip2=xairy2**(.25)
            gstr22=g2bd(j)*trm2bd(j)
            gstr2=gbd(j)*trmbd(j)
            s11=s1fbd(km1)*bdom16(km1)
            s21=s2fbd(km1)/bdom16(km1)
            s12=sqpie*(gstr22 + gstr2)*xaip2
            s22=sg*sqpie*(gstr22 - gstr2)/xaip2
            call lnfit(k,s11,s12,parma(k,5),parmb(k,5))
            call lnfit(k,s21,s22,parma(k,6),parmb(k,6))
         elseif((jc(km1) .eq. 2) .and. (jc(j) .eq. 1)) then
c: doublet changing to shadow:
            tav1=.5*(tbd(km1)+t2bd(km1))
            call lnfit(k,tav1,tbd(j),parma(k,1),parmb(k,1))
            phav1=.5*(trpbd(km1)+trp2bd(km1))+piequ
            call lnfit(k,phav1,trpbd(j),parma(k,2),parmb(k,2))
            gtf1=expo(bdf(km1),.5*(ebd(km1)+e2bd(km1)))
            gtf2=expo(bdf(j),ebd(j))
            call lnfit(k,gtf1,gtf2,parma(k,3),parmb(k,3))
            xairy1=(.75*abs(bdom(km1)*(t2bd(km1)-tbd(km1))+(trp2bd(km1)-
     .         trpbd(km1))))**(.66667) 
            xairy2=bdom23(j)*xaibd(j)
            call lnfit(k,xairy1,xairy2,parma(k,4),parmb(k,4))
c
            sg=sign(1.,t2bd(km1) - tbd(km1))
            xaip1=xairy1**(.25)
            gstr12=g2bd(km1)*trm2bd(km1)
            gstr1=gbd(km1)*trmbd(km1)
            s11=sqpie*(gstr12 + gstr1)*xaip1
            s21=sg*sqpie*(gstr12 - gstr1)/xaip1
            s12=s1fbd(j)*bdom16(j)
            s22=s2fbd(j)/bdom16(j)
            call lnfit(k,s11,s12,parma(k,5),parmb(k,5))
            call lnfit(k,s21,s22,parma(k,6),parmb(k,6))
         elseif((jc(km1) .eq. 1) .and. (jc(j) .eq. 0)) then
c: shadow changing to normal:
            call lnfit(k,tbd(km1),tbd(j),parma(k,1),parmb(k,1))
            call airy(-bdom23(km1)*xaibd(km1),ai,aip)
            ccre=ai*s1fbd(km1)*bdom16(km1)
            ccim=aip*s2fbd(km1)/bdom16(km1)
            ph1=trpbd(km1) + atan2(ccim,ccre)
            call lnfit(k,ph1,trpbd(j),parma(k,2),parmb(k,2))
            gtf1=expo(bdf(km1),ebd(km1))*sqrt(ccre**2+ccim**2)
            gtf2=gbd(j)*trmbd(j)*expo(bdf(j),ebd(j))
            call lnfit(k,gtf1,gtf2,parma(k,3),parmb(k,3))
            nparm(k)=3
         elseif((jc(km1) .eq. 0) .and. (jc(j) .eq. 1)) then
c: normal changing to shadow:
            call lnfit(k,tbd(km1),tbd(j),parma(k,1),parmb(k,1))
            call airy(-bdom23(j)*xaibd(j),ai,aip)
            ccre=ai*s1fbd(j)*bdom16(j)
            ccim=aip*s2fbd(j)/bdom16(j)
            ph2=trpbd(j) + atan2(ccim,ccre)
            call lnfit(k,trpbd(km1),ph2,parma(k,2),parmb(k,2))
            gtf1=gbd(km1)*trmbd(km1)*expo(bdf(km1),ebd(km1))
            gtf2=expo(bdf(j),ebd(j))*sqrt(ccre**2+ccim**2)
            call lnfit(k,gtf1,gtf2,parma(k,3),parmb(k,3))
            nparm(k)=3
         else
            print *,'unexpected jc(km1),jc(j) = ',jc(km1),jc(j)
         endif
10    continue
c
c: make multiplicative factor go to zero below kbdf1 and above kbdf2:
      if(kbdf1 .gt. 1) then
         kk=kbdf1+1
         gtf2=parma(kk,3)*bdf(kk) + parmb(kk,3)
         call lnfit(kbdf1,0.,gtf2,parma(kbdf1,3),parmb(kbdf1,3))
      endif
      if(kbdf2 .lt. nbdf) then
         kk=kbdf2+1
         gtf1=parma(kbdf2,3)*bdf(kbdf2) + parmb(kbdf2,3)
         call lnfit(kk,gtf1,0.,parma(kk,3),parmb(kk,3))
      endif
c
      return
      end
