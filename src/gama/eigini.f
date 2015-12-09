      subroutine eigini(kr,kl,xl,nray,jrec)
c
c: this subroutine initializes the eigenray list header and the
c: bar graph plot.
c
      implicit integer*4(i-n)
      include 'common/depth'
      include 'common/srlox'
      include 'common/svp'
      include 'common/gamaoptions'
      include 'common/vhfacs' 
      include 'common/pii'
      include 'common/caustix'
      include 'common/charcom'
      character*32 tlab,dblab,elab,silab,labb
      integer*4 kl(2,2500)
      real xl(7,2500)
c
      data tlab/'arrival time (sec)'/,labb/' '/
      data elab/'f-dep attn (db/khz) '/
      data silab/'receiver grazing angle (deg) (x)'/
      data dblab/'total attn(db) (o)'/
      data dimh/10.75/,dimv/8.0/
c
      emin=1.e13
      emax=-1.e13
      tmax=0.
      simin=1.e13
      jc=0
      do 10 nr=1,nray
         if(jc .eq. 1) then
            jc=0
            goto 10
         endif
         emin=min(emin,xl(4,nr))
         emax=max(emax,xl(4,nr))
         tmax=max(tmax,xl(3,nr))
         simin=min(simin,xl(1,nr))
         jc=mod(kl(1,nr)/100,10)
10    continue
c
      if(iieig .eq. 1) then
         write(59,102) zs,zr,.001*rng,rsrc,jzs,jrec
         write(59,103) kr,cs,cr,freq
         write(59,110)
         write(59,120)
         write(59,130)
         if(nray .eq. 0) then 
            write(59,140)
140   format('NO EIGENRAYS FOUND FOR THIS SRC/REC CONFIGURATION')
            dtau=0
         else
            dtau=tmax - tminsr
         endif
      endif
102   format(//8x,'ZS = ',f7.2,' M; ZR = ',f7.2,' M; RANGE = ',f9.2,
     .   ' KM; T/S = ',f9.2,'; JZS = ',i2,'; JREC = ',i2)
103   format(10x,'SRC SEQ # = ',i4,'; CS = ',f8.3,' M/S; CR = ',
     .   f8.3,' M/S; FREQUENCY = ',f10.2,' Hz')
110   format(/' RAY OCEAN PATH  BOTLAYTRAV   SRC ',
     .   '   REC   TURN    TRAVEL  FDEP REFL/TRAN GEOM ',
     .   'SURF SUBS   TOTAL  TOT  #')
120   format( '  #   ID #T #B    1 2 3 4 5  ANGLE',
     .   '  ANGLE   SV    TIME (S) ATTN  MAG   PH LOSS ',
     .   'SCAT SCAT   ATTEN   PH PS')
130   format( '---- ----------- ---------- ------',
     .   ' ------ ------ --------- ---- ---- ---- ---- ',
     .   '---- ---- ------- ---- --')
c
      if(iibar .eq. 1) then
         simin=acos(simin*cr)*piedeg
         dbmin=20.*log10(cutmag*dbsr)
         dbmx=20.*log10(dbsr)
         call axlabel(tminsr,tmax,tlo,thi,tdelt)
         call axlabel(dbmin,dbmx,dblo,dbhi,dbdel) 
         call axlabel(emin,emax,elo,ehi,edelt)
         call axlabel(-1.*simin,simin,silo,sihi,sdel)
         if(dbmin .eq. dblo) dblo=dblo-dbdel
c: draw axes on bar graph for the current range.
cxx      call pltlfn(l"raypic")
cxx      call pltdim(11.75,9.0,-1)
cxx      call pltorg(.5,.5)
cxx      call pltaxis(0.,0.,dimh,0.0,tlo,thi,tdelt/10.,tlab,-18,10)
cxx      call pltaxis(0.,dimv/2.,dimh,0.,tlo,thi,tdelt/10.,labb,1,0)
cxx      call plt(0.,dimv,3)
cxx      call plt(dimh,dimv,2)
cxx      call pltaxis(0.,0.,dimv/2.,90.,0.,ehi,edelt/2.,elab,23,4)
cxx      call pltaxis(dimh,0.,dimv,90.,silo,sihi,sdel/2.,silab,-32,2) 
cxx      call pltaxis(0.,dimv/2.,dimv/2.,90.0,dblo,dbhi,dbdel/2.,
cxx  .      dblab,18,4)
         hfac=dimh/(thi-tlo)
         vhfac=dimv/2./(dbhi-dblo)
         vlfac=dimv/2./ehi
         vlfac2=dimv/(sihi-silo)
c
cxx      call pltline(dimh/2.,dimv+.115,-.16)
cxx      write(8,200) ktitle(1:lktit),zs,zr,rng,rsrc,jxy,freq
cxx200      format(a,'; zs =',f7.2,'; zr =',f7.2,'; r =',f9.2,
cxc  .      '; t/s =',f9.2,'; jxy=',i2,'; f =',f8.2)
      endif
c
      return
      end 
