      SUBROUTINE powell(p,xi,xstep,xmin,n,np,mp,ftol,
     1     iter,fret,iforwloc,itmax)

c     p start point WITH n parameters
c     xi eigen vectors
c     xstep for normalization of point
c     np, mp dimesions
c     iter number of iterations
c     fret value of objective FUNCTION 
c     iforwloc

      INTEGER iter,n,np,mp,NMAX,ITMAX,itcnt,iforwloc
      REAL fret,ftol,p(np),xi(mp,mp),func
      PARAMETER (NMAX=500)
C      USES func,linmin
      INTEGER i,ibig,j,np1,ncom
      REAL del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX),xerr
      REAL xstep(np),xmin(np)
      COMMON /myiter/itcnt
      REAL pcom(NMAX),xicom(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
c      WRITE(*,*)'Entering powell',itmax
      ncom=n
      itcnt=0
      np1=np
      DO 11 j=1,n
         p(j)=(p(j)-xmin(j))/xstep(j)
         pt(j)=p(j)
c      WRITE(*,*)'p(j),xmin(j),xstep(j)',p(j),xmin(j),xstep(j),n,np,mp
 11   CONTINUE
      fp=func(p)
      if (fp.lt.10) then
         fret=fp
      endif
      WRITE(*,*)'initial fitness from powell',fret,fp

      iter=0
 1    iter=iter+1
      IF(iter.GT.ITMAX) THEN
         WRITE(*,*) 'test1, powell exceeding maximum iterations'
         GOTO 99
      ENDIF
      fp=fret
      ibig=0
      del=0.
      DO 13 i=1,n
         DO 12 j=1,n
            xit(j)=xi(j,i)
 12      CONTINUE
         fptt=fret
c       WRITE(*,*)'linesearch parameter',i,iter
         CALL linmin(p,xit,n,fret)
         IF(ABS(fptt-fret).GT.del)THEN
            del=ABS(fptt-fret)
            ibig=i
         ENDIF
c       WRITE(*,*)'  fptt,fret,del, #it', fptt,fret,del,itcnt
 13   CONTINUE
c      WRITE(*,*)'rel impr',2.*ABS(fp-fret)/(ABS(fp)+ABS(fret))
      
      xerr=0
      DO j=1,n
         xerr=MAX(xerr,ABS(p(j)-pt(j)))
      ENDDO
      WRITE(*,*)'Largest step,Improvement in costfunction',
     1     xerr,ABS(fp-fret)
      IF(2.*ABS(fp-fret).LE.ftol*(ABS(fp)+ABS(fret)) .OR. 
     1     xerr.LT.0.01 ) GOTO 99
      DO 14 j=1,n
         ptt(j)=2.*p(j)-pt(j)
         xit(j)=p(j)-pt(j)
         pt(j)=p(j)
 14   CONTINUE
      fptt=func(ptt)
      IF(fptt.GE.fp)GOTO 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      IF(t.GE.0.)GOTO 1
      WRITE(*,*)'linesearch stepest descend ',i,iter
      CALL linmin(p,xit,n,fret)
      DO 15 j=1,n
         xi(j,ibig)=xi(j,n)
         xi(j,n)=xit(j)
 15   CONTINUE
      IF(iter.EQ.ITMAX) THEN
         WRITE(*,*) 'powell exceeding maximum iterations'
         GOTO 99
      ENDIF
      GOTO 1

 99   DO j=1,n
         p(j)=P(j)*xstep(j)+xmin(j)
      ENDDO
      iforwloc=itcnt
      END
c     
c     
c     
c     
      SUBROUTINE linmin(p,xi,n,fret)
      INTEGER n,NMAX
      REAL fret,p(n),xi(n),TOL
      PARAMETER (NMAX=500,TOL=1.e-4)
C     U    USES brentrec,f1dim,mnbrak
      INTEGER j,ncom
      REAL ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX),brentrec
      COMMON /f1com/ pcom,xicom,ncom
      REAL maxgain,xcormax
      EXTERNAL f1dim
      ncom=n
      xcormax=0.
      DO 11 j=1,n
         pcom(j)=p(j)
         xicom(j)=xi(j)
         xcormax=MAX(xcormax,ABS(xicom(j)))
 11   CONTINUE

      IF (xcormax.LT.0.1) THEN
         DO j=1,n
            xicom(j)=xicom(j)*0.1/xcormax
         ENDDO
      ENDIF
      ax=0.
      xx=1.
c      WRITE(*,*)'search direction', (xicom(j),j=1,n),xcormax
      fa=fret
      CALL mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
c      WRITE(*,*)'calling brent'
      maxgain=0.
      DO j=1,n
         maxgain=MAX(maxgain,ABS(xicom(j)*(ax-bx)))
      ENDDO
c      WRITE(*,*)'gain?', maxgain
      IF (maxgain.GT.0.05) THEN
         fret=brentrec(ax,xx,bx,fa,fx,fb,f1dim,TOL,xmin)
      ELSE
         xmin=(ax+bx)/2
      ENDIF
      IF (xmin.NE.0) THEN
         DO 12 j=1,n
            xi(j)=xmin*xi(j)
            p(j)=p(j)+xi(j)
 12      CONTINUE
      ENDIF
c      write(*,*)'exit linmin'
      END
c     
c     
c     
      FUNCTION brentrec(ax,bx,cx,fa,fb,fc,f,tol,xmin)
      INTEGER ITMAX
      REAL brentrec,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS,fa,fb,fc
      EXTERNAL f
c     pg      PARAMETER (ITMAX=10,CGOLD=.3819660,ZEPS=1.0e-10)
      PARAMETER (ITMAX=5,CGOLD=.3819660,ZEPS=1.0e-10)
      INTEGER iter
      REAL a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
c     WRITE(*,*)'ax,bx,cx',ax,bx,cx
      a=MIN(ax,cx)
      b=MAX(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
c     fx=f(x)
c     WRITE(*,*)'x,bx,fx,fb',x,bx,fx,fb
c     peters additions 
      fx=fb
      w=a
      fw=fa
      v=b
      fv=fc
c     fv=fx
c     fw=fx

      DO 11 iter=1,ITMAX
         xm=0.5*(a+b)
         tol1=tol*ABS(x)+ZEPS
         tol2=2.*tol1
         tol2=2*(tol*ABS(x)+ZEPS+0.01)
         IF(ABS(x-xm).LE.(tol2-.5*(b-a))) GOTO 3
         IF(ABS(e).GT.tol1) THEN
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.*(q-r)
            IF(q.GT.0.) p=-p
            q=ABS(q)
            etemp=e
            e=d
            IF(ABS(p).GE.ABS(.5*q*etemp)
     *           .OR.p.LE.q*(a-x).OR.p.GE.q*(b-x)) 
     *           GOTO 1
            d=p/q
            u=x+d
            IF(u-a.LT.tol2 .OR. b-u.LT.tol2) d=SIGN(tol1,xm-x)
            GOTO 2
         ENDIF
 1       IF(x.GE.xm) THEN
            e=a-x
         ELSE
            e=b-x
         ENDIF
         d=CGOLD*e
 2       IF(ABS(d).GE.tol1) THEN
            u=x+d
         ELSE
            u=x+SIGN(tol1,d)
         ENDIF
         fu=f(u)
         IF(fu.LE.fx) THEN
            IF(u.GE.x) THEN
               a=x
            ELSE
               b=x
            ENDIF
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
         ELSE
            IF(u.LT.x) THEN
               a=u
            ELSE
               b=u
            ENDIF
            IF(fu.LE.fw .OR. w.EQ.x) THEN
               v=w
               fv=fw
               w=u
               fw=fu
            ELSE IF(fu.LE.fv .OR. v.EQ.x .OR. v.EQ.w) THEN
               v=u
               fv=fu
            ENDIF
         ENDIF
 11   CONTINUE
c     WRITE(*,*) 'brentrec exceed maximum iterations'
 3    xmin=x
      brentrec=fx
      END
c     
c     
c     
      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      REAL ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
      REAL dum,fu,q,r,u,ulim
c      WRITE(*,*)'mnbrak fa entry',fa
c      fa=func(ax)
c      WRITE(*,*)'mnbrak fb calle',fa,bx
      fb=func(bx)
c      WRITE(*,*)'mnbrak fb called',fb,fa
      IF(fb.GT.fa)THEN
         dum=ax
         ax=bx
         bx=dum
         dum=fb
         fb=fa
         fa=dum
      ENDIF
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
c      WRITE(*,*)'mnbrak fc called',fa,fb,fc,ax,ab,ac
 1    IF(fb.GE.fc)THEN
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*SIGN(MAX(ABS(q-r),TINY),q-r))
         ulim=bx+GLIMIT*(cx-bx)
         IF((bx-u)*(u-cx).GT.0.)THEN
            fu=func(u)
            IF(fu.LT.fc)THEN
               ax=bx
               fa=fb
               bx=u
               fb=fu
               RETURN
            ELSE IF(fu.GT.fb)THEN
               cx=u
               fc=fu
               RETURN
            ENDIF
            u=cx+GOLD*(cx-bx)
            fu=func(u)
         ELSE IF((cx-u)*(u-ulim).GT.0.)THEN
            fu=func(u)
            IF(fu.LT.fc)THEN
               bx=cx
               cx=u
               u=cx+GOLD*(cx-bx)
               fb=fc
               fc=fu
               fu=func(u)
            ENDIF
         ELSE IF((u-ulim)*(ulim-cx).GE.0.)THEN
            u=ulim
            fu=func(u)
         ELSE
            u=cx+GOLD*(cx-bx)
            fu=func(u)
         ENDIF
         ax=bx
         bx=cx
         cx=u
         fa=fb
         fb=fc
         fc=fu
         GOTO 1
      ENDIF
      RETURN
      END
c     
c     

      FUNCTION f1dim(x)
      REAL f1dim,func,x
      INTEGER j,ncom, NMAX
      PARAMETER (NMAX=500)
      REAL pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      DO 11 j=1,ncom
         xt(j)=pcom(j)+x*xicom(j)
c     WRITE(*,*),xt(j),pcom(j),x,xicom(j)
 11   CONTINUE
      f1dim=func(xt)
      END

      FUNCTION func(xt)
      INCLUDE 'comopt.h'
      REAL func
      INTEGER i,nmax
      PARAMETER (NMAX=500)
      INTEGER j,ncom,itcnt
      REAL pcom(NMAX),xicom(NMAX),xt(NMAX),xtt(Nmax)
      COMMON /f1com/ pcom,xicom,ncom
      COMMON /myiter/itcnt
      DO i=1,ncom
         xtt(i)=xt(i)*df(i)+fmin(i)
         IF ((xtt(i).LT.fmin(i)).OR.
     1        (xtt(i).GT.fmax(i))) THEN
c     IF ((xtt(i).LT.fmin(i)*xonel).OR.
c     1        (xtt(i).GT.fmax(i)*xoneu)) THEN
c     WRITE(*,*)'Powell: model rejected, parameter',i
c     WRITE(*,*)'xtt(i),fmin(i),fmax(i)'
c     WRITE(*,*)xtt(i),xt(i),fmin(i),fmax(i)
            func=100000000
            RETURN
         ENDIF
      ENDDO
      CALL setmodelreal(xtt)
c      WRITE(80,*) func,itcnt
c      WRITE(80,*) (xtt(i),i=1,ncom)
c      CALL flush(80)
      CALL forw2
       CALL cost(func) 
      itcnt=itcnt+1
c     WRITE(80,*) func,itcnt
c     WRITE(80,*) (xtt(i),i=1,ncom)
c     CALL flush(80)

      END
