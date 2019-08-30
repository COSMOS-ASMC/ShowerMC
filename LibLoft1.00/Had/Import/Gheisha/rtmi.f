*cmz :  3.14/16 13/03/89  14.48.41  by  nick van eijndhoven (cern)
*-- author :
      subroutine rtmi(x,f,fct,xli,xri,eps,iend,ier)
c
c *** nve 16-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (27-oct-1983)
c copied from r01utl.ssp.s  23.4.82
c
      external fct
c --- prepare iteration ---
      ier=0
      xl=xli
      xr=xri
      x=xl
      tol=x
      f=fct(tol)
      if(f)1,16,1
    1 fl=f
      x=xr
      tol=x
      f=fct(tol)
      if(f)2,16,2
    2 fr=f
      if(sign(1.,fl)+sign(1.,fr))25,3,25
c
c     basic assumption fl*fr less than 0 is satisfied.
c     generate tolerance for function values.
    3 i=0
      tolf=100.*eps
c
c
c     start iteration loop
    4 i=i+1
c
c     start bisection loop
      do 13 k=1,iend
      x=.5*(xl+xr)
      tol=x
      f=fct(tol)
      if(f)5,16,5
    5 if(sign(1.,f)+sign(1.,fr))7,6,7
c
c     interchange xl and xr in order to get the same sign in f and fr
    6 tol=xl
      xl=xr
      xr=tol
      tol=fl
      fl=fr
      fr=tol
    7 tol=f-fl
      a=f*tol
      a=a+a
      if(a-fr*(fr-fl))8,9,9
    8 if(i-iend)17,17,9
    9 xr=x
      fr=f
c
c     test on satisfactory accuracy in bisection loop
      tol=eps
      a=abs(xr)
      if(a-1.)11,11,10
   10 tol=tol*a
   11 if(abs(xr-xl)-tol)12,12,13
   12 if(abs(fr-fl)-tolf)14,14,13
   13 continue
c     end of bisection loop
c
c     no convergence after iend iteration steps followed by iend
c     successive steps of bisection or steadily increasing function
c     values at right bounds. error return.
      ier=1
   14 if(abs(fr)-abs(fl))16,16,15
   15 x=xl
      f=fl
   16 return
c
c     computation of iterated x-value by inverse parabolic interpolation
   17 a=fr-f
      dx=(x-xl)*fl*(1.+f*(a-tol)/(a*(fr-fl)))/tol
      xm=x
      fm=f
      x=xl-dx
      tol=x
      f=fct(tol)
      if(f)18,16,18
c
c     test on satisfactory accuracy in iteration loop
   18 tol=eps
      a=abs(x)
      if(a-1.)20,20,19
   19 tol=tol*a
   20 if(abs(dx)-tol)21,21,22
   21 if(abs(f)-tolf)16,16,22
c
c     preparation of next bisection loop
   22 if(sign(1.,f)+sign(1.,fl))24,23,24
   23 xr=x
      fr=f
      go to 4
   24 xl=x
      fl=f
      xr=xm
      fr=fm
      go to 4
c     end of iteration loop
c
c
c     error return in case of wrong input data
   25 ier=2
      return
      end
