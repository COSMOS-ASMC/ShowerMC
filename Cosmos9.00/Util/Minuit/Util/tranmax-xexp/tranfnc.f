      subroutine tranfnc(npar,gin,f,paramx,iflag)
      implicit none
      include "Zfit.h"
* npar: input. number of current variable  parameters
* gin:  optional output of gradient
*  f:  output. function value to be minimized 
* paramx: input. vector of const and variable parameters
*  iflag: input. depending on this value, what to do 
*        is determined.
*       1--> some preparation for computing f
*       2--> compute gin
*       3--> fittng finished
*       for all other cases: we  must compute f
*
      integer  npar
      real*8   gin(*)
      real*8   f
      real*8   paramx(*)
      integer iflag
c
      integer i
      real*8 fval,  xx
      real*8  a, b, c, d, z, x0
      save
      if (iflag .eq. 1) then
c         write(0,*) 'current npar=',npar,  ' npoint=',npoint
c         write(0,*) ' param=',(paramx(i), i=1,npar)
      endif
      if(iflag .eq. 2 ) then
         write(0,*) ' no grad computed'
         stop
      endif
c        compute f

      chisq = 0.
c         f(x)=  a*(x/x0)**b * exp(-c*(x/x0)**d)
c            z = x/x0
c         f'(z) = ab z**(b-1)* exp(-c*z**d) - acdz**b*exp(-cz**d)*z**(d-1)
c          =0 
c         --> b - cdz**d=0 ;  z=(b/cd)**(1./d);
c               x = x0 (b/cd)**(1./d)
c ccc not used:        if we use x0=max pos. cd=b
c ccc                  fmax = a*exp(-c)   

      do  i= 1, npoint
         x0 = paramx(2)
         z = x(i)/x0
         a = paramx(1)
         b =  paramx(3)
         c = paramx(4)
         d =  paramx(5)
         fval =   a*z**b*exp(-c*z**d)
c         chisq= chisq + (fval-y(i))**2
c         if(y(i) .lt. 3.) then
c            chisq= chisq + ( fval/y(i) - 1.0 )**2
c         else
c            chisq= chisq + ( y(i)/fval - 1.0 )**2
c            chisq= chisq + ( fval/y(i) - 1.0 )**2
            chisq= chisq + (fval-y(i))**2/y(i)
c         endif
      enddo
      f = chisq

      if(iflag .eq. 3) then
c         do i = 1, npoint
c            fval =   paramx(1)*x(i)**(paramx(2) +
c            write(*,*) x(i), y(i), fval
c         enddo
         do i = 1, npar
            oparam(i) = paramx(i)
         enddo
      endif
      end
