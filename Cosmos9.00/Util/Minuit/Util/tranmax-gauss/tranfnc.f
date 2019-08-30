      subroutine tranfnc(npar,gin,f,paramx,iflag)
      implicit none
      include "Zfit.h"
* npar: input. number of current variable  parameters
* gin:  optional output of gradient
*  f:  output. kaisq
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
!
      integer i
      real*8 ff,  xx
      real*8  a, s, x0
      save
      if (iflag .eq. 1) then
!         write(0,*) 'current npar=',npar,  ' npoint=',npoint
!         write(0,*) ' param=',(paramx(i), i=1,npar)
      endif
      if(iflag .eq. 2 ) then
         write(0,*) ' no grad computed'
         stop
      endif
!        compute f

      chisq = 0.
!         f(x)=  a*(x/x0)**b * exp(-c*(x/x0)**d)
!            z = x/x0
!         f'(z) = ab z**(b-1)* exp(-c*z**d) - acdz**b*exp(-cz**d)*z**(d-1)
!          =0 
!         --> b - cdz**d=0 ;  z=(b/cd)**(1./d);
!               x = x0 (b/cd)**(1./d)
! ccc not used:        if we use x0=max pos. cd=b
! ccc                  fmax = a*exp(-c)   

      do  i= 1, npoint
         a = paramx(1)
         s = paramx(2)
         x0 = paramx(3)
         ff = a*exp(-((x(i)-x0)/s)**2/2.)
!         chisq= chisq + (ff-y(i))**2
!         if(y(i) .lt. 3.) then
!            chisq= chisq + ( ff/y(i) - 1.0 )**2
!         else
!            chisq= chisq + ( y(i)/ff - 1.0 )**2
!            chisq= chisq + ( ff/y(i) - 1.0 )**2
            chisq= chisq + (ff-y(i))**2/y(i)

!         endif
      enddo
      f = chisq

      if(iflag .eq. 3) then
!         do i = 1, npoint
!            ff =   paramx(1)*x(i)**(paramx(2) +
!            write(*,*) x(i), y(i), ff
!         enddo
         do i = 1, npar
            oparam(i) = paramx(i)
         enddo
      endif
      end
