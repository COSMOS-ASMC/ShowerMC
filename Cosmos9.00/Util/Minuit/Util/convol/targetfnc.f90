subroutine targetfnc(npar, gin, f, paramx, iflag)
  implicit none
  include "Zlatfit.h"
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
!
      integer i
      real*8 fval,  xx
      real*8 temp
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
      maxdiff = 0.
      answer = 0.
      pw = paramx(4)
      do  i= 1, npoint
         if(badindex(i) .lt. 0) cycle
         fval =   paramx(1)/
!     *   x(i)**(paramx(2) + paramx(3)* log(x(i)) )
     *   x(i)**(paramx(2) + paramx(3)* x(i)**pw )
!         chisq= chisq + (fval-y(i))**2
!         if(y(i) .lt. 3.) then
!            chisq= chisq + ( fval/y(i) - 1.0 )**2
!         else
!ccc           chisq= chisq + ( fval/y(i) - 1.0 )**2
!          chisq= chisq + ( y(i)/fval - 1.0 )**2  NG
            
         if(y(i)%gt. 0.) then
            temp = ( log(fval/y(i)))**2

!            if(temp .gt. limit) answer = answer + 1.

            if(temp .gt. maxdiff)then
               maxdiff= temp
               maxindex = i
            endif
!               chisq= chisq + ( y(i)/fval + fval/y(i) - 2.0 )**2 
            chisq= chisq + temp
         endif
!           chisq= chisq + ( log(fval/y(i)))**2  NG
!cc              next is NG
!cc            chisq= chisq + ( fval - y(i) )**2/y(i)

!         endif
      enddo
!      if(limit .gt. 0.)  then
!         if(maxdiff .gt. limit) then
!            chisq = chisq - maxdiff
!         endif
!      endif
      f = chisq

      if(iflag .eq. 3) then
!         do i = 1, npoint
!            fval =   paramx(1)*x(i)**(paramx(2) +
!            write(*,*) x(i), y(i), fval
!         enddo
         do i = 1, npar
            oparam(i) = paramx(i)
         enddo
      endif
      end
