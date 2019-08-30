      subroutine cmBremE(up, x)
!
!
      implicit none
!
      real*8 up !  input. Upsilon value. = Ee/Me *  B/Bcr
!
      real*8 x  !  output.  sampled fractional gamma energy

      real*8 cmBremI1, cmBremI2, sum, u, xs1
      integer nc

!         select 1st term or 2nd term in Eq.(5) of Brainerd and Petrosian.
!       APJ, 320, 1987
!
       xs1 =  cmBremI1(up, 0.d0)
       sum = xs1 + cmBremI2(up)
       call rndc(u)
       if(u .lt. xs1/sum) then
!               use first term
          call cmBremE1(up, x)
       else
          call cmBremE2(up, x, nc)
       endif
       end
