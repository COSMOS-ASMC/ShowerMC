      subroutine intePrim2(comp, i1, i2, ans)
      implicit none
#include "Zptcl.h"
#include "Zprimary.h"
!          primary flux at E       
      type (component)::comp
      integer i1, i2
      real*8 ans

      real*8 sum

      integer i
      sum = 0.

      do i = i1, i2
         if(comp%cut .ne. 0.) then
            if(comp%emin .le. comp%energy(i)) then
               sum =  comp%norm_inte(i)
               goto 100
            endif
         else
            sum = comp%norm_inte(i)
            goto 100
         endif
      enddo
 100  continue
      ans = sum * comp%inte_value
      end
