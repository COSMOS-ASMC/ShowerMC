!         3 D vector's direction is obtained
      subroutine c3DV2DDCos(a,  b, len)
      implicit none
#include "Zcoord.h"
      type(coord)::a   ! input coord (only 3 vect is used)
      type(coord)::b   ! output coord.  in direction cos. b can be a.
      real*8 len         ! output. length of a.
      
         len = sqrt( a%r(1)**2 + a%r(2)**2 + a%r(3)**2)
         if(len .eq. 0.) then
            b%r(1) = 0.
            b%r(2) = 0.
            b%r(3) = 1.
         else
            b%r(1) = a%r(1)/len
            b%r(2) = a%r(2)/len
            b%r(3) = a%r(3)/len
         endif
         end
