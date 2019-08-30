!       get scaler product of two 3d vectors in cartesian coord.
!
      subroutine cscalerProd(a, b, c)

#include  "Zcoord.h"
      type(coord)::a 
      type(coord)::b
      real*8  c
!       c =  a . b
      c = 0.
      do i = 1, 3
         c = c + a%r(i) * b%r(i)
      enddo
      end

      
