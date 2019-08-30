!       get vector product of two 3d vectors in cartesian coord.
!
      subroutine cvecProd(a, b, c)
!----      include 'Zcoord.h'
#include  "Zcoord.h"
      type(coord)::a
      type(coord)::b
      type(coord)::c
!       c <-  a x b
      c%r(1) = a%r(2) * b%r(3) - a%r(3) * b%r(2)
      c%r(2) = a%r(3) * b%r(1) - a%r(1) * b%r(3)
      c%r(3) = a%r(1) * b%r(2) - a%r(2) * b%r(1)
      end

      
