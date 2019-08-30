      subroutine epdrawbox(copm, n, from, to)
      implicit none
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp
      integer n
       type(epPos):: from(*)
       type(epPos)::  to(*)

      from(1)%x = 0.
      from(1)%y = 0.
      from(1)%z = 0.

      to(1)%x = comp%vol(boxa)
      to(1)%y = 0.
      to(1)%z= 0.

      from(2) = to(1)
      to(2) =  to(1)
      to(2)%y = comp%vol(boxb)

      from(3) = to(2)
      to(3) = to(2)
      to(3)%x = 0.

      from(4) = to(3)
      to(4)= from(1)
!
      from(5) = from(1)
      from(5)%z= comp%vol(boxc)
      to(5) = to(1)
      to(5)%z = comp%vol(boxc)

      from(6) = to(5)
      to(6) = to(2)
      to(6)%z = comp%vol(boxc)
      
      from(7) = to(6)
      to(7) = to(3)
      to(7)%z = comp%vol(boxc)

      from(8) = to(7)
      to(8) = to(4)
      to(8)%z = comp%vol(boxc)
      end

      

