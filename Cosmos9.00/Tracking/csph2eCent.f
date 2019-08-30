      subroutine csph2eCent(a, b)
!         spherical to   rectangular coord
      implicit none
#include  "Zglobalc.h"
#include  "Zcoord.h"
      type(coord)::a
      type(coord)::b
!
      type(coord)::temp
      real*8  t
!
!      t =a.theta*Torad
      t =a%r(1)*Torad
!      temp.r(1) = a.radius*sin(t)*cos(a.phi*Torad)
!      temp.r(2) = a.radius*sin(t)*sin(a.phi*Torad)
!      temp.r(3) = a.radius* cos(t)
      temp%r(1) = a%r(3)*sin(t)*cos(a%r(2)*Torad)
      temp%r(2) = a%r(3)*sin(t)*sin(a%r(2)*Torad)
      temp%r(3) = a%r(3)* cos(t)
      temp%sys='xyz'
      b = temp
      end
