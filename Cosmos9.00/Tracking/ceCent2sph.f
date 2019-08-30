      subroutine ceCent2sph(a, bb)
!            rectangular coord to spherical coord
      implicit none

#include  "Zglobalc.h"

#include  "Zcoord.h"
      type(coord)::a
      type(coord)::bb
!
      type(coord)::b
      real*8 sint, cost, sinphi, cosphi
!
!      b.radius = sqrt(a.r(1)**2 + a.r(2)**2 + a.r(3)**2)
      b%r(3) = sqrt(a%r(1)**2 + a%r(2)**2 + a%r(3)**2)
!      if(b.radius .gt. 0.) then
      if(b%r(3) .gt. 0.) then
!           sint =sqrt( (a.r(1)/b.radius)**2 + (a.r(2)/b.radius)**2)
           sint =sqrt( (a%r(1)/b%r(3))**2 + (a%r(2)/b%r(3))**2)
!           cost = a.r(3)/b.radius
           cost = a%r(3)/b%r(3)
!           b.theta = atan2(sint, cost)*Todeg
           b%r(1) = atan2(sint, cost)*Todeg
           if(sint .ne. 0.d0) then
!              cosphi = a.r(1) / b.radius*sint
!              sinphi = a.r(2) / b.radius *sint
              cosphi = a%r(1) / b%r(3)*sint
              sinphi = a%r(2) / b%r(3) *sint
!              b.phi = atan2(sinphi, cosphi)*Todeg
              b%r(2) = atan2(sinphi, cosphi)*Todeg
           else
!              b.phi = 0.d0
              b%r(2) = 0.d0
           endif   
       else
!           b.theta = 0.d0
!           b.phi = 0.d0
           b%r(1) = 0.d0
           b%r(2) = 0.d0
       endif
       b%sys='sph'
       bb = b
       end
