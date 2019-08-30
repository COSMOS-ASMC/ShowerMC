      subroutine episoAngle( dir ) 
      implicit none
      real(8)::  dir(3)

      real*8  cost, cs, sn, sint

      call rndc(cost)
      cost = 2.0d0*cost-1.0d0
      call kcossn(cs,sn)
      sint = sqrt(1.d0-cost**2)
      dir(1)= cs*sint
      dir(2)= sn*sint
      dir(3)= cost
      end
