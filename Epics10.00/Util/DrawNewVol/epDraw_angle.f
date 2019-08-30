      subroutine epDraw_angle(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !  angle in local coordnate.
                              ! (x,y,z)= gpsep is a separator
                              ! to be converted to a blank line

      integer  n              ! output.  number of (x,y,z) data
                              ! put in p.  



       integer ia, ib, ic, id, ih, ita, itb, inow
       parameter( ita=1,  itb = 2,  ih=3,   ia=4,
     *   ib=5,  ic=6, id=7, inow=8)

       real*8 a, b, c, d, h
!
       n = 0
       a = Volat( comp%vol + ia)
       b = Volat( comp%vol + ib)
       h = Volat( comp%vol + ih)
       c = Volat( comp%vol + ic)
       d = Volat( comp%vol + id)

!           

       n = n + 1
       p(n)%x = 0.
       p(n)%y = 0.
       p(n)%z = 0.

       n = n + 1
       p(n)%x = a
       p(n)%y = 0.
       p(n)%z = 0.


       n = n + 1
       p(n)%x = a
       p(n)%y = d
       p(n)%z = 0.

       n = n + 1
       p(n)%x = c
       p(n)%y = d
       p(n)%z = 0.

       n = n + 1
       p(n)%x = c
       p(n)%y = b
       p(n)%z = 0.

       n = n + 1
       p(n)%x = 0
       p(n)%y = b
       p(n)%z = 0.

       n = n + 1
       p(n)%x = 0
       p(n)%y = 0
       p(n)%z = 0


!----------------  upper part
       n = n+1
       p(n)%x=gpsep

       n = n + 1
       p(n)%x = 0.
       p(n)%y = 0.
       p(n)%z = h

       n = n + 1
       p(n)%x = a
       p(n)%y = 0.
       p(n)%z = h


       n = n + 1
       p(n)%x = a
       p(n)%y = d
       p(n)%z = h

       n = n + 1
       p(n)%x = c
       p(n)%y = d
       p(n)%z = h

       n = n + 1
       p(n)%x = c
       p(n)%y = b
       p(n)%z = h

       n = n + 1
       p(n)%x = 0
       p(n)%y = b
       p(n)%z = h

       n = n + 1
       p(n)%x = 0
       p(n)%y = 0
       p(n)%z = h

       n= n+ 1
       p(n)%x=gpsep
       n= n+ 1
       p(n)%x=gpsep
!          floor
       n = n + 1
       p(n)%x = 0.
       p(n)%y = b
       p(n)%z = 0.
       
       n = n + 1
       p(n)%x = 0.
       p(n)%y = 0.
       p(n)%z = 0.

       n = n + 1
       p(n)%x = a
       p(n)%y = 0.
       p(n)%z = 0.

       n= n+ 1
       p(n)%x=gpsep
       
       n = n + 1
       p(n)%x = c
       p(n)%y = b
       p(n)%z = 0.

       n = n + 1
       p(n)%x = c
       p(n)%y = d
       p(n)%z = 0.

       n = n + 1
       p(n)%x = a
       p(n)%y = d
       p(n)%z = 0.

       n= n+ 1
       p(n)%x=gpsep
       n= n+ 1
       p(n)%x=gpsep

!          ceil
       n = n + 1
       p(n)%x = 0.
       p(n)%y = b
       p(n)%z = h
       
       n = n + 1
       p(n)%x = 0.
       p(n)%y = 0.
       p(n)%z = h

       n = n + 1
       p(n)%x = a
       p(n)%y = 0.
       p(n)%z = h

       n= n+ 1
       p(n)%x=gpsep
       
       n = n + 1
       p(n)%x = c
       p(n)%y = b
       p(n)%z = h

       n = n + 1
       p(n)%x = c
       p(n)%y = d
       p(n)%z = h

       n = n + 1
       p(n)%x = a
       p(n)%y = d
       p(n)%z = h

       n= n+ 1
       p(n)%x=gpsep
       n= n+ 1
       p(n)%x=gpsep

      end
