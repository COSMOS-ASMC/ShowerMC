      subroutine epDraw_scyl(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe scyl

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      integer ir, ih, in1x, in1y, in1z, in2x, in2y, in2z
      integer maxz, minz
      parameter (ir = 1,  ih = 2,  in1x= 3, in1y=4, in1z=5)
      parameter (in2x= 6, in2y=7, in2z=8, maxz=9, minz=10)


      logical kdgtest

      integer n1, n2, i
      real*8 cost, sint, k2, r

      r = Volat( comp%vol + ir)
      n = 0
      call epdrawCcl( r, 0.d0,  thetamax, thetamin, p, n1)
!         adjust z
      do i = 1, n1 - 1
         cost = p(i)%x/r
         sint = p(i)%y/r
         p(i)%z = -r*( Volat( comp%vol + in1x)*cost  +
     *                 Volat( comp%vol + in1y)*sint )
     *            /Volat( comp%vol + in1z)
      enddo

      n = n + n1
      call epdrawCcl(r, Volat( comp%vol + ih), 
     *      thetamax, thetamin, p(n+1), n2)
      k2 = Volat( comp%vol + ih) * Volat( comp%vol + in2z)
      do i = n+1, n + n2 -1
         cost = p(i)%x/r
         sint = p(i)%y/r
         p(i)%z = (k2-r*( Volat( comp%vol + in2x)*cost  +
     *                    Volat( comp%vol + in2y)*sint ))
     *          /Volat( comp%vol + in2z)
      enddo
      n = n+ n2
      n = n + 1
      p(n)%x = gpsep
!   
      if(kdgtest(howcyl, 1)) then
         call epdrawCylEdg(p, n1, 0.d0, p(n+1), n2)
         n = n + n2 
      endif
      if(kdgtest(howcyl, 2)) then
         call epdrawCylEdg( p(n1+1), n1, Volat( comp%vol + ih),
     *       p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif
      end
