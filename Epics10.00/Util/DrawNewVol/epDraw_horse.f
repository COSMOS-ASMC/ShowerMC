      subroutine epDraw_horse(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(42)     ! output. (x,y,z) to describe
                               !  horse in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      logical kdgtest

      integer ia, ib, ih, ix0, iy0,  iap, ibp
      parameter( ia = 1,  ib = 2,  ih = 3,  ix0=4, iy0=5,
     *            iap=6, ibp=7 )

       real*8 a, b, h, x0, y0,  ap, bp
!
       n = 0
       a = Volat( comp%vol + ia)
       b = Volat( comp%vol + ib)
       h = Volat( comp%vol + ih)
       x0= Volat( comp%vol + ix0)
       y0= Volat( comp%vol + iy0)
       ap = Volat( comp%vol + iap)
       bp = Volat( comp%vol + ibp)
!           use the same surface idex as the box.
!
!          |                                  1 x-y z=0
!          |     ******************           6 x-y z=c
!          |   * |          5    **           2 x-z y=0           
!          | *   |   6/        *  *           5 x-z y=b
!        c |********/ ********* 4 *           3 y-z x=0
!          | 3 b |/+++++++++++*+++*           4 y-z x=a 
!          |    /    1        *  *
!          |   /              * *
!          | /    2           *
!          |------------------------------
!                             a
!
 
!          we follow the  logic used in box drawing
      if( kdgtest(how, 1) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
     
      if( kdgtest(how, 2) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = x0
         p(n)%y = y0
         p(n)%z = h

         n = n + 1
         p(n)%x = x0+ap
         p(n)%y = y0
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      if( kdgtest(how, 3) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = x0
         p(n)%y = y0
         p(n)%z = h

         n = n + 1
         p(n)%x = x0
         p(n)%y = y0+bp
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      if( kdgtest(how, 4) )then
         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = x0+ap
         p(n)%y = y0
         p(n)%z = h

         n = n + 1
         p(n)%x = x0+ap
         p(n)%y = y0+bp
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif


      if( kdgtest(how, 5) )then
         n = n + 1
         p(n)%x = 0
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = x0
         p(n)%y = y0+bp 
         p(n)%z = h

         n = n + 1
         p(n)%x = x0+ap
         p(n)%y = y0+bp
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif


      if( kdgtest(how, 6) )then
         n = n + 1
         p(n)%x = x0
         p(n)%y = y0
         p(n)%z = h

         n = n + 1
         p(n)%x = x0+ap
         p(n)%y = y0
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = x0
         p(n)%y = y0+bp
         p(n)%z = h

         n = n + 1
         p(n)%x = x0+ap
         p(n)%y = y0+bp
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
      end
