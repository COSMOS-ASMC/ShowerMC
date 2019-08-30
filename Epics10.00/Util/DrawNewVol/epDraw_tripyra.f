      subroutine epDraw_tripyra(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(42)     ! output. (x,y,z) to describe
                               !  tripyra in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      logical kdgtest

      integer ia, ib, ic, id, ie, ih
      parameter( ia = 1,  ib = 2,  ic=3, id=4, ie=5, ih = 6)
       real*8 a, b, c, d, e,  h
!
       n = 0
       a = Volat( comp%vol + ia)
       b = Volat( comp%vol + ib)
       c = Volat( comp%vol + ic)
       d = Volat( comp%vol + id)
       e = Volat( comp%vol + ie)
       h = Volat( comp%vol + ih)
!
!           use the same surface idex as the box.
!
!          |                                  1 x-y z=0
!          |     ******************           6 not exists
!          |   * |          5    **           2 x-z y=0           
!          | *   |   6/        *  *           5 not exists
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
         p(n)%x = b
         p(n)%y = c
         p(n)%z = 0.

         n = n + 1
         p(n)%x = b
         p(n)%y = c
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
         p(n)%x = d
         p(n)%y = e
         p(n)%z = h

         n = n + 1
         p(n)%x = d 
         p(n)%y = e
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
         p(n)%x = b
         p(n)%y = c
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = d
         p(n)%y = e
         p(n)%z = h

         n = n + 1
         p(n)%x = d
         p(n)%y = e
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
         p(n)%x = b
         p(n)%y = c
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = d
         p(n)%y = e
         p(n)%z = h

         n = n + 1
         p(n)%x = d
         p(n)%y = e
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
      end

