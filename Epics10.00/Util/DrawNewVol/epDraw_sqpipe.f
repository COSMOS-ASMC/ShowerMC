      subroutine epDraw_sqpipe(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)      ! output. (x,y,z) to describe
                               !  sqpipe in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      logical kdgtest
 
      integer ia, ib, ig, ix0, iy0, iap, ibp
      parameter( ia = 1,  ib = 2,  ig=3, ix0=4, iy0= 5, iap=6, ibp=7)

      real*8 a, b, g, x0, y0, ap, bp
       type(epPos)::  p1, p2, p3, p4
       type(epPos)::  q1, q2, q3, q4
       type(epPos)::  r1, r2, r3, r4
       type(epPos)::  s1, s2, s3, s4

      a = Volat( comp%vol + ia)
      b = Volat( comp%vol + ib)
      g = Volat( comp%vol + ig)
      x0 = Volat( comp%vol + ix0)
      y0 = Volat( comp%vol + iy0)
      ap = Volat( comp%vol + iap)
      bp = Volat( comp%vol + ibp)

      p1%x = 0.
      p1%y = 0.
      p1%z = 0.

      p2%x = a
      p2%y = 0.
      p2%z = 0.

      p3%x = a
      p3%y = b
      p3%z = 0.

      p4%x = 0.
      p4%y = b
      p4%z = 0.
      
      q1%x = 0.
      q1%y = 0.
      q1%z = g

      q2%x = a
      q2%y = 0.
      q2%z = g

      q3%x = a
      q3%y = b
      q3%z = g

      q4%x = 0.
      q4%y = b
      q4%z = g

      r1%x = x0
      r1%y = y0
      r1%z = 0.

      r2%x = x0+ap
      r2%y = y0
      r2%z = 0.

      r3%x = x0+ap
      r3%y = y0+bp
      r3%z = 0.

      r4%x = x0
      r4%y = y0+bp
      r4%z = 0.

      s1%x = x0
      s1%y = y0
      s1%z = g

      s2%x = x0+ap
      s2%y = y0
      s2%z = g

      s3%x = x0+ap
      s3%y = y0+bp
      s3%z = g

      s4%x = x0
      s4%y = y0+bp
      s4%z = g

 

!          we follow the  logic used in box drawing
!           1 is x-y at z=0
!           6 is x-y at z=c
!           2 is x-z at y=0
!           5 is x-z at y=b
!           3 is y-z at x=0
!           4 is y-z at x=a

      if( kdgtest(how, 1) )then
!           1 is x-y at z=0
         if(r1%y .gt. p1%y) then
            n = n + 1
            p(n) = r1
            p(n)%y = p1%y
            n = n + 1
            p(n) = p2
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = r1
            n = n + 1
            p(n) = p2
            p(n)%y = r2%y
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         endif
         if(r2%x .lt. p2%x ) then
            n = n + 1
            p(n) = r2
            p(n)%x = p2%x
            n = n + 1
            p(n) = p3
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = r2
            n = n + 1
            p(n) = p3
            p(n)%x = r3%x
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         endif
         if(r3%y .lt. p3%y) then
            n = n + 1
            p(n) = r3
            p(n)%y = p3%y
            n = n + 1
            p(n) = p4
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = r3
            n = n + 1
            p(n) = p4
            p(n)%y = r4%y
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         endif
         if(r4%x .gt. p4%x) then
            n = n + 1
            p(n) = r4
            p(n)%x = p4%x
            n = n + 1
            p(n) = p1
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = r4
            n = n + 1
            p(n) = p1
            p(n)%x = r1%x
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         endif
      endif

      if( kdgtest(how, 2) )then
!          2      is x-z at y=0
         if(r1%y .gt. p1%y ) then
!               outer wall
            n = n + 1
            p(n)= p1
            n = n + 1
            p(n) = p2
            n = n + 1
            p(n)%x = gpsep
            
            n = n + 1
            p(n) = q1
            n = n + 1
            p(n) = q2
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
!                inner wall
            n = n + 1
            p(n)= r1
            n = n + 1
            p(n) = r2
            n = n + 1
            p(n)%x = gpsep
            
            n = n + 1
            p(n) = s1
            n = n + 1
            p(n) = s2
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n)%x = gpsep
         else
            if(r1%x .gt. p1%x) then
!                left 
               n = n + 1
               p(n)= p1
               n = n + 1
               p(n) = r1
               p(n)%y = 0.
               n = n + 1
               p(n)%x = gpsep
            
               n = n + 1
               p(n) = q1
               n = n + 1
               p(n) = s1
               p(n)%y = 0.
               n = n + 1
               p(n)%x = gpsep

               n = n + 1
               p(n)%x = gpsep
             endif
             if(r2%x .lt. p2%x) then
!                right
               n = n + 1
               p(n)= r2
               p(n)%y = p2%y
               n = n + 1
               p(n) = p2
               n = n + 1
               p(n)%x = gpsep
            
               n = n + 1
               p(n) = s2
               p(n)%y = q2%y
               n = n + 1
               p(n) = q2
               n = n + 1
               p(n)%x = gpsep

               n = n + 1
               p(n)%x = gpsep
             endif
         endif
      endif

      if( kdgtest(how, 3) )then
!           3 is y-z at x=0
         if(r1%x .gt. p1%x) then
            n = n + 1
            p(n) = p4
            n = n + 1
            p(n) = p1
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = q4
            n = n + 1
            p(n) = q1
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
!                    inner wall
            n = n + 1
            p(n) = r4
            n = n + 1
            p(n) = r1
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = s4
            n = n + 1
            p(n) = s1
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         else
            if(r4%y .lt. p4%y) then
               n = n + 1
               p(n) = p4
               n = n + 1
               p(n) = r4
               p(n)%x = p4%x
               n = n + 1
               p(n)%x = gpsep

               n = n + 1
               p(n) = q4
               n = n + 1
               p(n) = s4
               p(n)%x = p4%x
               n = n + 1
               p(n)%x = gpsep
               n = n + 1
               p(n)%x = gpsep
            endif
            if(r1%y .gt. p1%y) then
               n = n + 1
               p(n) = r1
               p(n)%x = p1%x
               n = n + 1
               p(n) = p1
               n = n + 1
               p(n)%x = gpsep

               n = n + 1
               p(n) = s1
               p(n)%x = p1%x
               n = n + 1
               p(n) = q1
               n = n + 1
               p(n)%x = gpsep
               n = n + 1
               p(n)%x = gpsep
            endif
         endif
      endif

      if( kdgtest(how, 4) )then
!                  y-z at a
         if( r2%x .lt. p2%x ) then
!                 outer wall
            n = n + 1
            p(n) = p2
            n = n + 1
            p(n) = p3
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = q2
            n = n + 1
            p(n) = q3
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
!                   inner wall
            n = n + 1
            p(n) = r2
            n = n + 1
            p(n) = r3
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = s2
            n = n + 1
            p(n) = s3
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         else
            if(r2%y .gt. p2%y) then
!                  y=0 to y0
               n = n + 1
               p(n) = p2
               n = n + 1
               p(n) = r2
               p(n)%x = p2%x
               n = n + 1
               p(n)%x = gpsep

               n = n + 1
               p(n) =q2
               n = n + 1
               p(n) = s2
               p(n)%x = q2%x
               n = n + 1
               p(n)%x = gpsep
               n = n + 1
               p(n)%x = gpsep
            endif
            if(r3%y .lt. p3%y) then
               n = n + 1
               p(n) = r3
               p(n)%x = p3%x
               n = n + 1
               p(n) = p3
               n = n + 1
               p(n)%x = gpsep

               n = n + 1
               p(n) = s3
               p(n)%x = q3%x
               n = n + 1
               p(n) = q3
               n = n + 1
               p(n)%x = gpsep
               n = n + 1
               p(n)%x = gpsep
            endif
         endif               
      endif

      if( kdgtest(how, 5) )then
!           5       is x-z at y=b
         if(r3%y .lt. p3%y) then
!             outer wall
            n = n + 1
            p(n) = p3
            n = n + 1
            p(n) = p4
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = q3
            n = n + 1
            p(n) = q4
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
!             inner wall
            n = n + 1
            p(n) = r3
            n = n + 1
            p(n) = r4
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = s3
            n = n + 1
            p(n) = s4
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         else
            if(r3%x .lt. p3%x) then
!                   x=a to x0+ap
               n = n + 1
               p(n) = p3
               n = n + 1
               p(n) = r3
               p(n)%y = p3%y
               n = n + 1
               p(n)%x = gpsep

               n = n + 1
               p(n) = q3
               n = n + 1
               p(n) = s3
               p(n)%y = q3%y
               n = n + 1
               p(n)%x = gpsep
               n = n + 1
               p(n)%x = gpsep
            endif               
            if(r4%x .gt. p4%x) then
!                  x=0 to x0
               n = n + 1
               p(n) = r4
               p(n)%y = p4%y
               n = n + 1
               p(n) = p4
               n = n + 1
               p(n)%x = gpsep

               n = n + 1
               p(n) = s4
               p(n)%y = q4%y
               n = n + 1
               p(n) = q4
               n = n + 1
               p(n)%x = gpsep
               n = n + 1
               p(n)%x = gpsep
            endif               
         endif
      endif


      if( kdgtest(how, 6) ) then
!           6 is x-y at top
         if( s1%y .gt. q1%y ) then
            n = n + 1
            p(n) = s1
            p(n)%y = q1%y
            n = n + 1
            p(n) = q2
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = s1
            n = n + 1
            p(n) = q2
            p(n)%y = s2%y
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         endif
         if(s2%x .lt. q2%x) then
            n = n + 1
            p(n) = s2
            p(n)%x = q2%x
            n = n + 1
            p(n) = q3
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = s2
            n = n + 1
            p(n) = q3
            p(n)%x = s3%x
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         endif
         if( s3%y .lt. q3%y ) then
            n = n + 1
            p(n) = s3
            p(n)%y = q3%y
            n = n + 1
            p(n) = q4
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = s3
            n = n + 1
            p(n) = q4
            p(n)%y =s4%y
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         endif
         if(s4%x .gt. q4%x) then
            n = n + 1
            p(n) =  s4
            p(n)%x = q4%x
            n = n + 1
            p(n) = q1
            n = n + 1
            p(n)%x = gpsep

            n = n + 1
            p(n) = s4
            n = n + 1
            p(n) = q1
            p(n)%x = s1%x
            n = n + 1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         endif
      endif
      end
