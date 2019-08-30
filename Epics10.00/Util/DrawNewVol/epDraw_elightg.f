      subroutine epDraw_elightg(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

      include "Zelightg.h"

       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   a elightg in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line
                               ! dimension of p must be >+ (nvccl+2)*2
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  



      real*8 teta1, teta2

      integer Lu1, Lu2, Ll1, Ll2
!        starting location of upper ellipse front, back,
!                             lower ellipse front, back.
      integer Lfw1, Lfw2, Lbw1, Lbw2, Lfwv, Lbwv

      integer n1, n2
!         number of points for upper ellipse, lower ellipse    
      integer nearl,  nearr
!        number of front/back wall points at x<=b-t2
!          and x > b-t2

      integer i
      real*8  elightgf1, elightgf2, elightgf3, elightgf4
      external  elightgf1, elightgf2, elightgf3, elightgf4

      real*8 temp

!         set common values
      call epelightgset(comp) 

!       draw upper  curved surface
      teta1 = 270.d0
      teta2 = 360.d0
      Lu1 = 1
      call epdrawElps(b-t2, a-t1, 0.d0, teta1, teta2, p(Lu1), n1)
!        modify origin and plane
      do i = 0, n1-2
         p(Lu1+i)%z = p(Lu1+i)%y + a 
         p(Lu1+i)%y = elightgf3(p(Lu1+i)%x)
      enddo
      n = n1
!          p(n) is unchanged (seprator)
      Lu2 = n + 1
      do i = 0, n1-2
         n = n + 1
         p(n)%x = p(i+Lu1)%x
         p(n)%y = elightgf4( p(n)%x )
         p(n)%z = p(i+Lu1)%z
      enddo
      n = n + 1
      p(n)%x = gpsep
      n = n + 1
      p(n)%x = gpsep

!          draw lower curved surface

      Ll1 = n + 1
      call epdrawElps(b, a, 0.d0, teta1, teta2, p(Ll1), n2)

!         modify origin and plane ; search nearest point to
!         x = b-t2

      temp = 1.d30
      do i = 0,  n2-2
         p(i+Ll1)%z = p(i+Ll1)%y + a
         p(i+Ll1)%y = elightgf3(p(i+Ll1)%x)
         if( abs( p(i+Ll1)%x - b + t2 ) .lt. temp) then
            temp = abs( p(i+Ll1)%x - b + t2 )
            nearl = i + 1
         endif
      enddo

      nearr = n2-nearl   ! branch point must overlap so -1  is not needed

!        adjust the nearest point to b-t2
      p(nearl+Ll1-1)%x = b-t2
      p(nearl+Ll1-1)%y = -w1p
      p(nearl+Ll1-1)%z = zc

      n = n + n2

      Ll2 = n + 1
      do i = 0, n2-2
         n = n + 1
         p(n)%x = p(i+Ll1)%x
         p(n)%y = elightgf4( p(n)%x )
         p(n)%z = p(i+Ll1)%z
      enddo
      n = n + 1
      p(n)%x = gpsep
      n = n + 1
      p(n)%x = gpsep


!       wall at y = -w1  -w1p
      Lfw1 = n + 1
      do  i = 0,  n1-2
         n = n + 1
         p(n)%x = p(i + Lu1)%x
         p(n)%y = p(i + Lu1)%y
         p(n)%z = p(i + Lu1)%z
      enddo
      n = n + 1
      p(n)%x = gpsep

      Lfwv = n + nearl   ! branch poin overlaps for both regions so no +1
      do  i = 0,  n1 - 2
         n = n + 1
         if( i .lt. nearl ) then
            p(n)%x = p(i+Ll1)%x
            p(n)%y = p(i+Ll1)%y
            p(n)%z = p(i+Ll1)%z
         else
            p(n)%x = b-t2
            p(n)%y = -w1p
            p(n)%z =
     *      (p(i+Lu1)%z - p(i+Ll1)%z)/(p(i+Lu1)%x - p(i+Ll1)%x) *
     *      (b-t2 - p(i+Lu1)%x) + p(i+Lu1)%z
         endif
      enddo
      n = n + 1
      p(n)%x = gpsep
      n = n + 1
      p(n)%x = gpsep

!          wall at x > b-t2 
!
      Lfw2 = n + 1
      do i = 0, nearr -1
         n = n + 1
         p(n)%x = p(i+Lfwv)%x
         p(n)%y = p(i+Lfwv)%y
         p(n)%z = p(i+Lfwv)%z
      enddo
      n = n + 1
      p(n)%x = gpsep

      do i = 0, nearr-1
         n = n + 1
         p(n)%x =p(i+Ll1+nearl-1)%x
         p(n)%y =p(i+Ll1+nearl-1)%y
         p(n)%z =p(i+Ll1+nearl-1)%z
      enddo
      n = n + 1
      p(n)%x = gpsep
      n = n + 1
      p(n)%x = gpsep



!       wall at y = w2  w2p

      Lbw1 = n + 1
      do  i = 0,  n2-2
         n = n + 1
         p(n)%x = p(i + Lu2)%x
         p(n)%y = p(i + Lu2)%y
         p(n)%z = p(i + Lu2)%z
      enddo
      n = n + 1
      p(n)%x = gpsep

      Lbwv = n + nearl
      do  i = 0,  n2 - 2
         n = n + 1
         if( i .lt. nearl ) then
            p(n)%x = p(i+Ll2)%x
            p(n)%y = p(i+Ll2)%y
            p(n)%z = p(i+Ll2)%z
         else
            p(n)%x = b-t2
            p(n)%y = w2p
            p(n)%z =
     *      (p(i+Lu2)%z - p(i+Ll2)%z)/(p(i+Lu2)%x - p(i+Ll2)%x) *
     *      (b-t2 - p(i+Lu2)%x) + p(i+Lu2)%z
         endif
      enddo
      n = n + 1
      p(n)%x = gpsep
      n = n + 1
      p(n)%x = gpsep

!          wall at x > b-t2 
!
      Lbw2 = n + 1
      do i = 0, nearr -1
         n = n + 1
         p(n)%x = p(i+Lbwv)%x
         p(n)%y = p(i+Lbwv)%y
         p(n)%z = p(i+Lbwv)%z
      enddo
      n = n + 1
      p(n)%x = gpsep

      do i = 0, nearr-1
         n = n + 1
         p(n)%x =p(i+Ll2+nearl-1)%x
         p(n)%y =p(i+Ll2+nearl-1)%y
         p(n)%z =p(i+Ll2+nearl-1)%z
      enddo
      n = n + 1
      p(n)%x = gpsep
      n = n + 1
      p(n)%x = gpsep

      
!        left wall
      n = n + 1
      p(n)%x = 0.
      p(n)%y = -w1
      p(n)%z = 0.
      n = n +1
      p(n)%x = 0.
      p(n)%y = w2
      p(n)%z = 0.
      n = n +1
      p(n)%x = gpsep

      n = n + 1
      p(n)%x = 0.
      p(n)%y = -w1
      p(n)%z = t1
      n = n +1
      p(n)%x = 0.
      p(n)%y = w2
      p(n)%z = t1
      n = n +1
      p(n)%x = gpsep
      n = n +1
      p(n)%x = gpsep

!        top wall
      n = n + 1
      p(n)%x = b-t2
      p(n)%y = -w1p
      p(n)%z = a
      n = n + 1
      p(n)%x = b-t2
      p(n)%y = w2p
      p(n)%z = a
      n = n +1
      p(n)%x = gpsep

      n = n + 1
      p(n)%x = b
      p(n)%y = -w1p
      p(n)%z = a
      n = n + 1
      p(n)%x = b
      p(n)%y = w2p
      p(n)%z = a
      n = n +1
      p(n)%x = gpsep

      n = n +1
      p(n)%x = gpsep
      end
