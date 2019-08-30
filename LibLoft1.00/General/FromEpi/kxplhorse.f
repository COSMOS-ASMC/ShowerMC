!
!        horse:                                  top view
!                                         Y     _______
!               side view                  | * |       | \
!              ------------                _*_ |       |____\_
!             /            \               |   |       | b'  |
!      h     /              \              |   |x0_____|     |
!           /                \             |  *    a'    \   | b
!          /                  \            | *           \   |
!         ---------------------            |*_______________\|  --> X
!                                                 a
!     
!      a: side length of the lower square >= 0
!      b: side length of  //           >=0  
!      h: hight  > 0
!   x0,y0: upper square origin
!     a': side length of the upper square   a'>=0
!     b': side length of the upper square   b'>=0.  (a=b=a'=b'=0 is not
!                                                    allowed)
      subroutine kxplhorse(a, b, h, x0, y0, ap, bp, 
     *                    pos, dir, length, icon, face)
       implicit none
!  #include "ZepTrackp.h"
#include "Zep3Vec.h"
! #include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
! #include "Zepdebug.h"

       real(8),intent(in):: a, b, h, x0, y0,  ap, bp  ! input homrse param.
       type(epPos)::  pos   ! input.  position.
       type(epDirec)::  dir  ! input. direction cosinse

       real(8),intent(out)::  length !  output length cm from pos to the boundary
       integer,intent(out)::icon  ! output 0: length obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume
       integer,intent(out):: face  ! when x-ing point exists, which face
!                          same as box case.
!                       face 1(bottom)  2: yz at a  3:xz at y-0  ; opposit face
!                       sum becomes 7.

       integer:: np, jcon
       
       type(epPos):: p1, p2, p3, p4

       real(8):: l, xpa(2)
       integer::  facea(2)
       np = 0
       if(a > 0. .and. b > 0.) then
!              bottom 
          p1 = epPos(0.d0, 0.d0, 0.d0)
          p2 = epPos(a, 0.d0, 0.d0)
          p3 = epPos(p2%x, b, 0.d0)
          p4 = epPos(0.d0, p3%y, 0.d0)
          call epxpLand4vp(p1, p2, p3, p4, pos, dir, l, icon, jcon)
          if(icon .le. 4 .and. l>0.) then
             np = np +1
             xpa(np) = l
             facea(np)= 1
          endif
       endif

       if(ap  > 0. .and. bp > 0.) then
!            top
          p1 = epPos(x0, y0, h)
          p2 = epPos(x0+ap, y0, h)
          p3 = epPos(p2%x, y0+bp, h)
          p4 = epPos(x0, p3%y, h)
          call epxpLand4vp(p1, p2, p3, p4, pos, dir, l, icon, jcon)

          if(icon .le. 4 .and. l > 0.) then
             np = np +1
             xpa(np) = l
             facea(np) = 6
          endif
          if(np .eq. 2) goto 100
       endif
!           side x-z 1       
       p1 = epPos(0.d0, 0.d0, 0.d0)
       p2 = epPos(a, 0.d0, 0.d0)
       p3 = epPos(x0+ap, y0, h)
       p4 = epPos(x0, y0, h)
       call epxpLand4vp(p1, p2, p3, p4, pos, dir, l, icon, jcon)

       if(icon .le. 4 .and. l>0.0) then
          np = np +1
          xpa(np) = l
          facea(np) = 3
       endif
       if(np .eq. 2) goto 100

!           side x-z 2
       p1 = epPos( 0.d0, b, 0.d0)
       p2 = epPos( x0, y0+bp, h)
       p3 = epPos(x0+ap, y0+bp, h)
       p4 = epPos(a, b, 0.d0)
       call epxpLand4vp(p1, p2, p3, p4, pos, dir, l, icon, jcon)

       if(icon .le. 4 .and. l > 0.) then
          np = np + 1
          xpa(np) = l
          facea(np) = 4
       endif
       if(np .eq. 2) goto 100
!           y-z 1       
       p1 = epPos(0.d0, 0.d0, 0.d0)
       p2 = epPos(x0, y0, h)
       p3 = epPos(x0, y0+bp, h)
       p4 = epPos(0.d0, b, 0.d0)
       call epxpLand4vp(p1, p2, p3, p4, pos, dir, l, icon, jcon)

       if(icon .le. 4 .and. l > 0. ) then
          np = np +1
          xpa(np) = l
          facea(np) = 5
       endif
       if(np .eq. 2) goto 100
!           y-z 2
       p1 = epPos(a, 0.d0, 0.d0)
       p2 = epPos(a, b, 0.d0)
       p3 = epPos(x0+ap, y0+bp, h)
       p4 = epPos(x0+ap, y0, h)
       call epxpLand4vp(p1, p2, p3, p4, pos, dir, l, icon, jcon)

       if(icon .le. 4 .and. l > 0. ) then
          np = np +1
          xpa(np) = l
          facea(np) = 2
       endif
 100   continue
       if(np == 2) then
!          outside
          if(xpa(1) >  0.d0  .and.  xpa(2) > 0.d0) then
!               outside
             icon = 1
             if( xpa(1) < xpa(2) ) then
                face = facea(1)
                length = xpa(1)
             else
                length = xpa(2)
                face = facea(2)
             endif
          else
             ! strage; 
             write(0,*) 'in  kxlphorse strange '
             write(0,*) 'a, b, h, x0, y0,  ap, bp'
             write(0,*) a, b, h, x0, y0,  ap, bp
             write(0,*) 'pos=',pos
             write(0,*) 'dir=',dir
             stop
          endif
       elseif( np == 1 ) then
          icon = 0
          length = xpa(1)
          face = facea(1)
       else
!           not cross with the volume
          icon = -1
       endif

       end
