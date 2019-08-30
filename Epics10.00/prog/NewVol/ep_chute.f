!
!  chute (slide).  
!      /.
!   z /   .                             
!   |.      .            x - a           z-h
!   | .        .       (------- )^2  + (-----)^2 =1
!   |   .         .         a             h
! h |  y  .            .
!   |  /      .              .
!   | /           .             b
!   |_____________________. x
!               a        
!     
!
!   Data format in config is:
!       ox oy oz  a b h [dir]
!
!      where (ox,oy,oz) is the origin in the world coord.
!           a,b,h>0  (we use h instead of c; c is used in Zglobal.h)
!
!      
      subroutine eprchute(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "chute"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer ia, ib, ih
       parameter( ia = 1,  ib = 2,  ih=3)

       real*8 a, b, h
!
!           read chute  data as 'new-*'
!           chute has 3 volume attributes and optional 6 direction cosines
!           of x,y (1~6)
!
!             next is mandatory
        call eprpst(comp, 3, 3, 1, 6)
!
        a = Volat( comp%vol + ia)
        b = Volat( comp%vol + ib)
        h = Volat( comp%vol + ih)

        if(a  .le. 0. .or. b .le. 0. .or. h .le. 0. ) then
           write(msg, *) comp%cn, '-th component: a=', a,
     *    ' b=', b, ' h=',h, 
     *    ' for chute;  invalid'
           call cerrorMsg(msg, 0)
        endif
       end
!   ***************************************
      subroutine epbchute(comp, pos, dir, length, icon)
       implicit none
#include "Zglobalc.h"
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"

!
!        find length to the boundary of 'comp' from 'pos'
!        with direction cos 'dir'
!     'pos' and 'dir' are given in this 'comp' local coordinate.
! 
 

       type(Component):: comp  ! input. you can extract volume parameters
                               !  by Volat( comp.vol + 1), etc
       type(epPos)::  pos   ! input.  position.
       type(epDirec)::  dir  ! input. direction cosinse

       real*8  length !  output length cm from pos to the boundary
       integer icon  ! output 0: length obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume

 
       integer ia, ib, ih
       parameter( ia = 1,  ib = 2,  ih=3)

       real*8 a, b, h

       real*8 aa, bb, cc, dd, f, xx,  zz, leng, x, y, z
       real*8 xpa(4)
       real*8 el

       integer nx, i

       real*8  eps
       save eps
       data eps/1.d-8/

       f(xx, zz) = ((xx-a)/a)**2 + ( (zz-h)/h)**2 - 1.0d0

       a = Volat( comp%vol + ia )
       b = Volat( comp%vol + ib )
       h = Volat( comp%vol + ih )
!          crossing point with a bounding box 
       call kxplbx(pos%x, pos%y, pos%z,  dir%x, dir%y, dir%z, 
     *    a, b,  h,  el,  icon)
       if(icon .eq. -1) then
!          no cross.
          return
       endif
!        x-point with the box 
       x = pos%x + el*dir%x
       y = pos%y + el*dir%y
       z = pos%z + el*dir%z



       if( abs(x) .lt. eps) then
!         at the left box wall.
          if(icon .eq. 1) then
             length = el
             return
          endif
       elseif( abs(z) .le. eps) then
!          at the bottom 
          if(icon .eq. 1) then
             length = el
             return
          endif
       endif
!         remaining cases are very comlex without systematic 
!         calc. which may sometimes be blue duck.
!
!         to obtain x-point with the  ellipse.
       aa = (dir%x/a)**2 +  (dir%z/h)**2
       bb = ((pos%x-a)*dir%x/a**2 +  (pos%z-h)*dir%z/h**2)
       cc = f(pos%x, pos%z)

       dd = bb**2 - aa*cc
 
       
       if(f(pos%x, pos%z) .ge. 0.) then
          if( (pos%z .le. h  .and. pos%z .ge. 0.) .and.
     *        (pos%x .le. a  .and. pos%x .ge. 0.) .and.
     *        (pos%y .le. b  .and. pos%y .ge. 0.)  ) then
!                inside
             icon = 0
          else
             icon = 1
          endif
       else
          icon = 1
       endif

       nx = 0
       x = pos%x + dir%x * el
       y = pos%y + dir%y * el
       z = pos%z + dir%z * el

       if( f(x, z) .ge. 0. ) then
          if(icon .eq. 1) then
             length = el
             return
          endif
          nx = nx +1
          xpa(nx) = el
       endif
 
       if(dd .ge. 0. .and. aa .gt. 0.d0) then
          dd = sqrt(dd)
          leng = (-bb+dd)/aa
          if(leng .ge. 0.) then
             x = pos%x + dir%x * leng
             y = pos%y + dir%y * leng
             z = pos%z + dir%z * leng
             if( x .le. a  .and. y .le. b .and. z .le. h .and.
     *           y .ge. 0.d0 ) then
                nx = nx + 1
                xpa(nx)= leng
             endif
          endif
          leng = (-bb -dd)/aa
          if( leng .ge. 0.) then
             x = pos%x + dir%x * leng
             y = pos%y + dir%y * leng
             z = pos%z + dir%z * leng
             if( x .le. a  .and. y .le. b .and. z .le. h .and.
     *           y .ge. 0.d0) then
                nx = nx + 1
                xpa(nx)= leng
             endif
          endif
       endif

       if(nx .eq. 0) then
          icon = -1
       else
          length = xpa(1)
          do i = 2, nx
             length = min(length, xpa(i))
          enddo
       endif
       end          

      subroutine epschute(comp, pos, icon)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
!           judge if a given 'pos' is inside 'comp'
!         
       type(Component)::  comp !input component
       type(epPos)::  pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside

 
       integer ia, ib, ih
       parameter( ia = 1,  ib = 2,  ih=3 )

       real*8 a, b, h

       a = Volat( comp%vol + ia )
       b = Volat( comp%vol + ib )
       h = Volat( comp%vol + ih )

       if(pos%z .gt. h) then
          icon = 1
       elseif(pos%z .lt. 0.) then
          icon = 1
       elseif(pos%x .lt. 0.) then
          icon = 1
       elseif(pos%x .gt. a ) then
          icon = 1
       elseif(pos%y .gt. b ) then
          icon =1 
       elseif(pos%y .lt. 0. ) then
          icon = 1
       elseif( ( (pos%x-a)/a )**2 + ( (pos%z-h)/h )**2 - 1.d0 
     *         .lt. 0.) then
          icon = 1
       else
          icon = 0
       endif
      end
!     **************************************
      subroutine epenvlpchute(comp, org, abc)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

!
!        give the envloping box of the component
!
       type(Component)::  comp  ! input.   component.
       type(epPos)::  org       ! output.  origin of the enveloping box
                               !          in local coord. 
       type(ep3Vec)::  abc      ! output.  a,b,c of the box

 
 
       integer ia, ib, ih
       parameter( ia = 1,  ib = 2,  ih=3)

       real*8 a, b, h


       a = Volat( comp%vol + ia)
       b = Volat( comp%vol + ib)
       h = Volat( comp%vol + ih)
       org%x = 0.
       org%y = 0.
       org%z = 0.
       abc%x = a
       abc%y = b
       abc%z = h
       NVTX = 6
       VTXx(1) = 0.
       VTXy(1) = 0.
       VTXz(1) = 0.

       VTXx(2) = a
       VTXy(2) = 0.
       VTXz(2) = 0

       VTXx(3) = a
       VTXy(3) = b
       VTXz(3) = 0

       VTXx(4) = 0.
       VTXy(4) = b
       VTXz(4) = 0.

       VTXx(5) = 0.
       VTXy(5) = 0.
       VTXz(5) = h

       VTXx(6) = 0.
       VTXy(6) = b
       VTXz(6) = h

      end
!     *************************************
      subroutine epatlocchute(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(3)
 
      integer i

      do i = 1, 3
         loc(i) = i
      enddo
      end

