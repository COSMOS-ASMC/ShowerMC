!
!  cut cylinder
!                                        
!  bottom circle center is at (0,0,0).  hight is directed to Z.
!                        +  ea
!                      +  
!                    + \*
!                       \  
!                 +      \ *
!                       * \
!                    *     \
!             +   *         \*
!              *          ^     sa
!           *_____________________   
!
!   Data format in config is:
!       ox oy oz  r  h  sa  ea (optional dir)
!
!      where (ox,oy,oz) is the origin in the world coord.
!            r: radius of the cylinder  cm
!            h: height of the //        cm
!           sa: starting angle (deg)  0 is the x-axis. counter clock wise.
!           ea: ending angle (deg).  sa=0 ea=360 means cyl.
!              sa may be > ea.
!      
      module modccyl
       integer,parameter::ir = 1,  ih = 2,  isa=3, iea=4,
     *        ix1=5, iy1=6, ix2=7, iy2=8
      end   module modccyl

      subroutine eprccyl(comp)
      use modccyl
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "ccyl"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 

       real*8 r, h, sa, ea
!
!           read cut cylinder data as 'new-*'
!           ccyl has 4 volume attributes and the direction cosines
!           of the  h (1~6)
!
!             next is mandatory
        call eprpst(comp, 4, 8, 1, 6)
!
!           next is optional
!           check some values
        r = Volat( comp%vol + ir)
        h = Volat( comp%vol + ih)
        sa= Volat( comp%vol + isa)
        ea= Volat( comp%vol + iea)
        if(r  .le. 0. .or. h .le. 0) then
           write(msg, *) comp%cn, '-th component: r=', r,
     *    ' h=', h, ' for ccyl;  invalid'
           call cerrorMsg(msg, 0)
        endif
!              to judge the pos. whether it is on the right
!         or legt side of line from x1(sa) to x2(ea)
        Volat( comp%vol + ix1) = r*cos( sa*Torad )
        Volat( comp%vol + iy1) = r*sin( sa*Torad )
        Volat( comp%vol + ix2) = r*cos( ea*Torad )
        Volat( comp%vol + iy2) = r*sin( ea*Torad )
       end
!   ***************************************
      subroutine epbccyl(comp, pos, dir, length, icon)
      use modccyl
       implicit none
#include "Zglobalc.h"
#include "ZepTrackp.h"
#include "ZepDirec.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"


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

!
       integer where, kcon, jcon
       type(epPos)::  p1, p2, p3, p4
       real*8 f, h, r, l, x, y, x1, y1, x2, y2
       type(epPos)::  xp
!              this gives wrong results
!       f(x,y) = (y-p1.y)*(p2.x-p1.x) -(p2.y-p1.y)*(x-p1.x)
       f(x,y) = (y-y1)*(x2 - x1) -(y2-y1)*(x-x1)

       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)

       call kxplcy(pos%x, pos%y, pos%z, dir%x, dir%y, dir%z, r, h,
     *  length,    icon,  where)
!  output:
!     length: x-ssing point is at pos + length*dir ( el>=0)
!   icon : output. 0 the point is in the cyl. length is obtained
!                  1 the point is out side of the cyl. length is
!                     obtained.
!                 -1 no x-ing point
!   where: output. 1  x-ing point is on x-y  top plane.
!                  2  //             on the side.
!                  6  //             on      bottom.
!                 -1  no x-ing point
       

       if( icon .eq. -1 ) then
!               no x-point
       else
          p1%x = Volat( comp%vol + ix1)
          p1%y = Volat( comp%vol + iy1)
          p2%x = Volat( comp%vol + ix2)
          p2%y = Volat( comp%vol + iy2)
          x1 = p1%x
          y1 = p1%y
          x2 = p2%x
          y2 = p2%y
          xp%x = pos%x + length*dir%x
          xp%y = pos%y + length*dir%y
          if(icon .eq. 1) then
!                         
!                 p3-------P4
!                 |        |
!                 |        |
!                 |        |
!                 p2------ p1
!
             if( f(xp%x, xp%y) .le. 0.) goto 100
!              check square part 
             p1%z = 0.
             p4%x = p1%x
             p4%y = p1%y
             p4%z = h
          
             p3%x = p2%x
             p3%y = p2%y
             p3%z = h
          
             p2%z = 0.
             call epxpLand4vp(p1, p4, p3, p2, pos, dir, l, kcon, jcon)
             if( kcon <= 4 .and. l > 0 ) then
                length = l
             else
                icon = -1
             endif
          else
!               icon = 0; pos is inside cyl.
             if( f(pos%x, pos%y) .le. 0.) then
                if( f(xp%x, xp%y) .le. 0.) goto 100
             endif
!              check square part 
             p1%z = 0.
             p4%x = p1%x
             p4%y = p1%y
             p4%z = h
          
             p3%x = p2%x
             p3%y = p2%y
             p3%z = h
          
             p2%z = 0.
             if( f(pos%x, pos%y) .le. 0.) then
                call epxpLand4vp(p2, p3, p4, p1, pos, dir, 
     *             l, kcon, jcon)
                if(kcon <= 4 .and. l > 0.0 ) then
                   length = l
                else
!                   icon = -1
!                   write(0,*) ' should not come here'
                endif
             else
                call epxpLand4vp(p1, p4, p3, p2, pos, dir,
     *              l, kcon, jcon)
                if(kcon <= 4 .and. l > 0.0 ) then
                   icon =  1
                   length = l
                else
                   icon = -1
                endif
             endif
          endif
       endif
 100   continue
       end          
!      **********************************
      subroutine epsccyl(comp, pos, icon)
      use modccyl
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


       real*8 r, h, x, y, x1, y1, x2, y2

       real*8 f

       f(x,y) = (y-y1)*(x2-x1) -(y2-y1)*(x-x1)


       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)

       if( pos%z .lt. 0.d0 ) then
          icon = 1
       elseif( pos%z .gt. h ) then
          icon = 1
       elseif(pos%x**2+ pos%y**2 .gt. r**2) then
          icon = 1
       else
          x1 = Volat( comp%vol + ix1)
          y1 = Volat( comp%vol + iy1)
          x2 = Volat( comp%vol + ix2)
          y2 = Volat( comp%vol + iy2)
          if( f(pos%x, pos%y) .le. 0.) then
             icon = 0
          else
             icon = 1
          endif
       endif
       end
!     **************************************
      subroutine epenvlpccyl(comp, org, abc)
      use modccyl
      implicit none
#include "Zglobalc.h"
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

       real*8 r, h, sa, ea
!       logical isinside
       real*8 x, xs, ys, xe, ye
       logical kphiinside !  external function
!       isinside(x) = mod(ea-sa+360.d0, 360.d0) .ge.
!     *               mod(x-sa+360.d0, 360.d0)


       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)
       xs = r*cos(sa*Torad)
       ys = r*sin(sa*Torad)
       xe = r*cos(ea*Torad)
       ye = r*sin(ea*Torad)

!      if(isinside(180.d0)) then
      if(kphiinside(180.d0, sa, ea)) then
          org%x = -r
       else
          org%x = min(xs, xe)
       endif
!       if(isinside(270.d0) )then
       if(kphiinside(270.d0, sa, ea) )then
          org%y = -r
       else
          org%y = min(ys, ye)
       endif
       org%z = 0.d0

       if( kphiinside(0.d0, sa, ea)) then
          abc%x = r - org%x
       else
          abc%x = max(xs, xe) - org%x
       endif

!       if(isinside(90.d0)) then
       if(kphiinside(90.d0, sa, ea)) then
          abc%y = r - org%y
       else
          abc%y = max(ys, ye) - org%y
       endif
       abc%z = h
       NVTX = 0
      end
!     *************************************
      subroutine epatlocccyl(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(*)

      loc(1) = 1
      loc(2) = 2
      loc(3) = 3
      loc(4) = 4
      end
