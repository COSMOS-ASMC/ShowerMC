!
!  cut elliptic cylinder
!                                        
!  bottom elliptic circle center is at (0,0,0).  hight is directed to Z.
!      
!
!
!   Data format in config is:
!       ox oy oz  ra rb  h  sa  ea
!
!      where (ox,oy,oz) is the origin in the world coord.
!            ra: x-radius of the cylinder  cm
!            rb: y-radius //
!            h: height of the //        cm
!           sa: starting angle (deg)  0 is the x-axis. counter clock wise.
!           ea: ending angle (deg).  sa=0 ea=360 means cyl.
!              sa may be > ea.
!      
      subroutine eprcecyl(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "cecyl"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer ira, irb, ih, isa, iea, 
     * ix0, iy0, ix1, iy1
       parameter( ira = 1,  irb=2, ih = 3,  isa=4, iea=5,
     *   ix0=6, iy0=7, ix1=8, iy1=9)

       real*8 ra, rb, h, sa, ea, r
!
        call eprpst(comp, 5, 9, 1, 6)
!
!           check some values
        ra = Volat( comp%vol + ira)
        rb = Volat( comp%vol + irb)
        h = Volat( comp%vol + ih)
        sa= Volat( comp%vol + isa)
        ea= Volat( comp%vol + iea)
        if(ra  .le. 0. .or. h .le. 0. .or. sa .lt. 0. 
     *      .or. ea .gt. 360. .or. rb .le. 0.) then
           write(msg, *) comp%cn, '-th component: r=', r,
     *    ' h=', h, ' sa=',sa, ' ea=',ea,' rb=',rb,
     *    ' for cecyl;  invalid'
           call cerrorMsg(msg, 0)
        endif
!             compute const for later use.
        r = 1.d0/sqrt((cos(sa*Torad)/ra)**2 + (sin(sa*Torad)/rb)**2)
        Volat( comp%vol + ix0) = r*cos(sa*Torad)
        Volat( comp%vol + iy0) = r*sin(sa*Torad)
        r = 1.d0/sqrt((cos(ea*Torad)/ra)**2 + (sin(ea*Torad)/rb)**2)
        Volat( comp%vol + ix1) = r*cos(ea*Torad)
        Volat( comp%vol + iy1) = r*sin(ea*Torad)
       end
!   ***************************************
      subroutine epbcecyl(comp, pos, dir, length, icon)
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


       integer ira, irb, ih, isa, iea, 
     * ix0, iy0, ix1, iy1
       parameter( ira = 1,  irb=2, ih = 3,  isa=4, iea=5,
     *   ix0=6, iy0=7, ix1=8, iy1=9)


       integer kcon, jcon
       type(epPos)::  p1, p2, p3, p4
       real*8 f, h, l, x, y, x1, y1, x2, y2
       type(epPos)::  xp
!              this gives wrong results
!       f(x,y) = (y-p1.y)*(p2.x-p1.x) -(p2.y-p1.y)*(x-p1.x)
       f(x,y) = (y-y1)*(x2 - x1) -(y2-y1)*(x-x1)


       
       call epbecyl(comp, pos, dir, length, icon)
       if( icon .eq. -1 ) then
!               no x-point
       else
          p1%x = Volat( comp%vol + ix0)
          p1%y = Volat( comp%vol + iy0)
          p2%x = Volat( comp%vol + ix1)
          p2%y = Volat( comp%vol + iy1)
          x1 = p1%x
          y1 = p1%y
          x2 = p2%x
          y2 = p2%y
          xp%x = pos%x + length*dir%x
          xp%y = pos%y + length*dir%y
          h = Volat( comp%vol + ih)
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
             if(kcon .le. 4 .and. l> 0.d0) then
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
                if(kcon .le. 4 .and. l > 0.d0 ) then
                   length = l
                else
!                   icon = -1
!                   write(0,*) ' should not come here'
                endif
             else
                call epxpLand4vp(p1, p4, p3, p2, pos, dir,
     *              l, kcon, jcon)
                if(kcon .le. 4 .and. l > 0.d0) then
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
      subroutine epscecyl(comp, pos, icon)
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

       integer ira, irb, ih, isa, iea, 
     * ix0, iy0, ix1, iy1
       parameter( ira = 1,  irb=2, ih = 3,  isa=4, iea=5,
     *   ix0=6, iy0=7, ix1=8, iy1=9)

       real*8 ra, rb, h, sa, ea

 
!

       real*8   ang 

       logical kphiinside
       real*8 x
!       isinside(x) = mod(ea-sa+360.d0, 360.d0) .ge.
!     *               mod(x-sa+360.d0, 360.d0)

       ra = Volat( comp%vol + ira)
       rb = Volat( comp%vol + irb)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)

       if( pos%z .lt. 0.d0 ) then
          icon = 1
       elseif( pos%z .gt. h ) then
          icon = 1
       elseif((pos%x/ra)**2+ (pos%y/rb)**2 .gt. 1.0d0) then
          icon = 1
       else
          ang = atan2(pos%y, pos%x)* Todeg
          if( kphiinside(ang, sa, ea)) then
             icon =0
          else
             icon = 1
          endif
       endif
      end
!     **************************************
      subroutine epenvlpcecyl(comp, org, abc)
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


       integer ira, irb, ih, isa, iea, 
     * ix0, iy0, ix1, iy1
       parameter( ira = 1,  irb=2, ih = 3,  isa=4, iea=5,
     *   ix0=6, iy0=7, ix1=8, iy1=9)

       real*8 ra, rb, h, sa, ea


       logical kphiinside
       real*8 x, xs, ys, xe, ye

       


       ra = Volat( comp%vol + ira)
       rb = Volat( comp%vol + irb)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)
       xs = Volat( comp%vol + ix0)
       ys = Volat( comp%vol + iy0)
       xe = Volat( comp%vol + ix1)
       ye = Volat( comp%vol + iy1)

       if(kphiinside(180.d0, sa, ea)) then
          org%x = - ra
       else
          org%x = min(xs, xe, 0.d0)
       endif
       if(kphiinside(270.d0, sa, ea) )then
          org%y = -rb
       else
          org%y = min(ys, ye, 0.d0)
       endif
       org%z = 0.d0

       if( kphiinside(0.d0, sa, ea)) then
          abc%x = ra - org%x
       else
          abc%x = max(xs, xe, 0.d0) - org%x
       endif

       if(kphiinside(90.d0, sa, ea)) then
          abc%y = rb - org%y
       else
          abc%y = max(ys, ye, 0.d0) - org%y
       endif
       abc%z = h
       NVTX = 0
      end
!     *************************************
      subroutine epatloccecyl(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
       type(Component)::  comp ! input.
      integer loc(*)
      integer i
      
      do i = 1, 5
         loc(i) = i
      enddo
      end


