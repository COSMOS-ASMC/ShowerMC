!
!  cut sphere
!                                        
!   sphere center is at (0,0,0).  hight is directed to Z.
!
!
!                    _
!               *         *
!            *               *
!           -------------------
!                   |          / 
!                   |        /
!                 h |      /
!                   |    / r
!                   |  /
!                    / 
!   Data format in config is:
!       ox oy oz  r  h   (optional dir)
!
!      where (ox,oy,oz) is the origin in the world coord.
!            r: radius of the sphere
!            h: height of the //        cm (may be < 0)
!      
      subroutine eprcsph(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "csph"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer ir, ih, irh
       parameter (ir = 1,  ih = 2, irh=3)

       real*8 r, h
!
!           read cut cylinder data as 'new-*'
!           ccyl has 2 volume attributes and the direction cosines
!           of the  h (7~9)
!
!             next is mandatory
        call eprpst(comp, 2, 3, 7, 9)
!
!           next is optional
!           check some values
        r = Volat( comp%vol + ir)
        h = Volat( comp%vol + ih)
        if(r  .le. 0. .or. abs(h) .gt. r) then
           write(msg, *) comp%cn, '-th component: r=', r,
     *    ' h=', h, ' for csph;  invalid'
           call cerrorMsg(msg, 0)
        endif
        Volat( comp%vol + irh) = sqrt( r**2 - h**2)
       end
!   ***************************************
      subroutine epbcsph(comp, pos, dir, length, icon)
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

       integer ir, ih, irh
       parameter (ir = 1,  ih = 2, irh=3)

       real*8  r, h 


       type(epPos)::  xp

       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)

       
       call kxplsph(pos%x, pos%y, pos%z, dir%x, dir%y, dir%z, r, 
     *    length, icon)
!  output:
!     length: x-ssing point is at pos + length*dir ( el>=0)
!   icon : output. 0 the point is in the sphere. length is obtained
!                  1 the point is out side of the sphere. length is
!                     obtained.
!                 -1 no x-ing point
       if( icon .eq. -1 ) then
!               no x-point
       elseif(icon .eq. 1) then
          xp%z = pos%z + dir%z*length
          if(xp%z .gt. h) goto 100
!              get x,y at h
          if(dir%z .gt. 0. .and. pos%z .lt. h) then
             xp%z = h
             length = (h - pos%z)/dir%z
             xp%x = pos%x + length *dir%x
             xp%y = pos%y + length *dir%y
             if(xp%x**2 + xp%y**2 .le.
     *                Volat( comp%vol + irh)**2) goto 100
          endif
          icon = -1
       elseif(icon .eq. 0) then
          if(pos%z .gt. h) then
             if(pos%z + length*dir%z .gt. h) goto 100
          endif
          if(dir%z .ne. 0.) then
             xp%z = h
             length = (h - pos%z)/dir%z
             xp%x = pos%x + length *dir%x
             xp%y = pos%y + length *dir%y
             if(xp%x**2 + xp%y**2 .le. 
     *            Volat( comp%vol + irh)**2 ) then
                if(pos%z .lt. h ) then
                   icon = 1
                endif
             else
                icon = -1
             endif
          else
             icon = -1
          endif
       endif
 100   continue
       end          
!      **********************************
      subroutine epscsph(comp, pos, icon)
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

 
       integer ir, ih, irh

       parameter (ir = 1,  ih = 2,  irh=3)

       real*8 r, h

       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)

       if( pos%z .lt. h ) then
          icon = 1
       elseif( pos%z .gt. r ) then
          icon = 1
       elseif(pos%x**2+ pos%y**2 + pos%z**2 .gt. r**2) then
          icon = 1
       else
          icon = 0
       endif
       end
!     **************************************
      subroutine epenvlpcsph(comp, org, abc)
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

 
       integer ir, ih, irh

       parameter (ir = 1,  ih = 2,  irh = 3)


       real*8 r, h, rh

       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
       rh = Volat( comp%vol + irh)

       if(h .gt. 0.) then
          org%x = -rh
       else
          org%x = -r
       endif
       if(h .gt. 0.) then
          org%y = -rh
       else
          org%y = -r
       endif
       org%z = h
       
       abc%x = -2.0 * org%x
       abc%y = abc%x
       abc%z = r - org%z
       NVTX = 0
      end
!     *************************************
      subroutine epatloccsph(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(2)
 
      loc(1) = 1
      loc(2) = 2
      end

