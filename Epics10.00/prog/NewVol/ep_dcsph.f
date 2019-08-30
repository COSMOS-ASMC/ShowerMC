!
!  duboule cut sphere
!                                        
!   sphere center is at (0,0,0).  hight is directed to Z.
!
!
!                    
!    h2         *---------*
!            *      |         *
!    h1     ------------------
!                   |          / 
!                   |        /
!                   |      /
!                   |    / r
!                   |  /
!                    / 
!   Data format in config is:
!       ox oy oz  r  h1 h2   (optional dir)
!
!      where (ox,oy,oz) is the origin in the world coord.
!            r: radius of the sphere
!            h1: height of the //     
!            h2: //
!               -r<=h1<= h2<=r
      subroutine eprdcsph(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "dcsph"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer ir, ih1, ih2, irh1, irh2
       parameter (ir = 1,  ih1 = 2, irh1=3, ih2=4, irh2=5 )

       real*8 r, h1, h2
!
!           read cut cylinder data as 'new-*'
!           ccyl has 2 volume attributes and the direction cosines
!           of the  h (7~9)
!
!             next is mandatory
        call eprpst(comp, 3, 5,  7, 9)
!
!           next is optional
!           check some values
        r = Volat( comp%vol + ir)
        h1 = Volat( comp%vol + 2)
        h2 = Volat( comp%vol + 3)
        Volat( comp%vol + ih1) = h1
        Volat( comp%vol + ih2) = h2
        if(r  .le. 0. .or. abs(h1) .gt. r .or. abs(h2) .gt. r) then
           write(msg, *) comp%cn, '-th component: r=', r,
     *    ' h1,2=', h1, h2, ' for dcsph;  invalid'
           call cerrorMsg(msg, 0)
        endif
        if(h1 .gt. h2) then
           write(msg, *) comp%cn, '-th component: h1=',h1,
     *     ' > h2=', h2 
           call cerrorMsg(msg,0)
        endif
        Volat( comp%vol + irh1) = sqrt( r**2 - h1**2)
        Volat( comp%vol + irh2) = sqrt( r**2 - h2**2)
       end
!   ***************************************
      subroutine epbdcsph(comp, pos, dir, length, icon)
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


       integer ir, ih1, ih2, irh1, irh2
       parameter (ir = 1,  ih1 = 2, irh1=3, ih2=4, irh2=5 )

       real*8 r, h1, h2

       type(epPos)::  xp

       r = Volat( comp%vol + ir)
       h1 = Volat( comp%vol + ih1)
       h2 = Volat( comp%vol + ih2)
       
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
          if(xp%z .ge. h1 .and. xp%z .le. h2 ) goto 100
!              get x,y at h
          if(dir%z .gt. 0. .and. pos%z .lt. h1) then
             xp%z = h1
             length = (h1 - pos%z)/dir%z
             xp%x = pos%x + length *dir%x
             xp%y = pos%y + length *dir%y
             if(xp%x**2 + xp%y**2 .le.
     *           Volat( comp%vol + irh1)**2) goto 100
          endif
          if(dir%z .lt. 0. .and. pos%z .gt. h2) then
             xp%z = h2
             length = (h2 - pos%z)/dir%z
             xp%x = pos%x + length *dir%x
             xp%y = pos%y + length *dir%y
             if(xp%x**2 + xp%y**2 .le.
     *        Volat( comp%vol + irh2)**2) goto 100
          endif
          icon = -1
       elseif(icon .eq. 0) then
          if(pos%z .ge. h1 .and. pos%z .le. h2) then
             xp%z = pos%z + length*dir%z 
             if(xp%z .ge.  h1 .and. xp%z .le. h2) goto 100
          endif
          if(dir%z .gt. 0. .and. pos%z .lt. h1) then
             xp%z = h1
             length = (h1 - pos%z)/dir%z
             xp%x = pos%x + length *dir%x
             xp%y = pos%y + length *dir%y
             if(xp%x**2 + xp%y**2 .le. 
     *      Volat( comp%vol + irh1)**2 ) then
                icon = 1
             else
                icon = -1 
             endif
          elseif(dir%z .lt. 0. .and. pos%z .gt. h2) then
             xp%z = h2
             length = (h2 - pos%z)/dir%z
             xp%x = pos%x + length *dir%x
             xp%y = pos%y + length *dir%y
             if(xp%x**2 + xp%y**2 .le. 
     *         Volat( comp%vol + irh2)**2 ) then
                icon = 1
             else
                icon = -1 
             endif
          elseif(dir%z .ne. 0. .and. pos%z .gt. h1 .and.
     *            pos%z .lt. h2 ) then
             if(dir%z .gt. 0.) then
                xp%z = h2
                length = (h2 - pos%z)/dir%z
                xp%x = pos%x + length *dir%x
                xp%y = pos%y + length *dir%y
                if(xp%x**2 + xp%y**2 .le.
     *             Volat( comp%vol + irh2)**2 ) then
                   icon = 0
                else
                   write(0,*) ' strange 1'
                   stop 134
                endif
             elseif( dir%z .lt. 0.) then    
                xp%z = h1
                length = (h1 - pos%z)/dir%z
                xp%x = pos%x + length *dir%x
                xp%y = pos%y + length *dir%y
                if(xp%x**2 + xp%y**2 .le.
     *              Volat( comp%vol + irh1)**2 ) then
                   icon = 0
                else
                   write(0,*) ' strange 2'
                   stop 222
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
      subroutine epsdcsph(comp, pos, icon)
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

 
       integer ir, ih1, ih2, irh1, irh2
       parameter (ir = 1,  ih1 = 2, irh1=3, ih2=4, irh2=5 )

       real*8 r, h1, h2



       r = Volat( comp%vol + ir)
       h1 = Volat( comp%vol + ih1)
       h2 = Volat( comp%vol + ih2)

       if( pos%z .lt. h1 ) then
          icon = 1
       elseif( pos%z .gt. h2 ) then
          icon = 1
       elseif(pos%x**2+ pos%y**2 + pos%z**2 .gt. r**2) then
          icon = 1
       else
          icon = 0
       endif
       end
!     **************************************
      subroutine epenvlpdcsph(comp, org, abc)
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

 
       integer ir, ih1, ih2, irh1, irh2
       parameter (ir = 1,  ih1 = 2, irh1=3, ih2=4, irh2=5 )

       real*8 r, h1, h2, rh1, rh2



       r = Volat( comp%vol + ir)
       h1 = Volat( comp%vol + ih1)
       h2 = Volat( comp%vol + ih2)
       rh1 = Volat( comp%vol + irh1)
       rh2 = Volat( comp%vol + irh2)

       if(h1 .gt. 0.) then
          org%x = -rh1
       elseif(h2 .lt. 0.) then
          org%x = -rh2
       else 
          org%x = -r
       endif
       org%y = org%x
       org%z = h1
       
       abc%x = -2.0 * org%x
       abc%y = abc%x
       abc%z = h2 - org%z
       NVTX = 0 
      end
!     *************************************
      subroutine epatlocdcsph(comp, loc)
!      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer  loc(3)
 


       integer ir, ih1, ih2, irh1, irh2
       parameter (ir = 1,  ih1 = 2, irh1=3, ih2=4, irh2=5 )

       loc(1) = ir
       loc(2) = ih1
       loc(3) = ih2
       end

