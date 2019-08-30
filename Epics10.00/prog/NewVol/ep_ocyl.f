!
!  open cylinder
!                                        
!  bottom circle center is at (0,0,0).  hight is directed to Z.
!      
!
!
!   Data format in config is:
!       ox oy oz  r  h  sa  ea
!
!      where (ox,oy,oz) is the origin in the world coord.
!            r: radius of the cylinder  cm
!            h: height of the //        cm
!           sa: starting angle (deg)  0 is the x-axis. counter clock wise.
!           ea: ending angle (deg).  sa=0 ea=360 means cyl.
!              sa may be > ea.
!      
      module modocyl
      implicit none
      integer,parameter::ir = 1,  ih = 2,  isa=3, iea=4,
     *   ircossa=5, irsinsa=6, ircosea=7, irsinea=8
      integer,save:: where
        ! 1 top,  3; cut part sa; 4; cut part ea; 6 bottom
      end  module modocyl

      subroutine eprocyl(comp)
      use modocyl
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "ocyl"
!
       type(Component)::  comp   ! output. to recieve the config data.
      character*150 msg
      real*8 r, h, sa, ea
!
!           read open cylinder data as 'new-*'
!           ocyl has 4 volume attributes and the direction cosines
!           of the  h (1~6)
!
!     next is mandatory
      call eprpst(comp, 4, 8, 1, 6)
!
!           next is optional
!           check some values
      r = Volat( comp%vol + ir)
      h = Volat( comp%vol + ih)
      sa= Volat( comp%vol + isa)
      ea= Volat( comp%vol + iea)
      if(r  .le. 0. .or. h .le. 0. .or. sa .lt. 0. 
     *     .or. ea .gt. 360.) then
         write(msg, *) comp%cn, '-th component: r=', r,
     *        ' h=', h, ' sa=',sa, ' ea=',ea,
     *        ' for ocyl;  invalid'
         call cerrorMsg(msg, 0)
      endif
!             compute const for later use.
      Volat( comp%vol + ircossa) = r*cos(sa*Torad)
      Volat( comp%vol + irsinsa) = r*sin(sa*Torad)
      Volat( comp%vol + ircosea) = r*cos(ea*Torad)
      Volat( comp%vol + irsinea) = r*sin(ea*Torad)
      end

      subroutine epbocyl(comp, pos, dir, length, icon)
      use modocyl
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

       real*8 r, h, sa, ea

       integer  jcon, kcon
!
       type(epPos)::  p1, p2, p3, p4, xp

       real*8   ltemp,  angs,  angx
       real*8 x
       integer:: wtemp
       integer:: insideflag
!          wrong function
!       isinside(x) = mod(ea-sa+360.d0, 360.d0) .ge.
!     *               mod(x-sa+360.d0, 360.d0)

       logical kphiinside  ! external function


       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)

       where = 0

       call kxplcy(pos%x, pos%y, pos%z, dir%x, dir%y, dir%z, r, h,
     *  length,    icon, where)
!  output:
!     length: x-sing point is at pos + length*dir ( el>=0)
!   icon : output. 0 the point is in the cyl. length is obtained
!                  1 the point is outside of the cyl. length is
!                     obtained.
!                 -1 no x-ing point
!   where: output. 1  x-ing point is on x-y  top plane.
!                  2  //             on the side.
!                  6  //             on      bottom.
!                 -1  no x-ing point
       
       if( sa .eq. 0. and. ea .eq. 360.d0 ) then
!           full cyl. 
          goto 100
       elseif( icon .eq. -1 ) then
!               no x-point
          goto 100
       else
          call epsocyl(comp, pos, insideflag)  ! pos is in/out ?
          call  ep_ocylSeeSqWall(comp, pos, dir,   ! xp two walls
     *            ltemp,   wtemp, jcon )
          xp%x = pos%x + length*dir%x
          xp%y = pos%y + length*dir%y
          angx = atan2(xp%y, xp%x)*Todeg
          if( kphiinside(angx, sa, ea) ) then
              ! direction to xp  is in sa~ea
             if( jcon == 0 ) then  ! so if there is xp with wall
                if( ltemp < length ) then  ! nearer one is real xp
                   length = ltemp
                   where = wtemp
                   icon = insideflag
                   goto 100
                else  ! use length, where as they are now
                   icon = insideflag
                   goto 100
                endif
             else ! use length, where as they are now
                icon = insideflag
             endif
          else  ! candidate xp  is outside direcion
             if( jcon == 0 ) then  ! so take wall xp 
                length = ltemp
                where = wtemp
                icon = insideflag
             else
                if( insideflag == 0 ) then 
                   write(0,*) 'ocyl: point is inside but no xp'
                   write(0,*) 
     *              'so strange;jcon,icon= ',jcon, icon
                   write(0,*) ' angx=',angx
                endif
                icon = -1
             endif
          endif
       endif
 100   continue
       end          
      subroutine  ep_ocylSeeSqWall(comp, pos, dir, 
     *    length,   wtemp, icon )
!       see the  crossing points with the two walls in ocyl
!       and get nearer one if any.
      use modocyl
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
       real(8),intent(out)::length !  length cm from pos to the boundary, if obtained

       integer,intent(out):: wtemp  ! xp's face #
       integer,intent(out):: icon   ! 0: xp  found. -1 not found
       

       integer:: kcon,jcon  ! output from  epxpLand4vp
       real(8):: ltemp
       integer:: nxp

       real(8)::p1(3), p2(3), p3(3), p4(3)

       nxp = 0

       p1(1) = Volat( comp%vol + ircossa)
       p1(2) = Volat( comp%vol + irsinsa)
       p1(3) = 0.
             
       p2(1:2) = p1(1:2)
       p2(3) = Volat( comp%vol + ih)

       p3(1:2) = 0.d0
       p3(3) = Volat( comp%vol + ih)
          
       p4(:) =  0.d0

       call epxpLand4vp(p1, p2, p3, p4, pos, dir, 
     *  ltemp, kcon, jcon)
       if(kcon <=  4 .and. ltemp >= 0. ) then
          nxp = nxp + 1
          length = ltemp
          wtemp = 3
          icon = 0
       endif
!               x.p  may be on the other cut part
       p1(1) = Volat( comp%vol + ircosea)
       p1(2) = Volat( comp%vol + irsinea)
       p1(3) = 0 

       p2(:) = 0.d0
          
       p3(1:2) = 0.d0
       p3(3) = Volat( comp%vol + ih)
          
       p4(1:2) = p1(1:2)
       p4(3) = Volat( comp%vol + ih)
       call epxpLand4vp(p1, p2, p3, p4, pos, dir, 
     *  ltemp, kcon, jcon)
       if(kcon <=  4 .and. ltemp >= 0.) then
          if( nxp > 0 ) then
             if( ltemp < length ) then
                length = ltemp
                wtemp = 4
                icon = 0
             else  ! use the one already obtained
             endif
          else  ! use the new one
             length = ltemp
             wtemp = 4
             icon = 0
          endif
       elseif( nxp == 0 ) then
          icon = -1
       else  ! use the one already obtained
       endif

       end subroutine  ep_ocylSeeSqWall

       subroutine epqocyl(whichpart)
       use modocyl
       implicit none
!         inquire the pos. of x.p 
       integer whichpart ! output. 1; top, 
                        !         2: side (full cyl)
                        !         3: cut sa:
                        !         4: cut ea;
                        !         6: bottom  

       
       whichpart = where
       end

!      **********************************
      subroutine epsocyl(comp, pos, icon)
      use modocyl
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

 

       real*8 r, h, sa, ea
!

       real*8   ang 

       logical kphiinside
       real*8 x
       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)

       if( pos%z .lt. 0.d0 ) then
          icon = 1
       elseif( pos%z .gt. h ) then
          icon = 1
       elseif(pos%x**2+ pos%y**2 .gt. r**2) then
          icon = 1
       else
          ang = atan2(pos%y, pos%x)* Todeg
          if(kphiinside(ang, sa, ea)) then
             icon =0
          else
             icon = 1
          endif
!//////////////
!          write(0,*) 'in check: pos=',pos.x, pos.y, pos.z
!          write(0,*) ' ang =',ang, ' kphi=', kphiinside(ang, sa, ea)
!          write(0,*) ' icon =',icon
!//////////////
       endif

      end
!     **************************************
      subroutine epenvlpocyl(comp, org, abc)
      use modocyl
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

 

       real*8 r, h, sa, ea

       real*8 x, xs, ys, xe, ye
       logical kphiinside





       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)
       xs = Volat( comp%vol + ircossa)
       ys = Volat( comp%vol + irsinsa)
       xe = Volat( comp%vol + ircosea)
       ye = Volat( comp%vol + irsinea)

       if( kphiinside(180.d0, sa, ea  )) then
          org%x = -r
       else
          org%x = min(xs, xe, 0.d0)
       endif
       if(kphiinside(270.d0, sa, ea) )then
          org%y = -r
       else
          org%y = min(ys, ye, 0.d0)
       endif
       org%z = 0.d0

       if( kphiinside(0.d0, sa,ea)) then
          abc%x = r - org%x
       else
          abc%x = max(xs, xe, 0.d0) - org%x
       endif

       if( kphiinside(90.d0, sa, ea)) then
          abc%y = r - org%y
       else
          abc%y = max(ys, ye, 0.d0) - org%y
       endif
       abc%z = h
       NVTX = 0
      end
!     *************************************
      subroutine epatlococyl(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(4)
 
      integer i

      do i = 1, 4
         loc(i) = i
      enddo
      end

