c
c     see Fig/NewVol2.pdf 
c      
      subroutine eproctagon(comp)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
c
c         interface to read configuration data for "octagon"
c
      record /Component/ comp   ! output. to recieve the config data.
      real(8),save:: a, b, c, d 
c           read octagon data as 'new-*'
c           octagon has 4 volume attributes and the direction cosines
c           of the 'x' and 'y'==> (1-6)
c
c             next is mandatory
      call eprpst(comp, 4, 4, 1, 6)
c
c           check some values
      a = Volat( comp.vol + 1)
      b = Volat( comp.vol + 2)
      c = Volat( comp.vol + 3) 
      d = Volat( comp.vol + 4) 
      if( a<=0 .or. b<=0 .or. c <=0 .or. d<=0 ) then
         write(0,*) ' a, b, c or d of octagon are invalid'
         write(0,*) a, b, c, d
         stop
      endif
      if( a < 2*d ) then
         write(0,*) ' a< 2d for octagon;  a',a, ' d=',d
         stop
      endif
      if( c < 2*d ) then
         write(0,*) ' c< 2d for octagon;  c',c, ' d=',d
         stop
      endif
      end
      subroutine epboctagon(comp, posl, dirl, el, icon)
      implicit none
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"

c
c        find length to the boundary of 'comp' from 'pos'
c        with direction cos 'dir'
c     'pos' and 'dir' are given in this 'comp' local coordinate.
c 
 

      record /Component/comp    ! input. you can extract volume parameters
                          !            by Volat( comp.vol + 1), etc
      record /epPos/ posl       ! input.  position.
      record /epDirec/ dirl      ! input. direction cosinse

      real(8):: el                !  output length cm from pos to the boundary
      integer:: icon              ! output 0: el obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume
 
      integer base, jcon
      real(8):: x, y, z
      real(8),parameter:: eps=1.d-7
      real(8):: a, b, c, d, ela(2) 
      integer:: nc, i

      base = comp.vol
      a = Volat( base + 1)
      b = Volat( base + 2)
      c = Volat( base + 3) 
      d = Volat( base + 4) 
!            find x-ing point with the box
      call kxplbx(
     *       posl.x, posl.y, posl.z, dirl.x, dirl.y, dirl.z,
     *       a, b, c,  el, jcon)
      if(jcon == -1) then
!             no x-ing.
         icon =  jcon
         return
      endif
!           x-ing point
      x = posl.x + el*dirl.x
      y = posl.y + el*dirl.y
      z = posl.z + el*dirl.z
      if( abs(x) <= eps .or. abs(x-a)<=eps ) then
          if( y  >= 0. .and. y <= b ) then
             if( z >= d .and.  z <= c-d ) then
!                 judge posl's in/out         
                call epsoctagon(comp, posl, icon)
                return
             endif
          endif
       elseif( abs(y) <= eps .or. abs(y-b) <= eps ) then
          if( x >= 0.  .and. x <= a ) then
             if( z >= 0.  .and. z <= c ) then
               if (z >= (d-x)  .and. z >= x-(a-d) .and.
     *             z <= x + (c-d) .and. z <= -x +a+c-d ) then
!                 judge posl's in/out         
                  call epsoctagon(comp, posl, icon)
                  return
               endif
            endif
         endif
      elseif( abs(z) <= eps .or. abs(z-c) <= eps ) then
         if( x >= d .and. x <= (a-d) ) then
            if(y >= 0. .and. y <= b) then
!              judge posl's in/out         
               call epsoctagon(comp, posl, icon)
               return
            endif
         endif
      endif
c        see Xsing point with inclined planes
c        left bottom
c            z=-x +d;   x + 0y + z = d
      nc = 0   ! clear counter for x-ing point
      call kxplp(posl.x, posl.y, posl.z, dirl.x, dirl.y, dirl.z,
     *   1.d0, 0.d0, 1.0d0,  d,     el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         x = posl.x + el *dirl.x
         if( x >= 0. .and. x <= d )  then
            y = posl.y + el *dirl.y
            if(y >=0. .and. y<= b ) then
               nc = nc + 1
               ela(nc) = el
            endif
         endif
      endif
c        right bottom
c            z=x -(a-d);   x + 0y - z = a-d
      call kxplp(posl.x, posl.y, posl.z, dirl.x, dirl.y, dirl.z,
     *   1.d0, 0.d0, -1.0d0, a-d,     el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         x = posl.x + el *dirl.x
         if( x >= a-d .and. x <= a )  then
            y = posl.y + el *dirl.y
            if(y >=0. .and. y<= b ) then
               nc = nc + 1
               ela(nc) = el
            endif
         endif
      endif
c        right top
c            z=-x + c + (a-d);   x + 0y + z = a-d +c
      call kxplp(posl.x, posl.y, posl.z, dirl.x, dirl.y, dirl.z,
     *   1.d0, 0.d0, 1.0d0, a-d+c,   el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         x = posl.x + el *dirl.x
         if( x >= a-d .and. x <= a )  then
            y = posl.y + el *dirl.y
            if(y >=0. .and. y<= b ) then
               nc = nc + 1
               ela(nc) = el
            endif
         endif
      endif
c        left top
c            z=x + (c-d);   x + 0y - z = d-c
      call kxplp(posl.x, posl.y, posl.z, dirl.x, dirl.y, dirl.z,
     *   1.d0, 0.d0, -1.0d0, d-c,     el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         x = posl.x + el *dirl.x
         if( x >=0.  .and. x <= d )  then
            y = posl.y + el *dirl.y
            if(y >=0. .and. y<= b ) then
               nc = nc + 1
               ela(nc) = el
            endif
         endif
      endif
      if( nc == 0) then
         icon = -1
      elseif(nc > 2 ) then
         write(0,*) ' error in epboctagon; x-point > 2 '
         stop
      else
         ! find minium distance
         el = 1.d20
         do i = 1, nc
            if( ela(i) < el ) then
               el = ela(i)
            endif
         enddo
!              judge posl's in/out         
         call epsoctagon(comp, posl, icon)
         
      endif
      end

c      **********************************
      subroutine epsoctagon(comp, pos, icon)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
c
c           judge if a given 'pos' is inside 'comp'
c         
      record /Component/ comp !input component
      record /epPos/ pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside
      real(8),save:: a, b, c, d 
      
      integer:: base
      base = comp.vol
      a = Volat( base + 1)
      b = Volat( base + 2)
      c = Volat( base + 3) 
      d = Volat( base + 4) 

      if( pos.z < 0.d0 .or. pos.z > c ) then
         icon = 1
         return
      elseif(pos.x < 0.d0 .or. pos.x > a) then
         icon = 1
         return
      elseif(pos.y < 0.d0 .or. pos.y > b) then
         icon = 1
         return
      else
!          left bot   \      z=-x + d
         if( pos.z < -pos.x + d ) then
            icon = 1
!           right bot   /   z = x -(a-d)
         elseif ( pos.z < pos.x  - (a-d) ) then
            icon = 1
!             right top \  z = -x + (a-d)+c
         elseif( pos.z > -pos.x + a - d  +c ) then
            icon = 1
!             left top   /  z= x +c-d
         elseif( pos.z > pos.x + c -d ) then
            icon = 1
         else
            icon = 0
         endif
      endif
      end
c     **************************************
      subroutine epenvlpoctagon(comp, org, abc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

c
c        give the envloping box of the component
c
      record /Component/ comp  ! input.   component.
      record /epPos/ org       ! output.  origin of the enveloping box
                               !          in local coord. 
      record /ep3Vec/ abc      ! output.  a,b,c of the box

      integer:: base
      base = comp.vol

      abc.x = Volat( base + 1)
      abc.y = Volat( base + 2)
      abc.z = Volat( base + 3) 

      NVTX = 0

      end
c    *************************************
      subroutine epatlococtagon(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

      record /Component/ comp ! input.
      integer loc(4)
 
      integer i

      do i = 1, 4
         loc(i) = i
      enddo
      end


