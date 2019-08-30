!   polygon.  lying on x-y plane at z=0  verteces upto Max given
!   in ZepMaxdef.h
!        thickness is h
!      
      module polygon
      integer,save::Compnum=-1000
      integer,save:: npoly  ! # polygon with npoly verteces. (
                      ! 1, 2, ... npoly.  the last vertex is connected
                      ! to the first one.
      real(8),save:: a, b, c, d    !  octagon consts. order is allways same as canonical one
    end module polygon

      subroutine eproctagon(comp)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
c
c         interface to read configuration data for "polygon"
c
      record /Component/ comp   ! output. to recieve the config data.
c           read polygon data as 'new-*'
c           A polygon has npoly+1  volume attributes and spcification
!          of the number of vertices(npoly)  and the direction cosines
c           of the 'x' and 'y'==> (1-6)
c
!            polygon Al 0 0 0 / x y z npoly h  v1 v2 ... vnpoly dirx diry
c
c             next is mandatory
      call eprpst(comp, 4, 4, 1, 6)
      call epoctagonCnst(comp)
      if( comp.struc == 'octagon') then
      elseif( comp.struc == 'octagon_x' ) then
      elseif(index(confdata, 'octagon_y') > 0) then
      elseif(index(confdata, 'octagon_z') > 0) then
      else
         write(0,*) 'structure=',comp.struc, ' not usable'
         stop
      endif
      end

      subroutine epoctagonCnst(comp)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
c
c         interface to read configuration data for "octagon"
c
      record /Component/ comp   ! output. to recieve the config data.

c           check some values
c      if( Compnum /= comp.cn ) then  ! see honeycomb
      if( Compnum /= comp.vol ) then
         a = Volat( comp.vol + 1)
         b = Volat( comp.vol + 2)
         c = Volat( comp.vol + 3) 
         d = Volat( comp.vol + 4) 
c         Compnum = comp.cn
         Compnum = comp.vol
         if( a<=0 .or. b<=0 .or. c <=0 .or. d<=0 ) then
            write(0,*) ' a, b, c or d of octagon are invalid'
            write(0,*) a, b, c, d
            stop
         endif
         if( b < 2*d ) then
            write(0,*) ' b< 2d for octagon;  b',b, ' d=',d
            stop
         endif
         if( c < 2*d ) then
            write(0,*) ' c< 2d for octagon;  c',c, ' d=',d
            stop
         endif
      endif
      end

      subroutine epboctagon(comp, posl, dirl, el, icon)
      use octagon
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

      record /epPos/ cposl
      record /epDirec/ cdirl

!      call epoctagonCnst(comp)
      call epv2c_octagon(comp, posl, cposl)
      call epv2cd_octagon(comp, dirl, cdirl)

      call epboctagn0(comp, cposl, cdirl, el, icon)
      end

c      **********************************
      subroutine epsoctagon(comp, pos, icon)
      use octagon
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
      record /epPos/ cpos

!      call epoctagonCnst(comp)
      call epv2c_octagon(comp, pos, cpos)
      call epsoctagn0(comp, cpos, icon)
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

      real(8):: temp
      call  epenvlpoctagn0(comp, org, abc)
      if(comp.struc == 'octagon' ) then
      elseif( comp.struc == 'octagon_x' ) then
      elseif( comp.struc == 'octagon_y' ) then
         !  x<->y 
         temp = abc.x
         abc = ep3Vec(abc.y, temp, abc.z)
      elseif( comp.struc == 'octagon_z' ) then
         ! x<->z exchange
         temp = abc.x
         abc = ep3Vec(abc.z, abc.y, temp)
      else
         write(0,*) ' comp.struc=',comp.struc
         write(0,*) ' error'
         stop
      endif
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



!*****************************
      subroutine epboctagn0(comp, posl, dirl, el, icon)
      use octagon
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
 
      integer jcon
      real(8):: x, y, z
      real(8):: ela(2)
      real(8),parameter:: eps=1.d-8

      integer:: nc, i
      call epoctagonCnst(comp)
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
      if( abs(y) <= eps .or. abs(y-b)<=eps ) then
          if( x  >= 0. .and. x <= a ) then
             if( z >= d .and.  z <= c-d ) then
!                 judge posl's in/out         
                call epsoctagn0(comp, posl, icon)
                return
             endif
          endif
       elseif( abs(x) <= eps .or. abs(x-a) <= eps ) then
          if( y >= 0.  .and. y <= b ) then
             if( z >= 0.  .and. z <= c ) then
               if (z >= (d-y)  .and. z >= y-(b-d) .and.
     *             z <= y + (c-d) .and. z <= -y +b+c-d ) then
!                 judge posl's in/out         
                  call epsoctagn0(comp, posl, icon)
                  return
               endif
            endif
         endif
      elseif( abs(z) <= eps .or. abs(z-c) <= eps ) then
         if( y >= d .and. y <= (b-d) ) then
            if(x >= 0. .and. x <= a) then
!              judge posl's in/out         
               call epsoctagn0(comp, posl, icon)
               return
            endif
         endif
      endif
c        see Xsing point with inclined planes
c        left bottom
c            z=-y +d;   0x + y + z = d
      nc = 0   ! clear counter for x-ing point
      call kxplp(posl.x, posl.y, posl.z, dirl.x, dirl.y, dirl.z,
     *   0.d0, 1.d0, 1.0d0,  d,     el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         y = posl.y + el *dirl.y
         if( y >= 0. .and. y <= d )  then
            x = posl.x + el *dirl.x
            if(x >=0. .and. x<= a ) then
               nc = nc + 1
               ela(nc) = el
            endif
         endif
      endif
c        right bottom
c            z=y -(b-d);   0x + y - z = b-d
      call kxplp(posl.x, posl.y, posl.z, dirl.x, dirl.y, dirl.z,
     *   0.d0, 1.d0, -1.0d0, b-d,     el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         y = posl.y + el *dirl.y
         if( y >= b-d .and. y <= b )  then
            x = posl.x + el *dirl.x
            if(x >=0. .and. x<= a ) then
               nc = nc + 1
               ela(nc) = el
            endif
         endif
      endif
c        right top
c            z=-y + c + (b-d);   0x + y + z = b-d +c
      call kxplp(posl.x, posl.y, posl.z, dirl.x, dirl.y, dirl.z,
     *   0.d0, 1.d0, 1.0d0, b-d+c,   el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         y = posl.y + el *dirl.y
         if( y >= b-d .and. y <= b )  then
            x = posl.x + el *dirl.x
            if(x >=0. .and. x<= a ) then
               nc = nc + 1
               ela(nc) = el
            endif
         endif
      endif
c        left top
c            z=y + (c-d);   0x + y - z = d-c
      call kxplp(posl.x, posl.y, posl.z, dirl.x, dirl.y, dirl.z,
     *   0.d0, 1.d0, -1.0d0, d-c,     el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         y = posl.y + el *dirl.y
         if( y >=0.  .and. y <= d )  then
            x = posl.x + el *dirl.x
            if(x >=0. .and. x<= a ) then
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
         call epsoctagn0(comp, posl, icon)
         
      endif
      end

c      **********************************
      subroutine epsoctagn0(comp, pos, icon)
      use octagon
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

      
      call epoctagonCnst(comp)


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
!          left bot   \      z=-y + d
         if( pos.z < -pos.y + d ) then
            icon = 1
!           right bot   /   z = y -(b-d)
         elseif ( pos.z < pos.y  - (b-d) ) then
            icon = 1
!             right top \  z = -y + (b-d)+c
         elseif( pos.z > -pos.y + b - d  +c ) then
            icon = 1
!             left top   /  z= y +c-d
         elseif( pos.z > pos.y + c -d ) then
            icon = 1
         else
            icon = 0
         endif
      endif
      end
c     **************************************
      subroutine epenvlpoctagn0(comp, org, abc)
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

      org =  epPos(0.d0, 0.d0, 0.d0)

      abc.x = Volat( base + 1)
      abc.y = Volat( base + 2)
      abc.z = Volat( base + 3) 

      NVTX = 0
      end

      subroutine epv2c_octagon(comp, posv, posc)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
      record /Component/ comp  ! input.   component.
      record /epPos/ posv     ! input
      record /epPos/ posc     ! output

      call epoctagonCnst(comp)
      if( comp.struc  == 'octagon' ) then
         posc = posv
      elseif( comp.struc  == 'octagon_x' ) then
         posc = posv
      elseif( comp.struc  == 'octagon_y' ) then
         posc = epPos( posv.y, b-posv.x, posv.z)
      elseif( comp.struc  == 'octagon_z' ) then
         posc = epPos( posv.z, b-posv.y, posv.x)
      else
         write(0,*) ' comp.struc=', comp.struc,
     *      ' to epboctagon is invalid'
         stop
      endif
      end   subroutine epv2c_octagon

      subroutine epc2v_octagon(comp, posc, posv)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
      record /Component/ comp  ! input.   component.
      record /epPos/ posc     ! input
      record /epPos/ posv     ! output

      call epoctagonCnst(comp)
!           (X,Y,Z)--->(x,y,z)
      if( comp.struc  == 'octagon' ) then
         posv = posc
      elseif( comp.struc  == 'octagon_x' ) then
         posv = posc
      elseif( comp.struc  == 'octagon_y' ) then
         posv = epPos( b-posc.y, posc.x,  posc.z)
      elseif( comp.struc  == 'octagon_z' ) then
         posv = epPos( posc.z, b-posc.y, posc.x)
      else
         write(0,*) ' comp.struc=', comp.struc,
     *      ' to epboctagon is invalid'
         stop
      endif
      end   subroutine epc2v_octagon
      subroutine epv2cd_octagon(comp, dirv, dirc)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
      record /Component/ comp  ! input.   component.
      record /epDirec/ dirv     ! input
      record /epDirec/ dirc     ! output

      call epoctagonCnst(comp)

      if( comp.struc  == 'octagon' ) then
         dirc = dirv
      elseif( comp.struc  == 'octagon_x' ) then
         dirc = dirv
      elseif( comp.struc  == 'octagon_y' ) then
         dirc = epDirec( dirv.y, -dirv.x, dirv.z)
      elseif( comp.struc  == 'octagon_z' ) then
         dirc = epDirec( dirv.z, -dirv.y, dirv.x)
      else
         write(0,*) ' comp.struc=', comp.struc,
     *      ' to epboctagon is invalid'
         stop
      endif
      end

      subroutine epc2vd_octagon(comp, dirc, dirv)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
      record /Component/ comp  ! input.   component.
      record /epDirec/ dirc     ! input
      record /epDirec/ dirv     ! output

      call epoctagonCnst(comp)
!           (X,Y,Z)<---(x,y,z)
      if( comp.struc  == 'octagon' ) then
         dirv = dirc
      elseif( comp.struc  == 'octagon_x' ) then
         dirv = dirc
      elseif( comp.struc  == 'octagon_y' ) then
         dirv = epDirec(-dirc.y, dirc.x,  dirc.z)
      elseif( comp.struc  == 'octagon_z' ) then
         dirv = epDirec( dirc.z,-dirc.y, dirc.z)
      else
         write(0,*) ' comp.struc=', comp.struc,
     *      ' to epboctagon is invalid'
         stop
      endif
      end   subroutine epc2vd_octagon
