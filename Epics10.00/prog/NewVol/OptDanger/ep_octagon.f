!         Octagon
!     see Fig/NewVol2.pdf 
!      
      module octagon
      integer,save::Compnum=-1000
      real(8),save:: a, b, c, d    !  octagon consts. order is allways same as canonical one
      end module octagon

      subroutine eproctagon(comp)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
!         interface to read configuration data for "octagon"
!
       type(Component)::  comp   ! output. to recieve the config data.
!           read octagon data as 'new-*'
!           octagon has 4 volume attributes and the direction cosines
!           of the 'x' and 'y'==> (1-6)
!
      character(len=MAX_STRUCCHR):: basename 
!             next is mandatory
      call eprpst(comp, 4, 4, 1, 6)
      call epoctagonCnst(comp)

      call epGetBaseStrucName(comp%struc, basename)

      if( basename /= 'octagon') then
         write(0,*) 'structure=',comp%struc, ' not usable'
         stop
      endif
      end

      subroutine epoctagonCnst(comp)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "octagon"
!
       type(Component)::  comp   ! output. to recieve the config data.

!           check some values
      if( Compnum /= comp%cn ) then
         a = Volat( comp%vol + 1)
         b = Volat( comp%vol + 2)
         c = Volat( comp%vol + 3) 
         d = Volat( comp%vol + 4) 
         Compnum = comp%cn
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

!
!        find length to the boundary of 'comp' from 'pos'
!        with direction cos 'dir'
!     'pos' and 'dir' are given in this 'comp' local coordinate.
! 
 

       type(Component):: comp    ! input. you can extract volume parameters
                          !            by Volat( comp.vol + 1), etc
       type(epPos)::  posl       ! input.  position.
       type(epDirec)::  dirl      ! input. direction cosinse

      real(8):: el                !  output length cm from pos to the boundary
      integer:: icon              ! output 0: el obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume

       type(epPos)::  cposl
       type(epDirec)::  cdirl

!      call epoctagonCnst(comp)
      call epv2c_octagon(comp, posl, cposl)
      call epv2cd_octagon(comp, dirl, cdirl)

      call epboctagn0(comp, cposl, cdirl, el, icon)
      end

!      **********************************
      subroutine epsoctagon(comp, pos, icon)
      use octagon
      implicit none
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
       type(epPos)::  cpos

!      call epoctagonCnst(comp)
      call epv2c_octagon(comp, pos, cpos)
      call epsoctagn0(comp, cpos, icon)
      end
!     **************************************
      subroutine epenvlpoctagon(comp, org, abc)
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

      real(8):: temp
      call  epenvlpoctagn0(comp, org, abc)
      if(comp%struc == 'octagon' .or.
     *   comp%struc == 'octagon_w' ) then
      elseif( comp%struc(1:9) == 'octagon_x' ) then
      elseif( comp%struc(1:9) == 'octagon_y' ) then
         !  x<->y 
         temp = abc%x
         abc = ep3Vec(abc%y, temp, abc%z)
      elseif( comp%struc(1:9) == 'octagon_z' ) then
         ! x<->z exchange
         temp = abc%x
         abc = ep3Vec(abc%z, abc%y, temp)
      else
         write(0,*) ' comp%struc=',comp%struc
         write(0,*) ' error'
         stop
      endif
      end
!    *************************************
      subroutine epatlococtagon(comp, loc)
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

!
!        find length to the boundary of 'comp' from 'pos'
!        with direction cos 'dir'
!     'pos' and 'dir' are given in this 'comp' local coordinate.
! 
 

       type(Component):: comp    ! input. you can extract volume parameters
                          !            by Volat( comp.vol + 1), etc
       type(epPos)::  posl       ! input.  position.
       type(epDirec)::  dirl      ! input. direction cosinse

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
     *       posl%x, posl%y, posl%z, dirl%x, dirl%y, dirl%z,
     *       a, b, c,  el, jcon)
      if(jcon == -1) then
!             no x-ing.
         icon =  jcon
         return
      endif
!           x-ing point
      x = posl%x + el*dirl%x
      y = posl%y + el*dirl%y
      z = posl%z + el*dirl%z
      if( abs(y) <= eps .or. abs(y-b)<=eps ) then
          if(x  >= -eps .and. x <= a+eps ) then
             if( z >= d-eps .and.  z <= c-d+eps ) then
!                 judge posl's in/out         
                call epsoctagn0(comp, posl, icon)
                return
             endif
          endif
       elseif( abs(x) <= eps .or. abs(x-a) <= eps ) then
          if( y >= -eps  .and. y <= b+eps ) then
             if( z >= -eps  .and. z <= c+eps ) then
               if (z >= (d-y)-eps .and. z >= y-(b-d)-eps .and.
     *             z <= y + (c-d)+eps .and. z <= -y +b+c-d+eps ) then
!                 judge posl's in/out         
                  call epsoctagn0(comp, posl, icon)
                  return
               endif
            endif
         endif
      elseif( abs(z) <= eps .or. abs(z-c) <= eps ) then
         if( y >= d-eps .and. y <= (b-d)+eps ) then
            if(x >= -eps .and. x <= a+eps) then
!              judge posl's in/out         
               call epsoctagn0(comp, posl, icon)
               return
            endif
         endif
      endif
!        see Xsing point with inclined planes
!        left bottom
!            z=-y +d;   0x + y + z = d
      nc = 0   ! clear counter for x-ing point
      call kxplp(posl%x, posl%y, posl%z, dirl%x, dirl%y, dirl%z,
     *   0.d0, 1.d0, 1.0d0,  d,     el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         y = posl%y + el *dirl%y
         if( y >= -eps .and. y <= d +eps)  then
            x = posl%x + el *dirl%x
            if(x >=-eps .and. x<= a+eps ) then
               nc = nc + 1
               ela(nc) = el
            endif
         endif
      endif
!        right bottom
!            z=y -(b-d);   0x + y - z = b-d
      call kxplp(posl%x, posl%y, posl%z, dirl%x, dirl%y, dirl%z,
     *   0.d0, 1.d0, -1.0d0, b-d,     el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         y = posl%y + el *dirl%y
         if( y >= b-d-eps .and. y <= b+eps )  then
            x = posl%x + el *dirl%x
            if(x >=-eps .and. x<= a+eps ) then
               nc = nc + 1
               ela(nc) = el
            endif
         endif
      endif
!        right top
!            z=-y + c + (b-d);   0x + y + z = b-d +c
      call kxplp(posl%x, posl%y, posl%z, dirl%x, dirl%y, dirl%z,
     *   0.d0, 1.d0, 1.0d0, b-d+c,   el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         y = posl%y + el *dirl%y
         if( y >= b-d-eps .and. y <= b+eps )  then
            x = posl%x + el *dirl%x
            if(x >=-eps .and. x<= a+eps ) then
               nc = nc + 1
               ela(nc) = el
            endif
         endif
      endif
!        left top
!            z=y + (c-d);   0x + y - z = d-c
      call kxplp(posl%x, posl%y, posl%z, dirl%x, dirl%y, dirl%z,
     *   0.d0, 1.d0, -1.0d0, d-c,     el, jcon)
      if(jcon == 0 .and. el > 0. ) then
         y = posl%y + el *dirl%y
         if( y >=-eps  .and. y <= d+eps )  then
            x = posl%x + el *dirl%x
            if(x >=-eps .and. x<= a+eps ) then
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

!      **********************************
      subroutine epsoctagn0(comp, pos, icon)
      use octagon
      implicit none
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

      
      call epoctagonCnst(comp)


      if( pos%z < 0.d0 .or. pos%z > c ) then
         icon = 1
         return
      elseif(pos%x < 0.d0 .or. pos%x > a) then
         icon = 1
         return
      elseif(pos%y < 0.d0 .or. pos%y > b) then
         icon = 1
         return
      else
!          left bot   \      z=-y + d
         if( pos%z < -pos%y + d ) then
            icon = 1
!           right bot   /   z = y -(b-d)
         elseif ( pos%z < pos%y  - (b-d) ) then
            icon = 1
!             right top \  z = -y + (b-d)+c
         elseif( pos%z > -pos%y + b - d  +c ) then
            icon = 1
!             left top   /  z= y +c-d
         elseif( pos%z > pos%y + c -d ) then
            icon = 1
         else
            icon = 0
         endif
      endif
      end
!     **************************************
      subroutine epenvlpoctagn0(comp, org, abc)
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

      integer:: base
      base = comp%vol

      org =  epPos(0.d0, 0.d0, 0.d0)

      abc%x = Volat( base + 1)
      abc%y = Volat( base + 2)
      abc%z = Volat( base + 3) 

      NVTX = 0
      end

      subroutine epv2c_octagon(comp, posv, posc)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp  ! input.   component.
       type(epPos)::  posv     ! input
       type(epPos)::  posc     ! output

      call epoctagonCnst(comp)
      if( comp%struc  == 'octagon' .or.
     *    comp%struc  == 'octagon_w' ) then
         posc = posv
      elseif( comp%struc(1:9)  == 'octagon_x' ) then
         posc = posv
      elseif( comp%struc(1:9)  == 'octagon_y' ) then
         posc = epPos( posv%y, b-posv%x, posv%z)
      elseif( comp%struc(1:9)  == 'octagon_z' ) then
         posc = epPos( posv%z, b-posv%y, posv%x)
      else
         write(0,*) ' comp%struc=', comp%struc,
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
       type(Component)::  comp  ! input.   component.
       type(epPos)::  posc     ! input
       type(epPos)::  posv     ! output

      call epoctagonCnst(comp)
!           (X,Y,Z)--->(x,y,z)
      if( comp%struc  == 'octagon' .or.
     *    comp%struc  == 'octagon_w'  ) then
         posv = posc
      elseif( comp%struc(1:9)  == 'octagon_x' ) then
         posv = posc
      elseif( comp%struc(1:9)  == 'octagon_y' ) then
         posv = epPos( b-posc%y, posc%x,  posc%z)
      elseif( comp%struc(1:9)  == 'octagon_z' ) then
         posv = epPos( posc%z, b-posc%y, posc%x)
      else
         write(0,*) ' comp%struc=', comp%struc,
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
       type(Component)::  comp  ! input.   component.
       type(epDirec)::  dirv     ! input
       type(epDirec)::  dirc     ! output

      call epoctagonCnst(comp)

      if( comp%struc  == 'octagon' .or.
     *    comp%struc  == 'octagon_w' ) then
         dirc = dirv
      elseif( comp%struc(1:9)  == 'octagon_x' ) then
         dirc = dirv
      elseif( comp%struc(1:9)  == 'octagon_y' ) then
         dirc = epDirec( dirv%y, -dirv%x, dirv%z)
      elseif( comp%struc(1:9)  == 'octagon_z' ) then
         dirc = epDirec( dirv%z, -dirv%y, dirv%x)
      else
         write(0,*) ' comp%struc=', comp%struc,
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
       type(Component)::  comp  ! input.   component.
       type(epDirec)::  dirc     ! input
       type(epDirec)::  dirv     ! output

      call epoctagonCnst(comp)
!           (X,Y,Z)<---(x,y,z)
      if( comp%struc  == 'octagon' .or.
     *    comp%struc  == 'octagon_w' ) then
         dirv = dirc
      elseif( comp%struc(1:9)  == 'octagon_x' ) then
         dirv = dirc
      elseif( comp%struc(1:9)  == 'octagon_y' ) then
         dirv = epDirec(-dirc%y, dirc%x,  dirc%z)
      elseif( comp%struc(1:9)  == 'octagon_z' ) then
         dirv = epDirec( dirc%z,-dirc%y, dirc%z)
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epboctagon is invalid'
         stop
      endif
      end   subroutine epc2vd_octagon
