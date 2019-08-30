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
!
!   Data format in config is:
!       ox oy oz  a  b  h  x0 y0 a' b'  x.x x.y x.z  y.x y.y y.z
!
!      where (ox,oy,oz) is the origin in the world coord.
!            (a b h x0 y0 a' b') are the parameters to describe the horse
!            (x.x ..  y.z) optional. direction cosines of the
!      'a'-axis and 'b'-axis in the world coordinate. If not
!       given, (1,0,0) and (0,1,0) are assumed.
!      
      subroutine eprhorse(comp)
       implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "horse"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer ia, ib, ih, ix0, iy0,  iap, ibp
       parameter( ia = 1,  ib = 2,  ih = 3,  ix0=4, iy0=5,
     *            iap=6, ibp=7 )

       real*8 a, b, h, x0, y0,  ap, bp
!
!           read horse data as 'new-*'
!           horse has 7 volume attributes and the direction cosines
!           of the 'a' and 'b'==> (1-6)
!
!             next is mandatory
        call eprpst(comp, 7, 7, 1, 6)
!
!           next is optional
!           check some values
        a = Volat( comp%vol + ia)
        b = Volat( comp%vol + ib)
        h = Volat( comp%vol + ih)
        x0= Volat( comp%vol + ix0)
        y0= Volat( comp%vol + iy0)
        ap = Volat( comp%vol + iap)
        bp = Volat( comp%vol + ibp)
        if(a  .lt. 0. .or. b .lt. 0. .or. h .le. 0. 
     *      .or. ap .lt. 0. .or. bp .lt. 0.) then
           write(msg, *) comp%cn, '-th component: a=', a,
     *    ' b=', b, ' h=', h, ' x0=',x0, ' y0=',y0,
     *    " a'=",ap, " b'=",bp,
     *    ' for horse;  invalid(must>=0)'
           call cerrorMsg(msg, 0)
        endif
        if(a .eq. 0. .and. b .eq. 0. and. ap .eq. 0. and. 
     *     bp .eq. 0.) then
           write(msg, *)
     *     comp%cn, "-th component: a=b=a'=b'=0 for horse"
           call cerrorMsg(msg, 0)
        endif
        
       end
!   ***************************************
      subroutine epbhorse(comp, pos, dir, length, icon)
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
 

       type(Component):: comp  ! input. you can extract volume parameters
                               !  by Volat( comp.vol + 1), etc
       type(epPos)::  pos   ! input.  position.
       type(epDirec)::  dir  ! input. direction cosinse

       real*8  length !  output length cm from pos to the boundary
       integer icon  ! output 0: length obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume
 
       integer ia, ib, ih, ix0, iy0,  iap, ibp
       parameter( ia = 1,  ib = 2,  ih = 3,  ix0=4, iy0=5,
     *            iap=6, ibp=7 )

       real*8 a, b, h, x0, y0,  ap, bp

       integer np, jcon
!
       
       type(epPos):: p1, p2, p3, p4

       real*8   l, xpa(2)
       integer:: face 

!           check some values
       a = Volat( comp%vol + ia)
       b = Volat( comp%vol + ib)
       h = Volat( comp%vol + ih)
       x0= Volat( comp%vol + ix0)
       y0= Volat( comp%vol + iy0)
       ap = Volat( comp%vol + iap)
       bp = Volat( comp%vol + ibp)

!////       call kxplhorse(a, b, h, x0, y0, ap, bp,  length, icon, face)
       call kxplhorse(a, b, h, x0, y0, ap, bp, pos, dir,
     *  length, icon, face)
       end

!      **********************************
      subroutine epshorse(comp, pos, icon)
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

       integer ia, ib, ih, ix0, iy0,  iap, ibp
       parameter( ia = 1,  ib = 2,  ih = 3,  ix0=4, iy0=5,
     *            iap=6, ibp=7 )

       real*8 a, b, h, x0, y0,  ap, bp

       type(epPos)::  dir
       real*8 length

      a = Volat( comp%vol + ia)
      b = Volat( comp%vol + ib)
      h = Volat( comp%vol + ih)
      x0= Volat( comp%vol + ix0)
      y0= Volat( comp%vol + iy0)
      ap = Volat( comp%vol + iap)
      bp = Volat( comp%vol + ibp)

      if( pos%z .lt. 0.d0 ) then
         icon = 1
      elseif( pos%z .gt. h ) then
         icon = 1
      elseif(pos%x .lt. min(0.d0, x0) ) then
         icon = 1
      elseif(pos%x .gt.  max(a, x0+ap) ) then
         icon = 1
      elseif(pos%y .lt. min(0.d0, y0)) then
         icon = 1
      elseif(pos%y .gt. max(b, y0+bp)) then
         icon = 1
      else
!            draw half-line   with dir. (0,0,1)     

         dir%x = 0.
         dir%y = 0.
         dir%z = 1.d0
         call epbhorse(comp, pos, dir, length, icon)
         if(icon .ne. 0) then
            icon = 1
         endif
      endif
      end
!     **************************************
      subroutine epenvlphorse(comp, org, abc)
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

       integer ia, ib, ih, ix0, iy0,  iap, ibp
       parameter( ia = 1,  ib = 2,  ih = 3,  ix0=4, iy0=5,
     *            iap=6, ibp=7 )

       real*8 a, b, h, x0, y0,  ap, bp


      a = Volat( comp%vol + ia)
      b = Volat( comp%vol + ib)
      h = Volat( comp%vol + ih)
      x0= Volat( comp%vol + ix0)
      y0= Volat( comp%vol + iy0)
      ap = Volat( comp%vol + iap)
      bp = Volat( comp%vol + ibp)

 

      org%x = min(0.d0, x0)
      org%y = min(0.d0, y0)
      org%z = 0.d0
      abc%x = max(a, x0+ap)-org%x  ! -org.x   v9.157
      abc%y = max(b, y0+bp)-org%y  ! -org.y    //
      abc%z = h


      NVTX = 0   !///////////   v9.157 following commented
!      VTXx(1) = 0.     
!      VTXy(1) = 0.
!      VTXz(1) = 0.
!
!      VTXx(2) = a
!      VTXy(2) = 0.
!      VTXz(2) = 0.
!
!      VTXx(3) = a
!      VTXy(3) = b
!      VTXz(3) = 0.
!
!      VTXx(4) = 0.
!      VTXy(4) = b
!      VTXz(4) = 0.
!
!      VTXx(5) = x0
!      VTXy(5) = y0
!      VTXz(5) = h
!
!      VTXx(6) = x0 + ap
!      VTXy(6) = y0
!      VTXz(6) = h
!
!      VTXx(7) = x0 + ap
!      VTXy(7) = y0 + bp
!      VTXz(7) = h
!
!      VTXx(8) = x0 
!      VTXy(8) = y0 + bp
!      VTXz(8) = h

      end
!     *************************************
      subroutine epatlochorse(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(7)
 
      integer i

      do i = 1, 7
         loc(i) = i
      enddo
      end
