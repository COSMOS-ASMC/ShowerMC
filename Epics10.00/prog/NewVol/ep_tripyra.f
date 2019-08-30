!
!        triangular pyramid
!                                  
!
!      |y
!      |
!      |
!      |     (b,c) 
!      |    /\          
!      |   /   \
!      |  /      \
!      | / *(d,e,h)\        * is top vertex
!      _/____________\______x
!             a
!                                                    allowed)
!
!   Data format in config is:
!       ox oy oz  a  b  c, d, e, h  (optional direction cos).
!
!      
      subroutine eprtripyra(comp)
       implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "tripyra"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer ia, ib, ic, id, ie, ih
       parameter( ia = 1,  ib = 2, ic=3,  id=4, ie=5,  ih = 6)

       real*8 a, b, c, d, e, h
!
!           read tripyra data as 'new-*'
!           tripyra has 6 volume attributes and the direction cosines
!           of the 'x' and 'y'==> (1-6)
!
!             next is mandatory
        call eprpst(comp, 6,  6,  1, 6)
!
!           check some values
        a = Volat( comp%vol + ia)
        b = Volat( comp%vol + ib)
        c = Volat( comp%vol + ic)
        d = Volat( comp%vol + id)
        e = Volat( comp%vol + ie)
        h = Volat( comp%vol + ih)

        if(a  .eq. 0. .or. c .eq. 0. .or. h .eq. 0. ) then
           write(msg, *) comp%cn, '-th component: a=', a,
     *    ' c=',c,  ' h=', h,
     *    ' for tripyra;  invalid(must be !=0)'
           call cerrorMsg(msg, 0)
        endif
       end
!   ***************************************
      subroutine epbtripyra(comp, pos, dir, length, icon)
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
 
       integer ia, ib, ic, id, ie, ih
       parameter( ia = 1,  ib = 2, ic=3,  id=4, ie=5,  ih = 6)

       real*8 a, b, c, d, e, h
       integer jcon
!
       
       type(epPos):: p1, p2, p3

       real*8   leng, xpa(2)
       integer np

       a = Volat( comp%vol + ia)
       b = Volat( comp%vol + ib)
       c = Volat( comp%vol + ic)
       d = Volat( comp%vol + id)
       e = Volat( comp%vol + ie)
       h = Volat( comp%vol + ih)
       
       np = 0
!              bottom 
       p1%x = 0.
       p1%y = 0.
       p1%z = 0.

       p2%x = a
       p2%y = 0.
       p2%z = 0.
       
       p3%x = b
       p3%y = c
       p3%z = 0.
       call epxpLand3vp(p1, p2, p3, pos, dir, leng, icon, jcon)
       if(icon .le. 3) then
          np = np +1
          xpa(np) = leng
       endif
!           side x-z
       p1%x = 0.
       p1%y = 0.
       p1%z = 0.

       p2%x = a
       p2%y = 0.
       p2%z = 0.
       
       p3%x = d
       p3%y = e
       p3%z = h
       call epxpLand3vp(p1, p2, p3, pos, dir, leng, icon, jcon)
       if(icon .le. 3) then
          np = np +1
          xpa(np) = leng
       endif
       if(np .eq. 2) goto 100

!           side y-z 1
       p1%x = 0.
       p1%y = 0.
       p1%z = 0.

       p2%x = d
       p2%y = e
       p2%z = h
       
       p3%x = b
       p3%y = c
       p3%z = 0.
       call epxpLand3vp(p1, p2, p3, pos, dir,leng, icon, jcon)
       if(icon .le. 3) then
          np = np +1
          xpa(np) = leng
       endif
       if(np .eq. 2) goto 100
!           y-z 2
       p1%x = a
       p1%y = 0.
       p1%z = 0.

       p2%x = b
       p2%y = c
       p2%z = 0.
       
       p3%x = d
       p3%y = e
       p3%z = h
       call epxpLand3vp(p1, p2, p3, pos, dir,leng, icon, jcon)
       if(icon .le. 3) then
          np = np +1
          xpa(np) = leng
       endif
       if(np .eq. 2) goto 100
       icon = -1
       goto 200
 100   continue
       if(xpa(1) .ge. 0. .and.  xpa(2) .ge. 0.) then
!             outside
          icon = 1
          length = min(xpa(1), xpa(2))
       elseif(xpa(1) .lt. 0. .and.  xpa(2) .lt. 0.) then
!             outside
          icon = 1
          length = max(xpa(1), xpa(2))
       else
!           inside
          icon = 0 
          length = max(xpa(1), xpa(2))
       endif
 200   continue
       end

!      **********************************
      subroutine epstripyra(comp, pos, icon)
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

       integer ia, ib, ic, id, ie, ih
       parameter( ia = 1,  ib = 2, ic=3,  id=4, ie=5,  ih = 6)

       real*8 a, b, c, d, e, h

       type(epPos)::  dir
       real*8 length

      a = Volat( comp%vol + ia)
      b = Volat( comp%vol + ib)
      c = Volat( comp%vol + ic)
      d = Volat( comp%vol + id)
      e = Volat( comp%vol + ie)
      h = Volat( comp%vol + ih)

      if( pos%z .lt. min(0.d0, h) ) then
         icon = 1
      elseif( pos%z .gt. max(h, 0.d0) ) then
         icon = 1
      elseif(pos%x .lt. min(0.d0, a, b, d) ) then
         icon = 1
      elseif(pos%x .gt. max(0.d0, a,  b, d) ) then
         icon = 1
      elseif(pos%y .lt. min(0.d0, c, e)) then
         icon = 1
      elseif(pos%y .gt. max(0.d0, c, e)) then
         icon = 1
      else
!            draw half-line   with dir. (0,0,1)     
         dir%x = 0.
         dir%y = 0.
         dir%z = 1.d0
         call epbtripyra(comp, pos, dir, length, icon)
         if(icon .ne. 0) then
            icon = 1
         endif
      endif

      end
!     **************************************
      subroutine epenvlptripyra(comp, org, abc)
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

       integer ia, ib, ic, id, ie, ih
       parameter( ia = 1,  ib = 2, ic=3,  id=4, ie=5,  ih = 6)

       real*8 a, b, c, d, e, h


      a = Volat( comp%vol + ia)
      b = Volat( comp%vol + ib)
      c = Volat( comp%vol + ic)
      d = Volat( comp%vol + id)
      e = Volat( comp%vol + ie)
      h = Volat( comp%vol + ih)

      org%x = min(0.d0, a, b, d)
      org%y = min(0.d0, c, e)
      org%z = min(0.d0, h)
      abc%x = max(0.d0, a, b, d)- org%x
      abc%y = max(0.d0, c, e) - org%y
      abc%z = h

      NVTX = 4

      VTXx(1) = 0.
      VTXy(1) = 0.
      VTXz(1) = 0.

      VTXx(2) = a
      VTXy(2) = 0.
      VTXz(2) = 0.

      VTXx(3) = b
      VTXy(3) = c
      VTXz(3) = 0.

      VTXx(4) = d
      VTXy(4) = e
      VTXz(4) = h

      end
!     *************************************
      subroutine epatloctripyra(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(6)
 
      integer i

      do i = 1, 6
         loc(i) = i
      enddo
      end

