!
! square pipe
!                                        
!  
!      
!  
!
!   Data format in config is:
!       ox oy oz  a b c  x0 y0 a' b'
!
!      where (ox,oy,oz) is the origin in the world coord.
!         a, b, c,  a',b' > 0   
!         x0, y0 >=0
!         x0+a' <= a
!         y0+b' <= b  
!
!      
      subroutine eprsqpipe(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "sqpipe"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*170 msg
 
       integer ia, ib, ig, ix0, iy0, iap, ibp
       parameter( ia = 1,  ib = 2,  ig=3, ix0=4, iy0= 5, iap=6, ibp=7)

       real*8 a, b, g, x0, y0, ap, bp
!
!             next is mandatory
        call eprpst(comp, 7, 7, 1, 6)
!
        a = Volat( comp%vol + ia)
        b = Volat( comp%vol + ib)
        g = Volat( comp%vol + ig)
        x0 = Volat( comp%vol + ix0)
        y0 = Volat( comp%vol + iy0)
        ap = Volat( comp%vol + iap)
        bp = Volat( comp%vol + ibp)

!         check
!         a, b, c,  a',b' > 0   
!         x0, y0 >=0

        if(a  .le. 0. .or. b .le. 0. .or. g .le. 0. 
     *    .or. ap .le.0. .or. bp .le. 0. ) then
           write(msg, *) comp%cn, '-th component: a=', a,
     *    ' b=', b, ' c=',g, " a'=",ap, " b'=",bp,
     *    ' for sqpipe;  invalid(must be > 0)'
           call cerrorMsg(msg, 0)
        endif
        if(x0 .lt. 0. .or. y0 .lt. 0.) then
           write(msg, *) comp%cn, '-th component: x0=', x0,
     *    ' y0=', y0, ' for sqpipe;  invalid (must be >=0)'
           call cerrorMsg(msg, 0)
        endif
!           check
!         x0+a' <= a
!         y0+b' <= b  
!
        if(x0+ap .gt. a .or. y0+bp  .gt. b ) then
           write(msg, *) comp%cn,"-th component: x0+a'<=a",
     *     " or y0+b'<=b is violated for sqpipe"
           call cerrorMsg(msg, 0)
        endif           
       end
!   ***************************************
      subroutine epbsqpipe(comp, pos, dir, length, icon)
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

 

       integer ia, ib, ig, ix0, iy0, iap, ibp
       parameter( ia = 1,  ib = 2,  ig=3, ix0=4, iy0= 5, iap=6, ibp=7)

       real*8 a, b, g, x0, y0, ap, bp
!
       type(epPos)::  p1, p2, p3, p4, q1, q2, q3, q4
       type(epPos)::  r1, r2, r3, r4, s1, s2, s3, s4

       real*8 xpa(4), x, y, leng

       integer nx, i, jcon

       a = Volat( comp%vol + ia )
       b = Volat( comp%vol + ib )
       g = Volat( comp%vol + ig )
       x0 = Volat( comp%vol + ix0 )
       y0 = Volat( comp%vol + iy0 )
       ap = Volat( comp%vol + iap )
       bp = Volat( comp%vol + ibp )

       nx = 0

       p1%x = 0.
       p1%y = 0.
       p1%z = 0.

       p2%x = a
       p2%y = 0.
       p2%z = 0.
       
       p3%x = a
       p3%y = b
       p3%z = 0.
       
       p4%x = 0.
       p4%y = b
       p4%z = 0.


       q1%x = 0.
       q1%y = 0.
       q1%z = g

       q2%x = a
       q2%y = 0.
       q2%z = g
       
       q3%x = a
       q3%y = b
       q3%z = g
       
       q4%x = 0.
       q4%y = b
       q4%z = g

       r1%x = x0
       r1%y = y0
       r1%z = 0.
       
       r2%x = x0+ap
       r2%y = y0
       r2%z = 0.
       
       r3%x = x0+ap
       r3%y = y0+bp
       r3%z = 0.
       
       r4%x = x0
       r4%y = y0+bp
       r4%z = 0.

       s1%x = x0
       s1%y = y0
       s1%z = g

       s2%x = x0+ap
       s2%y = y0
       s2%z = g
       
       s3%x = x0+ap
       s3%y = y0+bp
       s3%z = g
       
       s4%x = x0
       s4%y = y0+bp
       s4%z = g
!          outer  x-z: y=0
       call epxpLand4vp(
     *    p1, p2, q2, q1, pos, dir, leng, icon, jcon)
       if(icon .le. 4 .and. leng > 0.d0) then
          if(y0 .gt. 0.) then
             nx = nx + 1
             xpa(nx) = leng
          else
             x = pos%x + leng*dir%x
             if(x .le. x0 .or. x .ge. x0+ap) then
                nx = nx + 1
                xpa(nx) = leng
             endif
          endif
       endif
!           inner x-z; y=y0
       if(y0 .gt. 0.) then
          call epxpLand4vp(
     *      r1, r2, s2, s1, pos, dir, leng, icon, jcon)
          if(icon .le. 4 .and. leng >  0.d0) then
             nx = nx + 1
             xpa(nx) = leng
          endif
       endif
!            outer  x-z: y=b
       call epxpLand4vp(
     *      p3, q3, q4, p4, pos, dir, leng, icon, jcon)
       if(icon .le. 4 .and. leng >  0.d0) then
          if(y0+bp .lt. b) then
             nx = nx + 1
             xpa(nx) = leng
          else
             x = pos%x + leng*dir%x
             if(x .le. x0 .or. x .ge. x0+ap) then
                nx = nx + 1
                xpa(nx) = leng
             endif
          endif
       endif
!           inner x-z; y=y0+b'
       if(y0+bp .lt. b) then
          call epxpLand4vp(
     *      r3, s3, s4, r4, pos, dir, leng, icon, jcon)
          if(icon .le. 4 .and. leng > 0.d0) then
             nx = nx + 1
             xpa(nx) = leng
             if(nx .eq. 4) goto 100
          endif
       endif

!          outer  y-z: x=0

       call epxpLand4vp(
     *    p1, q1, q4, p4, pos, dir, leng, icon, jcon)
       if(icon .le. 4 .and. leng > 0.d0) then
          if(x0 .gt. 0.) then
             nx = nx + 1
             xpa(nx) = leng
             if(nx .eq. 4) goto 100
          else
             y = pos%y + leng*dir%y
             if(y .le. y0 .or. y .ge. y0+bp) then
                nx = nx + 1
                xpa(nx) = leng
                if(nx .eq. 4) goto 100
             endif
          endif
       endif
!          inner y-z; x=x0
       if(x0 .ne. 0.) then
          call epxpLand4vp(
     *    r1, s1, s4, r4, pos, dir, leng, icon, jcon)
          if(icon .le. 4 .and. leng > 0.d0) then
             nx = nx + 1
             xpa(nx) = leng
             if(nx .eq. 4) goto 100
          endif
       endif

!          outer  y-z: x=a
       call epxpLand4vp(
     *    p2, p3, q3, q2, pos, dir, leng, icon, jcon)
       if(icon .le. 4 .and. leng >  0.d0) then
          if(x0+ap .lt. a) then
             nx = nx + 1
             xpa(nx) = leng
             if(nx .eq. 4) goto 100
          else
             y = pos%y + leng * dir%y
             if(y .le. y0 .or. y .ge. y0+bp) then
                nx = nx + 1
                xpa(nx) = leng
                if(nx .eq. 4) goto 100
             endif
          endif
       endif
!          inner y-z; x=x0+a'
       if( x0 + ap .lt. a) then
          call epxpLand4vp(
     *    r3, s3, s2, r2, pos, dir, leng, icon, jcon)
          if(icon .le. 4 .and. leng >  0.d0) then
             nx = nx + 1
             xpa(nx) = leng
             if(nx .eq. 4) goto 100
          endif
       endif
!        -------------
!         top        
       call epxpLand4vp(
     *    q1, q2, q3, q4, pos, dir, leng, icon, jcon)
       if(icon .le. 4 .and. leng > 0.d0) then
          x = pos%x + leng*dir%x
          y = pos%y + leng*dir%y
          if( x .le. x0 .or. x .ge. x0+ap .or.
     *        y .le. y0 .or. y .ge. y0+bp) then
             nx = nx + 1
             xpa(nx) = leng
             if(nx .eq. 4) goto 100
          endif
       endif
!         bottom
       call epxpLand4vp(
     *    p1, p2, p3, p4, pos, dir, leng, icon, jcon)
       if(icon .le. 4 .and. leng > 0.d0) then
          x = pos%x + leng*dir%x
          y = pos%y + leng*dir%y
          if( x .le. x0 .or. x .ge. x0+ap .or.
     *        y .le. y0 .or. y .ge. y0+bp) then
             nx = nx + 1
             xpa(nx) = leng
             if(nx .eq. 4) goto 100
          endif
       endif
 100   continue
       if(nx .eq. 0) then
          icon = -1
       else
          length = xpa(1)
          do i = 2, nx
             length = min(length, xpa(i))
          enddo
          call epssqpipe(comp, pos, icon)
       endif
       end          

      subroutine epssqpipe(comp, pos, icon)
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

       integer ia, ib, ig, ix0, iy0, iap, ibp
       parameter( ia = 1,  ib = 2,  ig=3, ix0=4, iy0= 5, iap=6, ibp=7)

       real*8 a, b, g, x0, y0, ap, bp
!

       a = Volat( comp%vol + ia )
       b = Volat( comp%vol + ib )
       g = Volat( comp%vol + ig )
       x0 = Volat( comp%vol + ix0 )
       y0 = Volat( comp%vol + iy0 )
       ap = Volat( comp%vol + iap )
       bp = Volat( comp%vol + ibp )
 

       if( pos%z .gt. g) then
          icon = 1
       elseif( pos%z .lt. 0.) then
          icon = 1
       elseif( pos%x .gt. a) then
          icon = 1
       elseif( pos%x .lt. 0.) then
          icon = 1
       elseif( pos%y .gt. b) then
          icon = 1
       elseif(pos%y .lt. 0.) then
          icon = 1
       elseif(pos%x .gt. x0 .and. pos%x .lt. x0+ap .and.
     *        pos%y .gt. y0 .and. pos%y .lt. y0+bp) then
          icon = 1
       else
          icon = 0
       endif
      end
!     **************************************
      subroutine epenvlpsqpipe(comp, org, abc)
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


      integer ia, ib, ig, ix0, iy0, iap, ibp
      parameter( ia = 1,  ib = 2,  ig=3, ix0=4, iy0= 5, iap=6, ibp=7)


!

       org%x =0.
       org%y =0.
       org%z =0.
       abc%x = Volat( comp%vol + ia )
       abc%y = Volat( comp%vol + ib )
       abc%z = Volat( comp%vol + ig )
       NVTX =0
      end
!     *************************************
      subroutine epatlocsqpipe(comp, loc)
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
