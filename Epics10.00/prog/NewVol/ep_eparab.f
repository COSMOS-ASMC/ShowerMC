!
!  eplliptic paraboloid
!                                        
!  ellipse center is at (0,0).  hight is directed to Z.
!      
!     (x/a)**2 + (y/b)**2 = z
!
!   Data format in config is:
!       ox oy oz  a b  h1 h2
!
!      where (ox,oy,oz) is the origin in the world coord.
!           h1: h2:   0 <= h1<= z <= h2
!
!      
      subroutine epreparab(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "eparab"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer ia, ib,  ih1, ih2
       parameter( ia = 1,  ib = 2,  ih1=3, ih2= 4)

       real*8 a, b,  h1,  h2
!
!
!             next is mandatory
        call eprpst(comp, 4, 4, 1, 6)
!
        a = Volat( comp%vol + ia)
        b = Volat( comp%vol + ib)
        h1 = Volat( comp%vol + ih1)
        h2 = Volat( comp%vol + ih2)

        if(a  .le. 0. .or. b .le. 0. .or. h1 .lt. 0. ) then
           write(msg, *) comp%cn, '-th component: a=', a,
     *    ' b=', b, ' h1=',h1, 
     *    ' for eparab;  invalid'
           call cerrorMsg(msg, 0)
        endif
        if(h1 .ge. h2) then
           write(msg, *) comp%cn, '-th component: h1=', h1,
     *    ' >=h2=', h2,
     *    ' for eparab;  invalid'
           call cerrorMsg(msg, 0)
        endif           
       end
!   ***************************************
      subroutine epbeparab(comp, pos, dir, length, icon)
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

 

       integer ia, ib,  ih1, ih2
       parameter( ia = 1,  ib = 2,  ih1=3, ih2= 4)

       real*8 a, b,  h1,  h2


       real*8 aa, bb, cc, dd, f, xx, yy, zz, leng, x, y, z
       real*8 xpa(4)


       integer nx, i

       f(xx, yy, zz) = (xx/a)**2 + (yy/b)**2 - zz

       a = Volat( comp%vol + ia )
       b = Volat( comp%vol + ib )
       h1 = Volat( comp%vol + ih1 )
       h2 = Volat( comp%vol + ih2 )

       aa = (dir%x/a)**2 + (dir%y/b)**2
       bb = 
     *  (pos%x*dir%x/a**2 + pos%y*dir%y/b**2 - dir%z/2.0)
       cc = (pos%x/a)**2 + (pos%y/b)**2 - pos%z

       dd = bb**2 - aa*cc
       
       call epseparab(comp, pos, icon)
       nx = 0

       if(aa .eq. 0. ) then
          if(bb .ne. 0.) then
             leng = -cc/bb/2.0
             if(leng .gt. 0.) then
                z = pos%z + dir%z*leng
                if(z .ge. h1 .and. z .le. h2) then
                   nx = nx + 1
                   xpa(nx)= leng
                endif   
             endif
          endif
       elseif(dd .ge. 0.) then
          dd = sqrt(dd)
          leng = (-bb+dd)/aa
          if(leng .ge. 0.) then
             z = pos%z + dir%z * leng
             if(z .ge. h1 .and. z .le. h2 ) then
                nx = nx + 1
                xpa(nx)= leng
             endif
          endif
          leng = (-bb -dd)/aa
          if( leng .ge. 0.) then
             z = pos%z + dir%z * leng
             if(z .ge. h1 .and. z .le. h2 ) then
                nx = nx + 1
                xpa(nx)= leng
             endif
          endif
       endif
       if(dir%z .ne. 0.) then
          leng = (h2 - pos%z)/dir%z
          if(leng .ge. 0.) then
             x = pos%x + dir%x*leng
             y = pos%y + dir%y*leng
             if( f(x, y, h2) .le. 0.) then
                nx = nx + 1
                xpa(nx) = leng
             endif
          endif
          leng =( h1 - pos%z )/ dir%z
          if(leng .gt. 0.) then
             x = pos%x + leng*dir%x
             y = pos%y + leng*dir%y
             if(f(x, y, h1) .le. 0.)then
                nx = nx +1
                xpa(nx) = leng
             endif
          endif
       endif
       if(nx .eq. 0) then
          icon = -1
       else
          length = xpa(1)
          do i = 2, nx
             length = min(length, xpa(i))
          enddo
       endif
       end          

      subroutine epseparab(comp, pos, icon)
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


       integer ia, ib,  ih1, ih2
       parameter( ia = 1,  ib = 2,  ih1=3, ih2= 4)

       real*8 a, b,  h1,  h2
 

       a = Volat( comp%vol + ia )
       b = Volat( comp%vol + ib )
       h1 = Volat( comp%vol + ih1 )
       h2 = Volat( comp%vol + ih2 )

       if(pos%z .gt. h2) then
          icon = 1
       elseif(pos%z .lt. h1) then
          icon = 1
       elseif(
     *    (pos%x/a)**2 + (pos%y/b)**2 -pos%z
     *      .gt. 0.)    then
          icon = 1
       else
          icon = 0
       endif
      end
!     **************************************
      subroutine epenvlpeparab(comp, org, abc)
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

 

       integer ia, ib,  ih1, ih2
       parameter( ia = 1,  ib = 2,  ih1=3, ih2= 4)

       real*8 a, b,  h1,  h2
 

       a = Volat( comp%vol + ia)
       b = Volat( comp%vol + ib)
       h1 = Volat( comp%vol + ih1)
       h2 = Volat( comp%vol + ih2)
       org%x =-sqrt(h2)*a
       org%y =-sqrt(h2)*b
       org%z = h1
       abc%x = -2*org%x
       abc%y = -2*org%y
       abc%z = h2 - h1
       NVTX = 0
      end
!     *************************************
      subroutine epatloceparab(comp, loc)
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
