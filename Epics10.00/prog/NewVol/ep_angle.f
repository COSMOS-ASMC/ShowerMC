!
!  angle
!                                        
!     b    c
!      |****
!      |   |     height is h    
!      |   |                   d <= b & c <=a
!      |   |___________ d      
!      |              | 
!      ---------------| a
!
!   Data format in config is:
!       ox oy oz  a b c d  h
!
!      where (ox,oy,oz) is the origin in the world coord.
!      
      subroutine eprangle(comp)
       implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "angle"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer ia, ib, ic, id, ih, ita, itb, inow
       parameter( ita=1,  itb = 2,  ih=3,   ia=4,
     *   ib=5,  ic=6, id=7, inow=8)

       real*8 a, b, c, d, h
!
!           read angle data as 'new-*'
!           angle has 5 volume attributes and the direction cosines
!           of the  h (1~6)
!
!             next is mandatory
        call eprpst(comp, 5, 8, 1, 6)
!
!           check some values
        a = Volat( comp%vol + 1)
        b = Volat( comp%vol + 2)
        c = Volat( comp%vol + 3)
        d = Volat( comp%vol + 4)
        h = Volat( comp%vol + 5)
        if(a  .le. 0. .or. b .le. 0. .or. h .le. 0. 
     *      .or. c .le. 0. .or. d .le. 0.) then
           write(msg, *) comp%cn, '-th component: ',
     *    '  a, b, c, d, h=', a, b, c, d, h,
     *    ' for angle;  invalid'
           call cerrorMsg(msg, 0)
        endif

        if(c .gt. a .or. d .gt. b) then
           write(msg, *) comp%cn, '-th component: ',
     *    '  c > a or d > b ', ' for angle;  invalid'
           call cerrorMsg(msg, 0)
        endif
        
        Volat( comp%vol + ia) = a
        Volat( comp%vol + ib) = b
        Volat( comp%vol + ic) = c
        Volat( comp%vol + id) = d
        Volat( comp%vol + ih) = h
        Volat( comp%vol + inow) = 0
        call epangleset(comp, 1)
       end
      subroutine epangleset(comp, idx)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "angle"
!
       type(Component)::  comp  ! output. to recieve the config data.
       integer idx

 

       integer ia, ib, ic, id, ih, ita, itb, inow
       parameter( ita=1,  itb = 2,  ih=3,   ia=4,
     *   ib=5,  ic=6, id=7, inow=8)



       if(Volat( comp%vol + inow) .ne. idx) then
          if(idx .eq. 1) then
             Volat( comp%vol + ita) = Volat( comp%vol + ia)
             Volat( comp%vol + itb) = Volat( comp%vol + id)
          else
             Volat( comp%vol + ita) = Volat( comp%vol + ic)
             Volat( comp%vol + itb) = Volat( comp%vol + ib)
          endif
          Volat( comp%vol + inow) = idx
       endif
       end
!   ***************************************
      subroutine epbangle(comp, pos, dir, length, icon)
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


       integer ia, ib, ic, id, ih, ita, itb, inow
       parameter( ita=1,  itb = 2,  ih=3,   ia=4,
     *   ib=5,  ic=6, id=7, inow=8)



       type(epPos):: xp

       real*8 leng2, eps
       real*8 a, b, c, d 
       integer jcon

       data eps/1.d-5/

       a = Volat( comp%vol + ia)
       b = Volat( comp%vol + ib)
       c = Volat( comp%vol + ic)
       d = Volat( comp%vol + id)

       call epangleset(comp, 1)
       call epbbox(comp, pos, dir, length, icon)
       if(icon .eq. 0)  then
          xp%x = pos%x + length*dir%x
          xp%y = pos%y + length*dir%y
          if(xp%x .gt. 0. .and. xp%x .lt. c  .and.
     *         abs(xp%y-d) .lt. eps ) then
             xp%z = pos%z + length*dir%z
             call epangleset(comp, 2)
             goto 100
          else
             goto 200
          endif   
       elseif(icon .eq. 1) then
          call epangleset(comp, 2)
          call epbbox(comp, pos, dir, leng2, jcon)
          if(jcon .eq. 1) then
             length = min(length, leng2)
             goto 200
          elseif(jcon .eq. 0) then
             xp%x = pos%x + leng2*dir%x
             xp%y = pos%y + leng2*dir%y
             if(xp%y .lt. d .and. xp%y .gt. 0. .and.
     *          abs(xp%x -c) .lt. eps) then
                length = leng2
                xp%z = pos%z + length*dir%z
                call epangleset(comp, 1)
                goto 100
             else
                icon = 0
                length  = leng2
                goto 200
             endif
          else
             goto 200
          endif
       else
!           icon = -1 really?
          call epangleset(comp, 2)
          call epbbox(comp, pos, dir, length, icon)
       endif
       goto 200

 100   continue
       call epbbox(comp, xp, dir, leng2, jcon)
       if(jcon .eq. 1) then
          write(0,*) ' xp%x=',xp%x,xp%y,xp%z,
     *      dir%x, dir%y, dir%z, leng2
       elseif(jcon .eq.  0) then
          length = length + leng2
       else
          write(0, *) ' **** x=', xp%x, xp%y, xp%z,
     *   ' pos=', pos%x, pos%y,
     *     pos%z, '  dir=', dir%x, dir%y, dir%z
       endif
       icon = 0
 200   continue
       end          
!      **********************************
      subroutine epsangle(comp, pos, icon)
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



       integer ia, ib, ic, id, ih, ita, itb, inow
       parameter( ita=1,  itb = 2,  ih=3,   ia=4,
     *   ib=5,  ic=6, id=7, inow=8)



       if( pos%z .lt. 0.d0 ) then
          icon = 1
       elseif( pos%z .gt. Volat( comp%vol + ih) ) then
          icon = 1
       elseif(pos%x .lt. 0.) then
          icon = 1
       elseif(pos%y .lt. 0.) then
          icon = 1
       elseif(pos%x .gt. Volat( comp%vol + ia) ) then
          icon = 1
       elseif(pos%y .gt. Volat( comp%vol + ib) )  then
          icon = 1
       elseif(pos%x .gt. Volat( comp%vol + ic) .and.
     *        pos%y .gt. Volat( comp%vol + id)) then
          icon = 1
       else
          icon = 0
       endif
      end
!     **************************************
      subroutine epenvlpangle(comp, org, abc)
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


       integer ia, ib, ic, id, ih, ita, itb, inow
       parameter( ita=1,  itb = 2,  ih=3,   ia=4,
     *   ib=5,  ic=6, id=7, inow=8)

      org%x = 0.
      org%y = 0.
      org%z = 0.
      abc%x = Volat( comp%vol + ia)
      abc%y = Volat( comp%vol + ib)
      abc%z = Volat( comp%vol + ih)

      NVTX = 10  !  concave part neglected

      VTXx(1) = 0.
      VTXy(1) = 0.
      VTXz(1) = 0.

      VTXx(2) = abc%x
      VTXy(2) = 0.
      VTXz(2) = 0.

      VTXx(3) = abc%x
      VTXy(3) = Volat( comp%vol + id)
      VTXz(3) = 0.

      VTXx(4) = Volat( comp%vol + ic)
      VTXy(4) = Volat( comp%vol + ib)
      VTXz(4) = 0.

      VTXx(5) = 0.
      VTXy(5) = Volat( comp%vol + ib)
      VTXz(5) = 0.

      VTXx(6) = 0.
      VTXy(6) = 0.
      VTXz(6) = abc%z

      VTXx(7) = abc%x
      VTXy(7) = 0.
      VTXz(7) = abc%z

      VTXx(8) = abc%x
      VTXy(8) = Volat( comp%vol + id)
      VTXz(8) = abc%z

      VTXx(9) = Volat( comp%vol + ic)
      VTXy(9) = Volat( comp%vol + ib)
      VTXz(9) = abc%z

      VTXx(10) = 0.
      VTXy(10) = Volat( comp%vol + ib)
      VTXz(10) = abc%z

      end
!     *************************************
      subroutine epatlocangle(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(*)
!          bug correction: sep. 09, 2000.
      integer ia, ib, ic, id, ih, ita, itb, inow
      parameter( ita=1,  itb = 2,  ih=3,   ia=4,
     *   ib=5,  ic=6, id=7, inow=8)

      loc(1) = ita  
      loc(2) = itb  
      loc(3) = ic   
      loc(4) = id   
      loc(5) = ih
      end
