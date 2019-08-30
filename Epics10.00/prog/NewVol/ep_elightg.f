!  elightg: elliptic light guide
!                                        
!  see Doc/Fig/elightg.sd
!
!   Data format in config is:
!       ox oy oz  a b t1 t2 w1 w1p w2 w2p [dir]
!
!      where (ox,oy,oz) is the origin in the world coord.
!         a>0, b>0, t1>=0, t2>=0 (but not  t1=t2=0 )
!         w1>=0 w1p>=0 (but not w1=w1p=0) 
!         w2>=0 w2p>=0 (but not w2=w2p=0) 
!
!
!      
      subroutine eprelightg(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
       include "Zelightg.h"
!
!         interface to read configuration data for "elightg"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
!
!           read elightg data as 'new-*'
!           elightg has 8 volume attributes and the direction cosines
!           of x,y (1~6)
!
!             next is mandatory
        call eprpst(comp, 8, 11, 1, 6)
!
        a = Volat( comp%vol + ia)
        b = Volat( comp%vol + ib)
        t1 = Volat( comp%vol + it1)
        t2 = Volat( comp%vol + it2)
        w1 = Volat( comp%vol + iw1)
        w1p = Volat( comp%vol + iw1p)
        w2 = Volat( comp%vol + iw2)
        w2p = Volat( comp%vol + iw2p)

        if(a  .le. 0. .or. b .le. 0. ) then
           write(msg, *) comp%cn, '-th component: a=', a,
     *    ' b=', b, ' for elightg;  invalid'
           call cerrorMsg(msg, 0)
        endif
        if(t1 .lt. 0 .or. t2 .lt. 0 .or.
     *    (t1 .eq. 0. .and. t2 .eq. 0.)) then
           write(msg, *) comp%cn, '-th component: t1=', t1,
     *   ' and t2=', t2,' for elightg;  invalid'
           call cerrorMsg(msg, 0)
        endif           
        if( w1 .lt. 0. .or. w1p .lt. 0 .or.
     *      (w1 .eq. 0. .and. w1p .eq. 0.) ) then
           write(msg, *) comp%cn, '-th component: w1=', w1,
     *   ' w1p=',w1p,' for elightg;  invalid'
           call cerrorMsg(msg, 0)
        endif
        if( w2 .lt. 0. .or. w2p .lt. 0 .or.
     *      (w2 .eq. 0. .and. w2p .eq. 0.) )then
           write(msg, *) comp%cn, '-th component: w2=', w2,
     *   ' w2p=',w2p,' for elightg;  invalid'
           call cerrorMsg(msg, 0)
        endif
        if( b .le. t2 ) then
           write(msg, *) comp%cn, '-th component: b<=t2'
           call cerrorMsg(msg, 0)
        endif
        Volat( comp%vol + itan1) = (w1-w1p)/(b-t2)
        Volat( comp%vol + itan2) =-(w2-w2p)/(b-t2)
        Volat( comp%vol + izc ) = -a*sqrt(1.0d0-( (b-t2)/b)**2)+a
       end
!   ***************************************
      subroutine epbelightg(comp, pos, dir, length, icon)
       implicit none
#include "Zglobalc.h"
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"
       include "Zelightg.h"
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


       real*8 el
       type(epPos)::  org
       type(ep3Vec)::  abc
       real*8 elightgf1, elightgf2, elightgf3, elightgf4
       external elightgf1, elightgf2, elightgf3, elightgf4

       real*8 tempx, tempy, tempz, temp, leng

       real*8 xpa(5)


       integer nx, i, jcon


       nx = 0
       
       call epenvlpelightg(comp, org, abc)
!         obtain the crossing point with the bounding box
       call kxplbx( pos%x-org%x, pos%y-org%y, pos%z - org%z,
     *            dir%x, dir%y, dir%z, abc%x, abc%y, abc%z,
     *            el, jcon)

       if(jcon .eq. -1) then
          icon = -1
          return ! *********
       endif

       if(jcon .eq. 0) then
!            see if pos is inside
          call epselightg(comp, pos, icon) ! common has been set
       else
          icon = 1
       endif
!
!           see upper curved surfafe
       call epxpellipse(pos%x, pos%z, dir%x, dir%z,  b-t2, a-t1,
     *   0.d0,  a,   leng,  jcon)

       if( jcon .ne. -1) then
          tempx = pos%x + leng*dir%x 
          if( tempx .ge. 0.d0 .and. tempx .le. b-t2) then
             tempy =  pos%y + leng*dir%y
             if( tempy .ge. elightgf3(tempx) .and.
     *           tempy .le. elightgf4(tempx) ) then
                tempz = pos%z + leng*dir%z
                if(tempz .ge. t1 .and. tempz .le. a) then
                   nx = nx + 1
                   xpa(nx) = leng
                endif
             endif
          endif
       endif
!            see lower ellipse
       call epxpellipse(pos%x, pos%z, dir%x, dir%z,   b, a,
     *    0.d0, a,  leng,  jcon)
       if(jcon .ne. -1) then
          tempx = pos%x + leng*dir%x 
          if( tempx .ge. 0.d0 .and. tempx .le. b) then
             tempy =  pos%y + leng*dir%y
             if( tempy .ge. elightgf3(tempx) .and.
     *           tempy .le. elightgf4(tempx) ) then
                tempz = pos%z + leng*dir%z
                if(tempz .le.  a) then
                   nx = nx + 1
                   xpa(nx) = leng
                endif
             endif
          endif
       endif

!               see -y wall (x<b-t2)
!          pos.y + leng*dir.y = -w1 + tan1(pos.x+leng*dir.x)
!          leng *( tan1*dir.x-dir.y) = w1+pos.y - tan1*pos.x
!          
       temp = tan1*dir%x-dir%y
       if( temp .ne. 0.) then
          leng = (w1+pos%y - tan1*pos%x)/temp
          if(leng .ge. 0.) then
             tempx = pos%x + leng*dir%x 
             if(  tempx .le. b-t2 .and. tempx .ge. 0.d0) then
                tempz = pos%z + leng*dir%z
                if(tempz .ge. elightgf2(tempx) .and. 
     *             tempz .le. elightgf1(tempx) ) then
                   nx = nx + 1
                   xpa(nx) = leng
                endif
             endif
          endif
       endif
!           see  -w1p wall
!         -w1p = pos.y + leng*dir.y
       if(dir%y .ne. 0.) then
          leng =( -w1p - pos%y )/dir%y
          if(leng .ge. 0.) then
             tempx = pos%x + leng*dir%x
             if(tempx .ge. b-t2 .and. tempx .le. b ) then
                tempz = pos%z + leng* dir%z
                if(tempz .le.  a .and.
     *             tempz .ge.  elightgf2(tempx) ) then
                   nx = nx + 1
                   xpa(nx) = leng
                endif
             endif
          endif
       endif
             
!               see y >0 wall (x<b-t2)
!          pos.y + leng*dir.y = w2 + tan2(pos.x+leng*dir.x)
!          leng *(-tan2*dir.x+dir.y) = w2-pos.y + tan2*pos.x
!          
       temp = -tan2*dir%x+dir%y
       if( temp .ne. 0.) then
          leng = (w2 - pos%y + tan2*pos%x)/temp
          if(leng .ge. 0.) then
             tempx = pos%x + leng*dir%x 
             if(  tempx .le. b-t2 .and. tempx .ge. 0.d0) then
                tempz = pos%z + leng*dir%z
                if(tempz .ge. elightgf2(tempx) .and. 
     *             tempz .le. elightgf1(tempx) ) then
                   nx = nx + 1
                   xpa(nx) = leng
                endif
             endif
          endif
       endif
!
!           see  w2p wall
!           w2p = pos.y + leng*dir.y
       if(dir%y .ne. 0.) then
          leng =( w2p - pos%y )/dir%y
          if(leng .ge. 0.) then
             tempx = pos%x + leng*dir%x
             if(tempx .ge. b-t2 .and. tempx .le. b ) then
                tempz = pos%z + leng* dir%z
                if(tempz .le.  a .and.
     *             tempz .ge.  elightgf2(tempx) ) then
                   nx = nx + 1
                   xpa(nx) = leng
                endif
             endif
          endif
       endif
!          
!           x  0 wall
!
!     0 = pos.x + leng*dir.x
       if( dir%x .ne. 0.) then
          leng = -pos%x/dir%x
          if(leng .ge. 0.) then
             tempy = pos%y + leng * dir%y
             if(tempy .ge. -w1 .and. tempy .le. w2) then
                tempz = pos%z + leng *dir%z
                if(tempz .ge. 0. .and. tempz .le. t1) then
                   nx = nx + 1
                   xpa(nx) = leng
                endif
             endif
          endif
       endif
! 
!           z = a wall
!            
       if(dir%z .ne. 0.) then
          leng = (a-pos%z)/dir%z
          if(leng .ge. 0.) then
             tempx = pos%x + leng*dir%x
             if(tempx .ge. b-t2 .and. tempx .le. b) then
                tempy = pos%y + leng*dir%y
                if(tempy .ge. -w1p .and. tempy .le. w2p) then
                   nx = nx + 1
                   xpa(nx) = leng
                endif
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

      subroutine epselightg(comp, pos, icon)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
      include "Zelightg.h"
!
!           judge if a given 'pos' is inside 'comp'
!         
       type(Component)::  comp !input component
       type(epPos)::  pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside

      real*8 elightgf1
      real*8 elightgf2
      real*8 elightgf3
      real*8 elightgf4
      external elightgf1, elightgf2, elightgf3, elightgf4

      call epelightgset(comp)

      icon  = 0
      if(pos%x .lt. 0.) then
         icon = 1
      elseif(pos%x .gt. b) then
         icon = 1
      elseif(pos%x .lt. (b-t2)) then
         if( pos%y .lt. elightgf3(pos%x)) then
            icon = 1
         elseif( pos%y .gt. elightgf4(pos%x)) then
            icon = 1
         elseif( pos%z .gt. elightgf1(pos%x) ) then
            icon = 1
         elseif( pos%z .lt. elightgf2(pos%x) ) then
            icon = 1
         endif
      else
         if(pos%y .lt. -w1p .or. pos%y .gt. w2p) then
            icon = 1
         elseif(pos%z .lt. elightgf2(pos%x) ) then
            icon = 1
         elseif(pos%z .gt. a) then
            icon = 1
         endif
      endif
      end
      subroutine epenvlpelightg(comp, org, abc)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
      include "Zelightg.h"
!
!        give the envloping box of the component
!
       type(Component)::  comp  ! input.   component.
       type(epPos)::  org       ! output.  origin of the enveloping box
                               !          in local coord. 
       type(ep3Vec)::  abc      ! output.  a,b,c of the box

 
      call epelightgset(comp)     ! set common value

      org%x = 0.
      org%y = -w1
      org%z = 0.
      abc%x = b
      abc%y = max(w1+w2, w1p+w2p)
      abc%z = a 
      NVTX = 0
      end
!     *************************************
      subroutine epatlocelightg(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(8)
 
      integer i

      do i = 1, 8
         loc(i) = i
      enddo
      end
      subroutine epelightgset(comp)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
      include "Zelightg.h"
!         set common variables

       type(Component)::  comp ! input.

      integer cnsave
      save cnsave
      data cnsave/-1/

      if(comp%cn .ne. cnsave) then
         cnsave = comp%cn
         a = Volat( comp%vol + ia)
         b = Volat( comp%vol + ib)
         t1 = Volat( comp%vol + it1)
         t2 = Volat( comp%vol + it2)
         w1 = Volat( comp%vol + iw1)
         w1p = Volat( comp%vol + iw1p)
         w2 = Volat( comp%vol + iw2)
         w2p = Volat( comp%vol + iw2p)
         tan1 = Volat( comp%vol + itan1)
         tan2 = Volat( comp%vol + itan2)
         zc =  Volat( comp%vol + izc )
      endif
      end
! 
      real*8  function  elightgf1(x)
      include "Zelightg.h"
      real*8 x

!            before calling this, common must be set
!           upper ellipse
      elightgf1 =  -(a-t1)*sqrt(abs(1.d0-(x/(b-t2))**2)) + a
      end
! 
      real*8  function  elightgf2(x)
      include "Zelightg.h"
!            before calling this, common must be set
      real*8 x
!           lower ellipse
      elightgf2 = -a*sqrt(1.-(x/b)**2) + a
      end

      real*8  function  elightgf3(x)
      include "Zelightg.h"
      real*8 x
!            before calling this, common must be set
!              right plane
      if(x .lt. b-t2) then
         elightgf3 =  -w1 + tan1*x
      else
         elightgf3 =  -w1p
      endif
      end

      real*8  function  elightgf4(x)
      include "Zelightg.h"
      real*8 x
!            before calling this, common must be set
!           left plane
      if(x .lt.  b-t2) then
         elightgf4 =  w2 + tan2*x
      else
         elightgf4 =  w2p
      endif 
      end
