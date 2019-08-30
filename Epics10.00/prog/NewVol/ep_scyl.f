!
!  sliced cylinder
!                                        
!  bottom circle center is at (0,0,0).  hight is directed to Z.
!                
!                |\        * n2
!                | \    *
!                |  \* (x,y,z)=(0,0,h)
!                |   \
!                |    \
!                |     |
!                |     |
!                |     | 
!                |     |
!                |     |
!                |     |
!                |    *
!                |  *\ (x,y,z)=(0,0,0)
!                |*   \ n1
!      
!       plain's eq:    r*n1 = k1.  take r=(0,0,0) -->  k1 = 0
!                      r*n2 = k2.  take r=(0,0,h) -->  k2 = h*n2z
!   Data format in config is:
!       ox oy oz  r  h  n1x n1y n1z   n2x n2y n2z
!           
!      where (ox,oy,oz) is the origin in the world coord.
!            r: radius of the cylinder  cm
!            h: height of the //        cm
!           n1: plain's direction cos.  (going outward)  passes z=0.
!           n2: plain's direction cos. (//)         passed      z=h.
      subroutine eprscyl(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
!
!         interface to read configuration data for "scyl"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer ir, ih, in1x, in1y, in1z, in2x, in2y, in2z
       integer maxz, minz
       parameter (ir = 1,  ih = 2,  in1x= 3, in1y=4, in1z=5)
       parameter (in2x= 6, in2y=7, in2z=8, maxz=9, minz=10)


       real*8 r, h
       type(epDirec)::  n1, n2
       real*8 eps/1.d-4/
!
!           read cut cylinder data as 'new-*'
!           scyl has 8 volume attributes and the direction cosines
!           of the  h (1~6)
!
!             next is mandatory
        call eprpst(comp, 8, 10, 1, 6)
!
!           check some values
        r = Volat( comp%vol + ir)
        h = Volat( comp%vol + ih)
        n1%x = Volat( comp%vol + in1x)
        n1%y = Volat( comp%vol + in1y)
        n1%z = Volat( comp%vol + in1z)

        n2%x = Volat( comp%vol + in2x)
        n2%y = Volat( comp%vol + in2y)
        n2%z = Volat( comp%vol + in2z)

        if(r  .le. 0. .or. h .le. 0) then
           write(msg, *) comp%cn, '-th component: r=', r,
     *    ' h=', h, ' for scyl;  invalid'
           call cerrorMsg(msg, 0)
        endif
        if(n2%z .le. 0.) then
           write(msg, *) comp%cn, '-th component: n2%z=', n2%z,
     *     ' must be > 0.'
           call cerrorMsg(msg, 0)
        endif
        if(n1%z .ge. 0.) then
           write(msg, *) comp%cn, '-th component: n1%z=', n1%z,
     *     ' must be < 0.'
           call cerrorMsg(msg, 0)
        endif
        if( abs(n1%x**2 +n1%y**2+n1%z**2 -1.d0) .gt. eps) then
           write(msg, *) comp%cn, '-th component:',
     *  ' n1 is not normalized=',n1%x, n1%y, n1%z
        else
!              normalize for safty
           call epnormvec(n1)
           Volat( comp%vol + in1x) = n1%x
           Volat( comp%vol + in1y) = n1%y
           Volat( comp%vol + in1z) = n1%z
        endif
        if( abs(n2%x**2 +n2%y**2+n2%z**2 -1.d0) .gt. eps) then
           write(msg, *) comp%cn, '-th component:',
     *  ' n2 is not normalized=',n2%x, n2%y, n2%z
        else
!              normalize for safty
           call epnormvec(n2)
           Volat( comp%vol + in2x) = n2%x
           Volat( comp%vol + in2y) = n2%y
           Volat( comp%vol + in2z) = n2%z
        endif
        
        Volat( comp%vol + minz) = r* sqrt(n1%x**2 + n1%y**2)/n1%z
        Volat( comp%vol + maxz) =
     *       ( h * n2%z + r* sqrt(n2%x**2+n2%y**2))/n2%z
        if(-Volat( comp%vol + minz) .lt.
     *          h-Volat( comp%vol + maxz)) then
           write(msg, *) comp%cn, '-th component:',
     *    ' sliced plane intersects each other'
           call cerrorMsg(msg, 0)
        endif  
      end
      
!     ****************************
      subroutine epbscyl(comp, pos, dir, length, icon)
       implicit none
#include "Zglobalc.h"
#include "ZepTrackp.h"
#include "ZepDirec.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"


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

       integer ir, ih, in1x, in1y, in1z, in2x, in2y, in2z
       integer maxz, minz
       parameter (ir = 1,  ih = 2,  in1x= 3, in1y=4, in1z=5)
       parameter (in2x= 6, in2y=7, in2z=8, maxz=9, minz=10)


       real*8 leng
       type(epPos)::  xp
       type(epDirec)::  n1, n2
       real*8 f1, f2, r, h,  x, y, z, n1x, n1y, n1z,
     *        n2x, n2y, n2z, k2
       integer kcon, where
!           if point is lower part f1 > 0
       f1(x,y,z) = x*n1x + y*n1y + z*n1z 
!             if point is at upper part, f2> 0
       f2(x,y,z) = x*n2x + y*n2y + z*n2z - k2

       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
       call kxplcy(pos%x, pos%y, pos%z-Volat( comp%vol + minz),
     *     dir%x, dir%y, dir%z, r,
     *      (Volat( comp%vol + maxz)-Volat( comp%vol + minz)),
     *     length,    icon,  where)
!  output:
!     length: x-ssing point is at pos + length*dir ( el>=0)
!   icon : output. 0 the point is in the cyl. length is obtained
!                  1 the point is out side of the cyl. length is
!                     obtained.
!                 -1 no x-ing point
!   where: output. 1  x-ing point is on x-y  top plane.
!                  2  //             on the side.
!                  6  //             on      bottom.
!                 -1  no x-ing point
       

       if( icon .eq. -1 ) then
!               no x-point
       else
          xp%x = pos%x + length*dir%x
          xp%y = pos%y + length*dir%y
          xp%z = pos%z + length*dir%z
          n1%x = Volat( comp%vol + in1x)
          n1%y = Volat( comp%vol + in1y)
          n1%z = Volat( comp%vol + in1z)

          n2%x = Volat( comp%vol + in2x)
          n2%y = Volat( comp%vol + in2y)
          n2%z = Volat( comp%vol + in2z)

          n1x = n1%x
          n1y = n1%y
          n1z = n1%z

          n2x = n2%x
          n2y = n2%y
          n2z = n2%z
          k2 = h* n2%z
          if(icon .eq. 1) then
!                         
             if( f1(xp%x, xp%y, xp%z) .le. 0. .and.
     *           f2(xp%x, xp%y, xp%z) .le. 0.) goto 100
!        
             if(f1(xp%x, xp%y, xp%z) .gt. 0.) then
!                check lower sliced part
                call epxpLandP(pos, dir, n1, 0.d0, leng, kcon)
!        get a crossing point of a half line with a given plane.
!      integer kcon      ! output. 0: crossing point obtained. (l>=0)
!                        !         1: crossing point is at the backside (l<0).
!                        !         2: line seems to be on the plane (l=0)
!                        !         3: line seems to be parallel to the plane
                if(kcon .ne. 0) then
                   icon  = -1
                elseif( ( pos%x + dir%x*leng)**2 + 
     *                  ( pos%y + dir%y*leng)**2  .lt. r**2) then
                   length = leng
                else
                   icon = -1
                endif
                goto 100
             else
!                check upper sliced part
                call epxpLandP(pos, dir, n2, k2, leng, kcon)
                if(kcon .ne. 0) then
                   icon = -1
                elseif( ( pos%x + dir%x*leng)**2 + 
     *                  ( pos%y + dir%y*leng)**2  .lt. r**2) then
                   length = leng
                else
                   icon = -1
                endif
                goto 100
             endif
          elseif(icon .eq. 0) then
             if(f1(pos%x, pos%y, pos%z) .le. 0.  .and.
     *          f2(pos%x, pos%y, pos%z) .le. 0. ) then
                if( f1(xp%x, xp%y, xp%z)  .le. 0. .and.
     *              f2(xp%x, xp%y, xp%z)  .le. 0. ) goto 100
!                     should cross at sliced part
                call epxpLandP(pos, dir, n1, 0.d0, leng, kcon)
                if(kcon .eq. 0) then
                   if( (pos%x+dir%x*leng)**2 
     *               + (pos%y+dir%y*leng)**2 .le. r**2) then
                      length = leng
                      goto 100
                   endif
                endif
                call epxpLandP(pos, dir, n2, k2, leng, kcon)
                if(kcon .eq. 0) then
                   length = leng
                   goto 100
                endif
                write(0,*) 'strange 1'
                icon= -1
             elseif(f1(pos%x, pos%y, pos%z) .gt. 0. ) then
                call epxpLandP(pos, dir, n1, 0.d0, leng, kcon)
                if(kcon .eq. 0) then
                   xp%x = pos%x + leng*dir%x
                   xp%y = pos%y + leng*dir%y
                   if(xp%x**2 + xp%y**2 .le. r**2) then
                      length = leng
                      icon = 1
                   else
                      icon = -1
                   endif
                else
                   icon = -1
                endif
             elseif(f2(pos%x, pos%y, pos%z) .gt. 0. ) then
                call epxpLandP(pos, dir, n2, k2, leng, kcon)
                if(kcon .eq. 0) then
                   xp%x = pos%x + leng*dir%x
                   xp%y = pos%y + leng*dir%y
                   if(xp%x**2 + xp%y**2 .le. r**2) then
                      length = leng
                      icon = 1
                   else
                      icon = -1
                   endif
                else
                   icon = -1
                endif
             else
                write(0,*)  ' stragne 0'
                icon = -1
             endif
          endif
       endif
 100   continue
       end          
!      **********************************
      subroutine epsscyl(comp, pos, icon)
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

       integer ir, ih, in1x, in1y, in1z, in2x, in2y, in2z
       integer maxz, minz
       parameter (ir = 1,  ih = 2,  in1x= 3, in1y=4, in1z=5)
       parameter (in2x= 6, in2y=7, in2z=8, maxz=9, minz=10)




       real*8 f1, f2, r, h, x, y, z, n1x, n1y, n1z,
     *        n2x, n2y, n2z, k2


!           if point is lower part f1 > 0
       f1(x,y,z) = x*n1x + y*n1y + z*n1z 
!             if point is at upper part, f2> 0
       f2(x,y,z) = x*n2x + y*n2y + z*n2z - k2

       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
 

       if( pos%z .lt. Volat( comp%vol + minz) ) then
          icon = 1
       elseif( pos%z .gt. Volat( comp%vol + maxz) ) then
          icon = 1
       elseif(pos%x**2+ pos%y**2 .gt. r**2) then
          icon = 1
       else

          n1x  = Volat( comp%vol + in1x)
          n1y  = Volat( comp%vol + in1y)
          n1z  = Volat( comp%vol + in1z)

          n2x  = Volat( comp%vol + in2x)
          n2y  = Volat( comp%vol + in2y)
          n2z  = Volat( comp%vol + in2z)

          k2 = h*n2z

          if( f1(pos%x, pos%y, pos%z) .le. 0. .and.
     *        f2(pos%x, pos%y, pos%z) .le. 0. ) then
             icon = 0
          else
             icon = 1
          endif
       endif
       end
!     **************************************
      subroutine epenvlpscyl(comp, org, abc)
      implicit none
#include "Zglobalc.h"
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


       integer ir, ih, in1x, in1y, in1z, in2x, in2y, in2z
       integer maxz, minz
       parameter (ir = 1,  ih = 2,  in1x= 3, in1y=4, in1z=5)
       parameter (in2x= 6, in2y=7, in2z=8, maxz=9, minz=10)


 


      org%x = -Volat( comp%vol + ir)
      org%y =org%x
      org%z = Volat( comp%vol + minz)
      abc%x = 2*Volat( comp%vol + ir)
      abc%y =  abc%x
      abc%z = Volat( comp%vol + maxz) - Volat( comp%vol + minz)
      NVTX = 0
      end
!     *************************************
      subroutine epatlocscyl(comp, loc)
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

