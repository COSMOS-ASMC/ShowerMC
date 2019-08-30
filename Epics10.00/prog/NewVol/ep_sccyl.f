!
!  sliced cut cylinder
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
!       ox oy oz  r  h  n1x n1y n1z   n2x n2y n2z sa ea
!              10 vol. attrib.
!      where (ox,oy,oz) is the origin in the world coord.
!            r: radius of the cylinder  cm
!            h: height of the //        cm
!           n1: plain's direction cos.  (going outward)  passes z=0.
!           n2: plain's direction cos. (//)         passed      z=h.
!           sa: starting angle (deg)
!           ea: ending angle (deg)  (ea >= sa)
      subroutine eprsccyl(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
!
!         interface to read configuration data for "sccyl"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer ir, ih, in1x, in1y, in1z, in2x, in2y, in2z, isa, iea
       integer maxz, minz, isapx, isapy, isapz1, isapz2,
     *        ieapx, ieapy, ieapz1, ieapz2, imaxz, iminz

       parameter (ir = 1,  ih = 2,  in1x= 3, in1y=4, in1z=5)
       parameter (in2x= 6, in2y=7, in2z=8, maxz=9, minz=10  )
       parameter (isa = 11, iea=12, isapx=13, isapy=14, isapz1=15,
     *      isapz2=16, ieapx=17, ieapy=18, ieapz1=19, ieapz2=20,
     *      imaxz=21, iminz=22)

       real*8 r, h, sa, ea
       type(epDirec)::  n1, n2

       real*8  k2, x, y, angatmax, angatmin 
       logical isinside
       real*8 xx
       real*8 eps/1.d-4/       
       isinside(xx) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(xx-sa+360.d0, 360.d0)



!
!           read sliced cut cylinder data as 'new-*'
!           sccyl has 10 volume attributes and the direction cosines
!           of the x,y
!
        call eprpst(comp, 10, 22, 1, 6)
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

        k2 = h* n2%z

        sa = Volat( comp%vol + in2z+1 )
        ea = Volat( comp%vol + in2z+2 )
        Volat( comp%vol + isa ) = sa
        Volat( comp%vol + iea ) = ea

        if(r  .le. 0. .or. h .le. 0) then
           write(msg, *) comp%cn, '-th component: r=', r,
     *    ' h=', h, ' for sccyl;  invalid'
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
!              normalize for safety
           call epnormvec(n1)
           Volat( comp%vol + in1x) = n1%x
           Volat( comp%vol + in1y) = n1%y
           Volat( comp%vol + in1z) = n1%z
        endif
        if( abs(n2%x**2 + n2%y**2 + n2%z**2 -1.d0) .gt. eps) then
           write(msg, *) comp%cn, '-th component:',
     *  ' n2 is not normalized=',n2%x, n2%y, n2%z
        else
!              normalize for safety
           call epnormvec(n2)
           Volat( comp%vol + in2x) = n2%x
           Volat( comp%vol + in2y) = n2%y
           Volat( comp%vol + in2z) = n2%z
        endif

!             this should be in the same pos. as for the case of
!             scyl.        
        Volat( comp%vol + minz) = r* sqrt(n1%x**2 + n1%y**2)/n1%z
        Volat( comp%vol + maxz) =
     *        ( k2 + r* sqrt(n2%x**2+n2%y**2))/n2%z

        

        x = r*cos(sa*Torad)
        y = r*sin(sa*Torad)
        Volat( comp%vol + isapx ) = x
        Volat( comp%vol + isapy ) = y
        Volat( comp%vol + isapz1 ) = (0.- n1%x*x - n1%y*y)/n1%z
        Volat( comp%vol + isapz2 ) = (k2- n2%x*x - n2%y*y)/n2%z

        x = r*cos(ea*Torad)
        y = r*sin(ea*Torad)
        Volat( comp%vol + ieapx ) = x
        Volat( comp%vol + ieapy ) = y
        Volat( comp%vol + ieapz1 ) = (0.- n1%x*x - n1%y*y)/n1%z
        Volat( comp%vol + ieapz2 ) = (k2- n2%x*x - n2%y*y)/n2%z

        angatmax = atan2(n2%y, n2%x)
!             cos(angatmax) should have opposit sign as n2.x so that
!             angatmax is the max z pos.
        if(cos(angatmax) * n2%x .gt. 0.) then
           angatmax = angatmax + pi
        endif
        angatmax = angatmax*Todeg
        angatmin = atan2(n1%y, n1%x)

!             cos(angatmin) should have same sign as n1.x
        if(cos(angatmin)* n1%x .lt. 0.) then
           angatmin = angatmin + pi
        endif

        angatmin = angatmin*Todeg
!           see if such points are inside the cut region.
        if(isinside(angatmax)) then
           Volat( comp%vol + imaxz ) = Volat( comp%vol + maxz)
        else
           Volat( comp%vol + imaxz ) =max(
     *       Volat( comp%vol + isapz2), Volat( comp%vol + ieapz2))
        endif
        if(isinside(angatmin)) then
           Volat( comp%vol + iminz ) = Volat( comp%vol + minz)
        else
           Volat( comp%vol + iminz ) =min(
     *       Volat( comp%vol + isapz1), Volat( comp%vol + ieapz1))
        endif

        if( Volat( comp%vol + iminz) .gt. 
     *             Volat( comp%vol + imaxz)) then
           write(msg, *) comp%cn, '-th component:',
     *    ' sliced planes intersect each other'
           call cerrorMsg(msg, 0)
        endif  
      end
      
!     ****************************
      subroutine epbsccyl(comp, pos, dir, length, icon)
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


       type(epPos)::  p1, p2, p3, p4, xp
       real*8 x1, y1, x2, y2, x, y, l, f, r, h
       integer kcon, jcon

       integer ir, ih, in1x, in1y, in1z, in2x, in2y, in2z, isa, iea
       integer maxz, minz, isapx, isapy, isapz1, isapz2,
     *        ieapx, ieapy, ieapz1, ieapz2, imaxz, iminz

       parameter (ir = 1,  ih = 2,  in1x= 3, in1y=4, in1z=5)
       parameter (in2x= 6, in2y=7, in2z=8, maxz=9, minz=10  )
       parameter (isa = 11, iea=12, isapx=13, isapy=14, isapz1=15,
     *      isapz2=16, ieapx=17, ieapy=18, ieapz1=19, ieapz2=20,
     *      imaxz=21, iminz=22)
       

       f(x,y) = (y-y1)*(x2 - x1) -(y2-y1)*(x-x1)


       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
       call epbscyl(comp, pos, dir, length, icon)
       if(icon .eq. -1) then
!            no x
       else
          p1%x = Volat( comp%vol + isapx)
          p1%y = Volat( comp%vol + isapy)

          p2%x = Volat( comp%vol + ieapx)
          p2%y = Volat( comp%vol + ieapy)

          x1 = p1%x
          y1 = p1%y
          x2 = p2%x
          y2 = p2%y
          xp%x = pos%x + length*dir%x
          xp%y = pos%y + length*dir%y
          if(icon .eq. 1) then
!                         
!                 p3-------P4
!                 |        |
!                 |        |
!                 |        |
!                 p2------ p1
!
             if( f(xp%x, xp%y) .le. 0.) goto 100
!              check square part 
             p1%z = Volat( comp%vol + isapz1 )
             p4%x = p1%x
             p4%y = p1%y
             p4%z = Volat( comp%vol + isapz2)
          
             p3%x = p2%x
             p3%y = p2%y
             p3%z = Volat( comp%vol + ieapz2 )
          
             p2%z = Volat( comp%vol + ieapz1 )
             call epxpLand4vp(p1, p4, p3, p2, pos, dir, l, kcon, jcon)
             if(kcon .le. 4 .and. l > 0.d0 ) then
                length = l
             else
                icon = -1
             endif
          else
!               icon = 0; pos is inside cyl.
             if( f(pos%x, pos%y) .le. 0.) then
                if( f(xp%x, xp%y) .le. 0.) goto 100
             endif
!              check square part 
             p1%z = Volat( comp%vol + isapz1 )
             p4%x = p1%x
             p4%y = p1%y
             p4%z = Volat( comp%vol + isapz2 )
          
             p3%x = p2%x
             p3%y = p2%y
             p3%z = Volat( comp%vol + ieapz2 )
          
             p2%z = Volat( comp%vol + ieapz1 )
             if( f(pos%x, pos%y) .le. 0.) then
                call epxpLand4vp(p2, p3, p4, p1, pos, dir, 
     *             l, kcon, jcon)
                if(kcon .le. 4 .and. l > 0.d0 ) then
                   length = l
                else
!                   icon = -1
!                   write(0,*) ' should not come here'
                endif
             else
                call epxpLand4vp(p1, p4, p3, p2, pos, dir,
     *              l, kcon, jcon)
                if(kcon .le. 4 .and. l > 0.d0 ) then
                   icon =  1
                   length = l
                else
                   icon = -1
                endif
             endif
          endif
       endif
 100   continue
       end          
!      **********************************
      subroutine epssccyl(comp, pos, icon)
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


       integer ir, ih, in1x, in1y, in1z, in2x, in2y, in2z, isa, iea
       integer maxz, minz, isapx, isapy, isapz1, isapz2,
     *        ieapx, ieapy, ieapz1, ieapz2, imaxz, iminz

       parameter (ir = 1,  ih = 2,  in1x= 3, in1y=4, in1z=5)
       parameter (in2x= 6, in2y=7, in2z=8, maxz=9, minz=10  )
       parameter (isa = 11, iea=12, isapx=13, isapy=14, isapz1=15,
     *      isapz2=16, ieapx=17, ieapy=18, ieapz1=19, ieapz2=20,
     *      imaxz=21, iminz=22)



       real*8 f, f1, f2, r, h, x, y, z, n1x, n1y, n1z,
     *        n2x, n2y, n2z, k2, x1, y1, x2, y2

!           if point is lower part f1 > 0
       f1(x,y,z) = x*n1x + y*n1y + z*n1z 
!             if point is at upper part, f2> 0
       f2(x,y,z) = x*n2x + y*n2y + z*n2z - k2
!             if point is right for x1-->x2;  f < 0
       f(x,y) = (y-y1)*(x2 - x1) -(y2-y1)*(x-x1)

       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
 

       if( pos%z .lt. Volat( comp%vol + iminz) ) then
          icon = 1
       elseif( pos%z .gt. Volat( comp%vol + imaxz) ) then
          icon = 1
       elseif(pos%x**2+ pos%y**2 .gt. r**2) then
          icon = 1
       else
          x1 =  Volat( comp%vol + isapx)
          y1 =  Volat( comp%vol + isapy)
          x2 =  Volat( comp%vol + ieapx)
          y2 =  Volat( comp%vol + ieapy)
          if( f(pos%x, pos%y) .gt. 0.) then
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
     *            f2(pos%x, pos%y, pos%z) .le. 0. ) then
                icon = 0
             else
                icon = 1
             endif
          endif
       endif
       end
!     **************************************
      subroutine epenvlpsccyl(comp, org, abc)
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
       integer ir, ih, in1x, in1y, in1z, in2x, in2y, in2z, isa, iea
       integer maxz, minz, isapx, isapy, isapz1, isapz2,
     *        ieapx, ieapy, ieapz1, ieapz2, imaxz, iminz

       parameter (ir = 1,  ih = 2,  in1x= 3, in1y=4, in1z=5)
       parameter (in2x= 6, in2y=7, in2z=8, maxz=9, minz=10  )
       parameter (isa = 11, iea=12, isapx=13, isapy=14, isapz1=15,
     *      isapz2=16, ieapx=17, ieapy=18, ieapz1=19, ieapz2=20,
     *      imaxz=21, iminz=22)

       real*8 r, h, sa, ea

       logical isinside
       real*8 x, xs, ys, xe, ye
       isinside(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)


       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)
       xs = Volat( comp%vol + isapx)
       ys = Volat( comp%vol + isapy)
       xe = Volat( comp%vol + ieapx)
       ye = Volat( comp%vol + ieapy)

       if(isinside(180.d0)) then
          org%x = -r
       else
          org%x = min(xs, xe)
       endif
       if(isinside(270.d0) )then
          org%y = -r
       else
          org%y = min(ys, ye)
       endif

       org%z = Volat( comp%vol + iminz )

       if(isinside(0.d0)) then
          abc%x = r - org%x
       else
          abc%x = max(xs, xe) - org%x
       endif

       if(isinside(90.d0)) then
          abc%y = r - org%y
       else
          abc%y = max(ys, ye) - org%y
       endif
       abc%z = Volat( comp%vol + imaxz) - org%z
       NVTX = 0
      end
!     *************************************
      subroutine epatlocsccyl(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(10)
 
      integer i


       integer ir, ih, in1x, in1y, in1z, in2x, in2y, in2z, isa, iea
       integer maxz, minz, isapx, isapy, isapz1, isapz2,
     *        ieapx, ieapy, ieapz1, ieapz2, imaxz, iminz

       parameter (ir = 1,  ih = 2,  in1x= 3, in1y=4, in1z=5)
       parameter (in2x= 6, in2y=7, in2z=8, maxz=9, minz=10  )
       parameter (isa = 11, iea=12, isapx=13, isapy=14, isapz1=15,
     *      isapz2=16, ieapx=17, ieapy=18, ieapz1=19, ieapz2=20,
     *      imaxz=21, iminz=22)


      do i = 1, 8
         loc(i) = i
      enddo
      loc(9) = isa
      loc(10) = iea
      end
