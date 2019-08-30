!
!  cut elliptic cone
!                                        
!  bottom elliptic circle center is at (0,0,0).  hight is directed to Z.
!      
!  cut surface is inclined. 4 vertexes are determined by
!  crossing point of sa and ea angles with the bottom and top ellipses.
!
!   Data format in config is:
!       ox oy oz  ra rb  h rap   sa  ea
!
!      where (ox,oy,oz) is the origin in the world coord.
!            ra: x-radius of the cylinder  cm
!            rb: y-radius //
!            h: height of the //        cm
!           rap: ra--> rap at h.
!           sa: starting angle (deg)  0 is the x-axis. counter clock wise.
!           ea: ending angle (deg).  sa=0 ea=360 means cyl.
!              sa may be > ea.
!      
      subroutine eprccone(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "ZepPos.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "ccone"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*180 msg

       integer ira, irb, ih, irap, ik, isa, iea, ix0, iy0, ix1, iy1
       integer inx, iny, inz, ikp
       parameter (inx = 12, iny=13, inz=14, ikp=15)
       parameter( ira = 1,  irb = 2,  ih = 3, irap=4, ik=5,
     *      isa=6,  iea=7, ix0=8, iy0=9, ix1=10, iy1=11)




       real*8 ra, rb, h,  rap, sa, ea, r, k
       type(epPos)::  p1, p2, p3, n
       real*8  kp
!
       call eprpst(comp, 6, 15, 1, 6)
!
!           check some values
       ra = Volat( comp%vol +  ira)
       rb = Volat( comp%vol +  irb)
       h = Volat( comp%vol +  ih)
       rap = Volat( comp%vol +  irap)
       sa= Volat( comp%vol +  5)
       ea= Volat( comp%vol +  6)
       Volat( comp%vol +  isa) = sa
       Volat( comp%vol +  iea) = ea 
       if(ra  .le. 0. .or. h .le. 0. .or. sa .lt. 0. 
     *      .or. ea .gt. 360. .or. rb .le. 0. .or. rap .le. 0.) then
          write(msg, *) comp%cn, '-th component: ra=', ra,
     *    ' h=', h, ' sa=',sa, ' ea=',ea,' rb=',rb, ' rap=',rap,
     *    ' for ccone;  invalid'
          call cerrorMsg(msg, 0)
       endif
!             compute const for later use.
        k = rap/ra
        Volat( comp%vol +  ik ) = k
        r = 1.d0/sqrt((cos(sa*Torad)/ra)**2 + (sin(sa*Torad)/rb)**2)
        Volat( comp%vol +  ix0) = r*cos(sa*Torad)
        Volat( comp%vol +  iy0) = r*sin(sa*Torad)
        r = 1.d0/sqrt((cos(ea*Torad)/ra)**2 + (sin(ea*Torad)/rb)**2)
        Volat( comp%vol +  ix1) = r*cos(ea*Torad)
        Volat( comp%vol +  iy1) = r*sin(ea*Torad)

!          define the cut plane
        p1%x =  Volat( comp%vol +  ix0)
        p1%y =  Volat( comp%vol +  iy0)
        p1%z =  0.
        p2%x =  p1%x*k
        p2%y =  p1%y*k
        p2%z =  h
        p3%x =  Volat( comp%vol + ix1)*k
        p3%y =  Volat( comp%vol + iy1)*k
        p3%z = h
        call ep3p2plane(p1, p2, p3, n, kp)
        Volat( comp%vol + inx) = n%x
        Volat( comp%vol + iny) = n%y
        Volat( comp%vol + inz) = n%z
        Volat( comp%vol + ikp) = kp
       end
!   ***************************************
      subroutine epbccone(comp, pos, dir, length, icon)
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
                               !  by Volat( comp.vol +  1), etc
       type(epPos)::  pos   ! input.  position.
       type(epDirec)::  dir  ! input. direction cosinse

       real*8  length !  output length cm from pos to the boundary
       integer icon  ! output 0: length obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume


       integer ira, irb, ih, irap, ik, isa, iea, ix0, iy0, ix1, iy1
       integer inx, iny, inz, ikp
       parameter (inx = 12, iny=13, inz=14, ikp=15)
       parameter( ira = 1,  irb = 2,  ih = 3, irap=4, ik=5,
     *      isa=6,  iea=7, ix0=8, iy0=9, ix1=10, iy1=11)

       real*8  h,  k
 
       integer  jcon, kcon
!
       type(epPos)::  p1, p2, p3, p4, xp

       real*8   leng
       real*8 x, y,  z, f, nx, ny, nz, kp

       f(x,y,z) = x*nx + y*ny + z*nz - kp

       call epbcone(comp, pos, dir, length, icon)

       if(icon .eq. -1) then
       else
          h = Volat( comp%vol +  ih)
          nx = Volat( comp%vol + inx)
          ny = Volat( comp%vol + iny)
          nz = Volat( comp%vol + inz)
          kp = Volat( comp%vol + ikp)

          xp%x = pos%x + length*dir%x
          xp%y = pos%y + length*dir%y
          xp%z = pos%z + length*dir%z
!          
          p1%x = Volat( comp%vol + ix0)
          p1%y = Volat( comp%vol + iy0)
          p2%x = Volat( comp%vol + ix1)
          p2%y = Volat( comp%vol + iy1)
          k =  Volat( comp%vol + ik)
          if(icon .eq. 1) then
!
!                 p3-------P4
!                 |        |
!                 |        |
!                 |        |
!                 p2------ p1
!
             if( f(xp%x, xp%y, xp%z) .le. 0.) goto 100
!              check square part
             p1%z = 0.
             p4%x = p1%x*k
             p4%y = p1%y*k
             p4%z = h

             p3%x = p2%x*k
             p3%y = p2%y*k
             p3%z = h
             
             p2%z = 0.

             call epxpLand4vp(p1, p4, p3, p2, pos, dir, 
     *       leng, kcon, jcon)
             if(kcon .le. 4 .and. leng> 0.) then
                length = leng
             else
                icon = -1
             endif
          else
!           icon = 0;  pos is inside of the cone
!          
             if(f(pos%x, pos%y, pos%z) .le. 0.) then
                if(f(xp%x, xp%y, xp%z) .le. 0.) goto 100
             endif
!              check square part
             p1%z = 0.
             p4%x = p1%x*k
             p4%y = p1%y*k
             p4%z = h

             p3%x = p2%x*k
             p3%y = p2%y*k
             p3%z = h
             p2%z = 0.

             call epxpLand4vp(p1, p4, p3, p2, pos, dir, 
     *       leng, kcon, jcon)
             if(kcon .le. 4 .and. leng > 0.) then
                length = leng
             else
                icon = -1
             endif
          endif
       endif
 100   continue
       end          

!      **********************************
      subroutine epsccone(comp, pos, icon)
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

       integer ira, irb, ih, irap, ik, isa, iea, ix0, iy0, ix1, iy1
       integer inx, iny, inz, ikp
       parameter (inx = 12, iny=13, inz=14, ikp=15)
       parameter( ira = 1,  irb = 2,  ih = 3, irap=4, ik=5,
     *      isa=6,  iea=7, ix0=8, iy0=9, ix1=10, iy1=11)

       real*8 ra, rb, h, k
 
       real*8 x, y, z, f, nx, ny, nz, kp
       real*8  a, b, alpha

       f(x,y,z) = x*nx + y*ny + z*nz - kp






       h = Volat( comp%vol +  ih)

       if( pos%z .lt. 0.d0 ) then
          icon = 1
       elseif( pos%z .gt. h ) then
          icon = 1
       else
          ra = Volat( comp%vol +  ira)
          rb = Volat( comp%vol +  irb)
          k = Volat( comp%vol +  ik)
          alpha = (k-1.d0)/h
          a  = ra*(1.d0 + alpha*pos%z)
          b  = rb*(1.d0 + alpha*pos%z)
          if(a .eq. 0. or. b .eq. 0.) then
             icon = 1
          elseif( (pos%x/a)**2 + (pos%y/b)**2 .gt. 1.d0) then
             icon =1
          else
             nx = Volat( comp%vol + inx )
             ny = Volat( comp%vol + iny )
             nz = Volat( comp%vol + inz )
             kp = Volat( comp%vol + ikp )
             if(f(pos%x, pos%y, pos%z) .le. 0.) then
                icon = 0
             else
                icon = 1
             endif
          endif
       endif
       end
!     **************************************
      subroutine epenvlpccone(comp, org, abc)
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


       integer ira, irb, ih, irap, ik, isa, iea, ix0, iy0, ix1, iy1
       parameter( ira = 1,  irb = 2,  ih = 3, irap=4, ik=5,
     *      isa=6,  iea=7, ix0=8, iy0=9, ix1=10, iy1=11)

       real*8 ra, rb, h, sa, ea, k, rap, rbp
 


       real*8 x, xsmx, xsmn, ysmx, ysmn,  xemx, xemn, yemx, yemn
       logical isinside
       isinside(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)





       ra = Volat( comp%vol +  ira)
       rb = Volat( comp%vol +  irb)
       sa = Volat( comp%vol +  isa)
       ea = Volat( comp%vol +  iea)
       h = Volat( comp%vol +  ih)
       k = Volat( comp%vol +  ik)
       rap =Volat( comp%vol +  irap)
       rbp = rb*k
       xsmx =
     *  max( Volat( comp%vol +  ix0), Volat( comp%vol +  ix0)*k)
       xsmn =
     *  min( Volat( comp%vol +  ix0), Volat( comp%vol +  ix0)*k)
       ysmx =
     *  max(Volat( comp%vol +  iy0), Volat( comp%vol +  iy0)*k)
       ysmn =
     *  min(Volat( comp%vol +  iy0), Volat( comp%vol +  iy0)*k)
       xemx =
     *  max(Volat( comp%vol +  ix1), Volat( comp%vol +  ix1)*k)
       xemn =
     *  min(Volat( comp%vol +  ix1), Volat( comp%vol +  ix1)*k)
       yemx = 
     *  max(Volat( comp%vol +  iy1), Volat( comp%vol +  iy1)*k)
       yemn = 
     *  min(Volat( comp%vol +  iy1), Volat( comp%vol +  iy1)*k)


       if(isinside(180.d0)) then
          org%x =  min(-ra, -rap)
       else
          org%x = min(xsmn, xemn)
       endif
       if(isinside(270.d0) )then
          org%y = min(-rb, -rbp)
       else
          org%y = min(ysmn, yemn)
       endif
       org%z = 0.d0

       if(isinside(0.d0)) then
          abc%x = max(ra, rap) - org%x
       else
          abc%x = max(xsmx, xemx) - org%x
       endif

       if(isinside(90.d0)) then
          abc%y =  max( rb, rbp) - org%y
       else
          abc%y = max(ysmx, yemx) - org%y
       endif
       abc%z = h
       NVTX = 0
      end
!     *************************************
      subroutine epatlocccone(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
       type(Component)::  comp ! input.
      integer loc(*)
      integer i
       integer ira, irb, ih, irap, ik, isa, iea, ix0, iy0, ix1, iy1
       integer inx, iny, inz, ikp
       parameter (inx = 12, iny=13, inz=14, ikp=15)
       parameter( ira = 1,  irb = 2,  ih = 3, irap=4, ik=5,
     *      isa=6,  iea=7, ix0=8, iy0=9, ix1=10, iy1=11)

       do i = 1, 4
          loc(i) = i
       enddo
       loc(5) = isa
       loc(6)=  iea

       end


