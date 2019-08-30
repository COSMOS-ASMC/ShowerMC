!
!  cut skewed cone
!                                        
!  bottom  circle center is at (0,0,0).  hight is directed to Z.
!  top                      at (x2, 0, h)    
!
!
!
!   Data format in config is:
!       ox oy oz  ra rb  h  x2   sa  ea
!
!      where (ox,oy,oz) is the origin in the world coord.
!            ra: radius of the circle at bottom  cm
!            rb: radius of the top circle
!            h: height of the //        
!            x2: center of the top circle at h
!           sa: starting angle (deg)  0 is the x-axis. counter clock wise.
!           ea: ending angle (deg).  sa=0 ea=360 means full circle
!           sa may be > ea.
!      
      subroutine eprcskcone(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "ocone"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg

       integer ira, irb, ih, ix2, isa, iea, ia, ib
       parameter( ira = 1,  irb = 2,  ih = 3, ix2=4,  ik = 5,
     *      isa=6,  iea=7, ia=8, ib=9)

       integer ira, irb, ih, irap, ik, isa, iea, ix0, iy0, ix1, iy1
       parameter( ira = 1,  irb = 2,  ih = 3, irap=4, ik=5,
     *         isa=6,  iea=7, ix0=8, iy0=9, ix1=10, iy1=11)

                          


       real*8 ra, rb, h,  x2, sa, ea, 
       real*8 r, a, b, aa, bb, cc
!
       call eprpst(comp, 6, 9, 1, 6)
!
!           check some values
       ra = Volat( comp%vol +  ira)
       rb = Volat( comp%vol +  irb)
       h = Volat( comp%vol +  ih)
       x2 = Volat( comp%vol +  ix2)
       sa= Volat( comp%vol +  5)
       ea= Volat( comp%vol +  6)
       Volat( comp%vol +  isa) = sa
       Volat( comp%vol +  iea) = ea 
        if(ra  .le. 0. .or. h .le. 0. .or. sa .lt. 0. 
     *      .or. ea .gt. 360. .or. rb .lt. 0. ) then
           write(msg, *) comp%cn, '-th component: ra=', ra,
     *    ' h=', h, ' sa=',sa, ' ea=',ea,' rb=',rb,
     *    ' for cskcone;  invalid'
           call cerrorMsg(msg, 0)
        endif
!             compute const for later use.
       Volat( comp%vol +  ia ) = x2/h
       Volat( comp%vol +  ib ) = (ra-rb)/h
       Volat( comp%vol +  ik) = rb/ra
        r = 1.d0/sqrt((cos(sa*Torad)/ra)**2 + (sin(sa*Torad)/rb)**2)
        Volat( comp%vol +  ix0) = r*cos(sa*Torad)
        Volat( comp%vol +  iy0) = r*sin(sa*Torad)
        r = 1.d0/sqrt((cos(ea*Torad)/ra)**2 + (sin(ea*Torad)/rb)**2)
        Volat( comp%vol +  ix1) = r*cos(ea*Torad)
        Volat( comp%vol +  iy1) = r*sin(ea*Torad)
       end
!   ***************************************
      subroutine epbocone(comp, pos, dir, length, icon)
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

       common /Zcskcone/  where
       integer where   ! 1 top, 3; cut part sa; 4; cut part ea; 6 bottom

       integer ira, irb, ih, irap, ik, isa, iea, ix0, iy0, ix1, iy1
       parameter( ira = 1,  irb = 2,  ih = 3, irap=4, ik=5,
     *      isa=6,  iea=7, ix0=8, iy0=9, ix1=10, iy1=11)

       real*8 ra, rb, h, sa, ea
 

       integer  jcon, kcon
!
       type(epPos)::  p1, p2, p3, p4, xp

       real*8   l, angs,  angx
       logical kphiinside
       real*8 x, aa, bb, cc
!       isinside(x) = mod(ea-sa+360.d0, 360.d0) .ge.
!     *               mod(x-sa+360.d0, 360.d0)


       ra = Volat( comp%vol +  ira)
       rb = Volat( comp%vol +  irb)
       h = Volat( comp%vol +  ih)
       sa= Volat( comp%vol +  isa)
       ea= Volat( comp%vol +  iea)

       where = 0

       

       if( sa .eq. 0. and. ea .eq. 360.d0 ) then
!           full cyl. 
          return     ! **********
       elseif( icon .eq. -1 ) then
!               no x-point
          return     ! ***********
       elseif(icon .eq. 1) then
          xp%x = pos%x + length*dir%x
          xp%y = pos%y + length*dir%y
          angx = atan2(xp%y, xp%x)*Todeg
          if( kphiinside(angx, sa, ea) ) goto 100

!               x.p  may be on the cut part
          p1%x = Volat( comp%vol +  ix0)
          p1%y = Volat( comp%vol +  iy0)
          p1%z = 0.
             
          p2%x = p1%x*Volat( comp%vol +  ik)
          p2%y = p1%y*Volat( comp%vol +  ik)
          p2%z = h

          p3%x = 0.d0
          p3%y = 0.d0
          p3%z = h

          p4%x = 0.d0
          p4%y = 0.d0
          p4%z = 0.d0
          call epxpLand4vp(p1, p2, p3, p4, pos, dir, l, kcon, jcon)
          if(kcon .le. 4 .and. l > 0.d0 ) then
             length = l
             where = 3
             goto 100
          endif
!               x.p  may be on the other cut part
          p1%x = Volat( comp%vol +  ix1)
          p1%y = Volat( comp%vol +  iy1)
          p1%z = 0.
             
          p2%x = 0.
          p2%y = 0.
          p2%z = 0.
          
          p3%x = 0.d0
          p3%y = 0.d0
          p3%z = h
          
          p4%x = p1%x*Volat( comp%vol +  ik)
          p4%y = p1%y*Volat( comp%vol +  ik)
          p4%z = h
          call epxpLand4vp(p1, p2, p3, p4, pos, dir, l, kcon, jcon)
          if(kcon .le. 4 .and. l > 0.d0 ) then
             length = l
             where = 4
             goto 100
          else
             icon = -1
          endif
       else
!           icon = 0;  pos is inside of the cone
!          
          angs = atan2(pos%y, pos%x)*Todeg
          if( kphiinside(angs, sa, ea)) then
             xp%x = pos%x + length*dir%x
             xp%y = pos%y + length*dir%y
             angx = atan2(xp%y, xp%x)*Todeg
             if( kphiinside(angx, sa, ea)) goto 100
          endif
          p1%x = Volat( comp%vol +  ix0)
          p1%y = Volat( comp%vol +  iy0)
          p1%z = 0.
             
          p2%x = p1%x*Volat( comp%vol +  ik)
          p2%y = p1%y*Volat( comp%vol +  ik)
          p2%z = h

          p3%x = 0.d0
          p3%y = 0.d0
          p3%z = h

          p4%x = 0.d0
          p4%y = 0.d0
          p4%z = 0.d0
          call epxpLand4vp(p1, p2, p3, p4, pos, dir, l, kcon, jcon)
          if(kcon .le. 4 .and. l > 0.d0) then
             length = l
             where = 3
             if(.not. kphiinside(angs, sa, ea)) then
                icon = 1
             endif
             goto 100
          endif
!               x.p  may be on the cut part
          p1%x = Volat( comp%vol +  ix1)
          p1%y = Volat( comp%vol +  iy1)
          p1%z = 0.
             
          p2%x = 0.
          p2%y = 0.
          p2%z = 0.

          p3%x = 0.d0
          p3%y = 0.d0
          p3%z = h

          p4%x = p1%x*Volat( comp%vol +  ik)
          p4%y = p1%y*Volat( comp%vol +  ik)
          p4%z = h
          call epxpLand4vp(p1, p2, p3, p4, pos, dir, l, kcon, jcon)
          if(kcon .le. 4 .and. l > 0.d0) then
             length = l
             where = 4
             if(.not. kphiinside(angs, sa, ea)) then
                icon = 1
             endif
          else
             icon = -1
          endif                      
       endif
 100   continue
       end          
       subroutine epqocone(whichpart)
       implicit none
!         inquire the pos. of x.p 
       integer whichpart ! output. 1; top, 
                        !         2: side (full cyl)
                        !         3: cut sa:
                        !         4: cut ea;
                        !         6: bottom  
        common /Zocone/  where
        integer where   ! 1 top, 3; cut part sa; 4; cut part ea; 6 bottom
       
       whichpart = where
       end

!      **********************************
      subroutine epsocone(comp, pos, icon)
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
       parameter( ira = 1,  irb = 2,  ih = 3, irap=4, ik=5,
     *      isa=6,  iea=7, ix0=8, iy0=9, ix1=10, iy1=11)

       real*8 ra, rb, h, sa, ea, k
 

!

       real*8   ang, alfa, a0, b0, a, b

       logical kphiinside
       real*8 x

       ra = Volat( comp%vol +  ira)
       rb = Volat( comp%vol +  irb)
       h = Volat( comp%vol +  ih)
       sa= Volat( comp%vol +  isa)
       ea= Volat( comp%vol +  iea)
       h = Volat( comp%vol +  ih)

       if( pos%z .lt. 0.d0 ) then
          icon = 1
       elseif( pos%z .gt. h ) then
          icon = 1
       else
          k = Volat( comp%vol +  ik)
          alfa = (k-1.d0)/h
          a0 =  Volat( comp%vol +  ira)
          b0 =  Volat( comp%vol +  irb)
          a  = a0*(1.d0 + alfa*pos%z)
          b  = b0*(1.d0 + alfa*pos%z)
          if(a .eq. 0. or. b .eq. 0.) then
             icon = 1
          elseif( (pos%x/a)**2 + (pos%y/b)**2 .gt. 1.d0) then
             icon =1
          else
             ang = atan2(pos%y, pos%x)* Todeg
             if(kphiinside(ang, sa, ea)) then
                icon =0
             else
                icon = 1
             endif
          endif
       endif
       end
!     **************************************
      subroutine epenvlpocone(comp, org, abc)
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

       real*8 ra, rb, h, sa, ea, k, rap
 

       logical kphiinside
       real*8 x, xsmx, xsmn, ysmx, ysmn,  xemx, xemn, yemx, yemn


       ra = Volat( comp%vol +  ira)
       rap = Volat( comp%vol +  irap)
       rb = Volat( comp%vol +  irb)
       h = Volat( comp%vol +  ih)
       sa= Volat( comp%vol +  isa)
       ea= Volat( comp%vol +  iea)
       k = Volat( comp%vol +  ik)
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



       if(kphiinside(180.d0, sa, ea)) then
          org%x =  min(-ra, -rap)
       else
          org%x = min(xsmn, xemn, 0.d0)
       endif
       if(kphiinside(270.d0, sa, ea) )then
          org%y = min(-rb, -rb*k)
       else
          org%y = min(ysmn, yemn, 0.d0)
       endif
       org%z = 0.d0

       if(kphiinside(0.d0, sa, ea)) then
          abc%x = max(ra, rap) - org%x
       else
          abc%x = max(xsmx, xemx, 0.d0) - org%x
       endif

       if( kphiinside(90.d0, sa, ea)) then
          abc%y =  max( rb, rb*k) - org%y
       else
          abc%y = max(ysmx, yemx, 0.d0) - org%y
       endif
       abc%z = h
       NVTX = 0
      end
!     *************************************
      subroutine epatlococone(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(6)
 
      integer i



       integer ira, irb, ih, irap, ik, isa, iea
       parameter( ira = 1,  irb = 2,  ih = 3, irap=4, ik=5,
     *      isa=6,  iea=7)

      do i = 1, 4
         loc(i) = i
      enddo
      loc(5) = isa
      loc(6) = iea
      end
