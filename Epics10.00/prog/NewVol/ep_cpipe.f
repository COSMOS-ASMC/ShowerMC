!
!  cut pipe; 
!    #news  new-2  ccyl
! is needed beside   #news  new-1 cpipe in the config file
!                                        
!  bottom circle center is at (0,0,0).  hight is directed to Z.
!
!   Data format in config is:
!       ox oy oz  ir or  h  sa  ea (optional dir)
!
!      where (ox,oy,oz) is the origin in the world coord.
!            ir: inner radius of the cylinder  cm
!            or: outer radius of the cylinder
!            h: height of the //        cm
!           sa: starting angle (deg)  0 is the x-axis. counter clock wise.
!           ea: ending angle (deg).  sa=0 ea=360 means cyl.
!              sa may be > ea.  they are for the outer cylinder.
!      
      subroutine eprcpipe(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "cpipe"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer ih,  isa, iea, iir, ior
       integer ix1, iy1, ix2, iy2
       parameter (ih = 2,  isa=3, iea=4, iir=9, ior=1)
       parameter (ix1=5, iy1=6, ix2=7, iy2=8)


       real*8  h, sa, ea, ri, ro
!
!           read cut pipe data as 'new-*'
!           cpipe has 5 volume attributes and the direction cosines
!           of x,y (1~6)
!
!             next is mandatory
        call eprpst(comp, 5, 9, 1, 6)
!
        ri = Volat( comp%vol + 1)
        ro = Volat( comp%vol + 2)
        h = Volat( comp%vol + 3)
        sa= Volat( comp%vol + 4)
        ea= Volat( comp%vol + 5)
        
        Volat( comp%vol + iir) = ri
        Volat( comp%vol + ior) = ro
        Volat( comp%vol + ih) = h
        Volat( comp%vol + isa) = sa
        Volat( comp%vol + iea) = ea
        if(ri  .le. 0. .or. h .le. 0 .or. ro .le. ri) then
           write(msg, *) comp%cn, '-th component: ri=', ri,
     *    ' h=', h, ' ro=', ro, ' for cpipe;  invalid'
           call cerrorMsg(msg, 0)
        endif
!             compute consts for later use.(+ is counter clockwise).
        Volat( comp%vol + ix1) =  ro*cos(sa*Torad)
        Volat( comp%vol + iy1) =  ro*sin(sa*Torad)
        Volat( comp%vol + ix2) =  ro*cos(ea*Torad)
        Volat( comp%vol + iy2) =  ro*sin(ea*Torad)
        end
        subroutine epangcpipe(comp, inout, sa, ea)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
       type(Component)::  comp
       integer inout
       real*8 sa, ea
!
 
       integer ih,  isa, iea, iir, ior
       integer ix1, iy1, ix2, iy2
       parameter (ih = 2,  isa=3, iea=4, iir=9, ior=1)
       parameter (ix1=5, iy1=6, ix2=7, iy2=8)

       real*8 bma, cosba, sap

       sa = Volat( comp%vol + isa)
       ea = Volat( comp%vol + iea)
       if(inout .eq. 1) then
          bma =Volat( comp%vol + ior) *
     *      cos(  mod(ea-sa+360.d0, 360.d0)/2.d0 *Torad)
          cosba = bma/Volat( comp%vol + iir)
          if(abs(cosba) .gt.  1.d0) then
             if(mod(ea-sa+360.d0, 360.d0) .gt. 180.d0) then
                sa = 0
                ea = 359.9999999d0
             else
                sa = 0.
                ea =0.
             endif
          else
             sap = (sa+ ea)/2.d0 - acos(cosba)*Todeg
             ea =(sa + ea)- sap
             sa = sap
             ea = mod(ea+360.d0, 360.d0)
             sa = mod(sa+360.d0, 360.d0)
          endif
       endif
       end
!   ***************************************
      subroutine epbcpipe(comp, pos, dir, length, icon)
       implicit none
#include "Zglobalc.h"
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"


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

 
       integer ih,  isa, iea, iir, ior
       integer ix1, iy1, ix2, iy2
       parameter (ih = 2,  isa=3, iea=4, iir=9, ior=1)
       parameter (ix1=5, iy1=6, ix2=7, iy2=8)


       integer jcon
       real*8 sa, ea,  leng2
       type(epPos)::  xp


      
       call epbccyl(comp, pos, dir, length, icon)
       if(icon .eq. 1) then
          xp%x = pos%x + length*dir%x
          xp%y = pos%y + length*dir%y
          if(xp%x**2 + xp%y**2 .ge.
     *             Volat( comp%vol + iir)**2) goto 100
!              true x.p may be on the curved part of inner cyl.
          call epangcpipe(comp, 1, sa, ea)
          call epxpcylwall(Volat( comp%vol + iir),
     *     Volat( comp%vol + ih), sa, ea, pos, dir, leng2, jcon)
          if(jcon .eq. -1) then
             icon = -1
          else
             length = leng2
             icon = 1
          endif
       elseif(icon .eq. 0) then
          xp%x = pos%x + length* dir%x
          xp%y = pos%y + length* dir%y
          call epangcpipe(comp, 1, sa, ea)
          if(pos%x**2 + pos%y**2 .gt. 
     *         Volat( comp%vol + iir)**2) then
             call epxpcylwall(Volat( comp%vol + iir), 
     *        Volat( comp%vol + ih), sa, ea, pos, dir, leng2, jcon)
             if(jcon .eq. 0 ) then
                write(0,*) ' strange 1'
                stop 0000
             elseif(jcon .eq. -1) then
             else
                length = leng2
             endif
          else
             call epxpcylwall(Volat( comp%vol + iir),
     *          Volat( comp%vol + ih), sa, ea, pos, dir, leng2, jcon)
             if(jcon .eq. 0) then
                icon = 1
                length = leng2
             elseif(jcon .eq. 1) then
                write(0,*) ' strange 2'
                stop 999
             else
                icon = -1
             endif
          endif
       endif
 100   continue
       end          
!      **********************************
       subroutine epscpipe(comp, pos, icon)
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

       integer ih,  isa, iea, iir, ior
       integer ix1, iy1, ix2, iy2
       parameter (ih = 2,  isa=3, iea=4, iir=9, ior=1)
       parameter (ix1=5, iy1=6, ix2=7, iy2=8)


       real*8 f, x1, y1,  x2, y2, x, y

       f(x,y) = (y-y1)*(x2 - x1) -(y2-y1)*(x-x1)


       if( pos%z .lt. 0.d0 ) then
          icon = 1
       elseif( pos%z .gt. Volat( comp%vol + ih) ) then
          icon = 1
       elseif(pos%x**2+ pos%y**2 .gt. Volat( comp%vol + ior)**2) then
          icon = 1
       elseif(pos%x**2+ pos%y**2 .lt. Volat( comp%vol + iir)**2) then
          icon = 1
       else
          x1 = Volat( comp%vol + ix1)
          y1 = Volat( comp%vol + iy1)
          x2 = Volat( comp%vol + ix2)
          y2 = Volat( comp%vol + iy2)
          if( f(pos%x, pos%y) .gt. 0. ) then
             icon = 1
          else
             icon = 0 
          endif
       endif
       end
!     **************************************
      subroutine epenvlpcpipe(comp, org, abc)
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

 
       integer ih,  isa, iea, iir, ior
       integer ix1, iy1, ix2, iy2
       parameter (ih = 2,  isa=3, iea=4, iir=9, ior=1)
       parameter (ix1=5, iy1=6, ix2=7, iy2=8)



       real*8 r, h, sa, ea
       logical isinside
       real*8 x, xs, ys, xe, ye
       isinside(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)


       r = Volat( comp%vol + ior)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)
       xs =Volat( comp%vol + ix1)
       ys =Volat( comp%vol + iy1)
       xe =Volat( comp%vol + ix2)
       ye =Volat( comp%vol + iy2)

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
       org%z = 0.d0

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
       abc%z = h
       NVTX = 0
      end
!     *************************************
      subroutine epatloccpipe(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
       type(Component)::  comp ! input.
      integer loc(*)
 
       integer ih,  isa, iea, iir, ior
       parameter (ih = 2,  isa=3, iea=4, iir=9, ior=1)

       
       loc(1) = iir
       loc(2) = ior
       loc(3) = ih
       loc(4) = isa
       loc(5) = iea
       end
