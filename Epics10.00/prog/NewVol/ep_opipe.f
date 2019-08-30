!
!  open pipe
!   #news  new-1  ocyl
!  is needed beside
!   #news  new-2  opipe
!                                        
!  bottom circle  center is at (0,0,0).  hight is directed to Z.
!
!                         
!   Data format in config is:
!       ox oy oz  ir or  h  sa  ea
!
!      where (ox,oy,oz) is the origin in the world coord.
!            ir: inner radius of the pipe  cm
!            or: outer radius of  //
!            h: height of the //        cm
!           sa: starting angle (deg)  0 is the x-axis. counter clock wise.
!           ea: ending angle (deg).  sa=0 ea=360 means pipe
!              sa may be > ea.
!      
      subroutine epropipe(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "opipe"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg
 
       integer iir, ior, ih, isa, iea,
     *         ircossa, irsinsa, ircosea, irsinea
       parameter( iir = 9, ior=1, ih = 2,  isa=3, iea=4, 
     *            ircossa=5, irsinsa=6, ircosea=7, irsinea=8  )


       real*8 ir, or, h, sa, ea
!
!           read general pipe data as 'new-*'
!           opipe has 5 volume attributes and the direction cosines
!           of x,y (1~6)
!
!             next is mandatory
       call eprpst(comp, 5, 9, 1, 6)
!          rearranage the position so that we can use
!          gcyl routine 
       ir = Volat( comp%vol + 1)
       or = Volat( comp%vol + 2)
       h = Volat( comp%vol + 3)
       sa= Volat( comp%vol + 4)
       ea= Volat( comp%vol + 5)
       Volat( comp%vol + ior) = or
       Volat( comp%vol + iir) = ir
       Volat( comp%vol + ih) = h
       Volat( comp%vol + isa) = sa
       Volat( comp%vol + iea) = ea


       if(ir .ge. or) then
          write(msg,*)  comp%cn,'-th comp invalid', 
     *         'inner radius must > outer radius for opipe'
          call cerrorMsg(msg, 0)
       endif
       if(ir  .lt. 0. .or. h .le. 0. .or.  or  .le. 0.) then
          write(msg, *) comp%cn, '-th component: ir=', ir,
     *         ' h=', h,  ' or=',or,
     *         ' for opipe;  invalid'
          call cerrorMsg(msg, 0)
       endif
       Volat( comp%vol + ircossa) = or*cos(sa*Torad)
       Volat( comp%vol + irsinsa) = or*sin(sa*Torad)
       Volat( comp%vol + ircosea) = or*cos(ea*Torad)
       Volat( comp%vol + irsinea) = or*sin(ea*Torad)
      end
!   ***************************************
      subroutine epbopipe(comp, pos, dir, length, icon)
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
 
       integer iir, ior, ih, isa, iea,
     *         ircossa, irsinsa, ircosea, irsinea
       parameter( iir = 9, ior=1, ih = 2,  isa=3, iea=4, 
     *            ircossa=5, irsinsa=6, ircosea=7, irsinea=8  )


!       real*8 sa, ea
 

       type(epPos)::  xp

       integer kcon


       real*8  leng2

       
       call epbocyl(comp, pos, dir, length, icon)
       xp%x = pos%x + length* dir%x
       xp%y = pos%y + length* dir%y
       if(icon .eq. 1) then
          if(xp%x**2 + xp%y**2 .ge. 
     *       Volat( comp%vol + iir)**2) goto 100  ! side/up/bot
!             true x.p may be on the curved part of inner cyl.
!             get the x.p there
          call epxpcylwall(Volat( comp%vol + iir), 
     *     Volat( comp%vol + ih), Volat( comp%vol + isa),
     *     Volat( comp%vol + iea), pos, dir, length, icon)
          if(icon .ne. -1) then
             icon = 1
          endif
       elseif(icon .eq. 0) then
          if( pos%x**2 + pos%y**2 .ge. 
     *         Volat( comp%vol + iir)**2 ) then
!                       see if x with inner cyl wall.
             call epxpcylwall(Volat( comp%vol + iir),
     *         Volat( comp%vol + ih),
     *         Volat( comp%vol + isa), Volat( comp%vol + iea), 
     *         pos, dir, leng2, kcon)
             if(kcon .eq. 1) then
                icon = 0
                length = leng2
             elseif(kcon .eq. 0) then
                write(0,*) ' strange 2'
                stop 999
             else
                icon = 0
             endif
          else
             call epxpcylwall(Volat( comp%vol + iir), 
     *         Volat( comp%vol + ih),
     *         Volat( comp%vol + isa), Volat( comp%vol + iea),
     *          pos, dir, leng2, kcon)
             if( kcon .eq. 0 ) then
                icon = 1
                length =leng2
             elseif(kcon .eq. 1) then
                write(0,*) 'strange 3'
                stop 999
             else
                icon = -1
             endif
          endif                
       endif
 100   continue
       end          
!      **********************************
      subroutine epsopipe(comp, pos, icon)
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

       integer iir, ior, ih, isa, iea,
     *         ircossa, irsinsa, ircosea, irsinea
       parameter( iir = 9, ior=1, ih = 2,  isa=3, iea=4, 
     *            ircossa=5, irsinsa=6, ircosea=7, irsinea=8  )

       real*8 sa, ea

       logical isinside
       real*8 x
       isinside(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)


 
      if(pos%z .lt. 0. .or. pos%z .gt. Volat( comp%vol + ih)) then
         icon = 1
      elseif(pos%x**2 + pos%y**2 .gt. Volat( comp%vol + ior)**2) then
         icon = 1
      elseif(pos%x**2 + pos%y**2 .lt. Volat( comp%vol + iir)**2) then
         icon = 1
      else
         sa = Volat( comp%vol + isa)
         ea = Volat( comp%vol + iea)
         if(.not. isinside(atan2(pos%y, pos%x)*Todeg) ) then
            icon = 1
         else
            icon = 0
         endif
      endif
      end
!     **************************************
      subroutine epenvlpopipe(comp, org, abc)
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

 
       integer iir, ior, ih, isa, iea,
     *         ircossa, irsinsa, ircosea, irsinea
       parameter( iir = 9, ior=1, ih = 2,  isa=3, iea=4, 
     *            ircossa=5, irsinsa=6, ircosea=7, irsinea=8  )


       real*8 ir, or, sa, ea, xs, ys, xe, ye, xs2, ys2, 
     *            xe2, ye2


       logical isinside
       real*8 x
       isinside(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)


       sa = Volat( comp%vol + isa)
       ea = Volat( comp%vol + iea)
       or = Volat( comp%vol + ior)
       ir = Volat( comp%vol + iir)

       xs = or*cos(sa*Torad)
       ys = or*sin(sa*Torad)
       xe = or*cos(ea*Torad)
       ye = or*sin(ea*Torad)


       xs2 = ir*cos(sa*Torad)
       ys2 = ir*sin(sa*Torad)
       xe2 = ir*cos(ea*Torad)
       ye2 = ir*sin(ea*Torad)

       if(isinside(180.d0)) then
          org%x = -or
       else
          org%x = min(xs, xs2, xe, xe2)
       endif
       if(isinside(270.d0) )then
          org%y = -or
       else
          org%y = min(ys, ys2, ye, ye2)
       endif
       org%z = 0.d0

       if(isinside(0.d0)) then
          abc%x = or - org%x
       else
          abc%x =  max( xs, xs2, xe, xe2) - org%x
       endif

       if(isinside(90.d0)) then
          abc%y = or - org%y
       else
          abc%y = max( ys, ys2, ye, ye2)  - org%y
       endif
       abc%z = Volat( comp%vol + ih)
       NVTX =  0
      end
!     *************************************
      subroutine epatlocopipe(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(5)

 
       integer iir, ior, ih, isa, iea,
     *         ircossa, irsinsa, ircosea, irsinea
       parameter( iir = 9, ior=1, ih = 2,  isa=3, iea=4, 
     *            ircossa=5, irsinsa=6, ircosea=7, irsinea=8  )

       loc(1) = iir 
       loc(2) = ior
       loc(3) = ih
       loc(4) = isa
       loc(5) = iea
      end
