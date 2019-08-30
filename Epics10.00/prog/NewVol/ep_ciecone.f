!
!  cut inclined eliptic cone
!                                        
!  bottom ellipse center is at (0,0,0).  hight is directed to z.
!  top    ellipse  center is at (xt, 0, h)
!
!
!   Data format in config is:
!       ox oy oz  ra rb  h  xt  rap sa ea
!
!      where (ox,oy,oz) is the origin in the world coord.
!            ra: bottom x radius of the ellispe cm
!            rb:        y
!             h: heigth  cm
!            xt: center of the top ellipse. i.e, (xt, 0, h)
!                (bottom one is at (0, 0, 0)
!           rap: top  x radius of the ellipse
!                All ellipses inbetween 0 ~h have the same rb/ra 
!            sa: starting angl of the cut plane (deg) (0~360)
!            ea: ending //                            (0~360)
!                sa may be >ea.  sa,ea are measured at the large
!                ellipse.
!      
      subroutine eprciecone(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
       include "Zciecone.h"
!
!         interface to read configuration data 
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg


       real*8 ra, rb, h, xt, rap, rbp, sa, ea
       real*8  r, x, y, s, shift
!            8 attribute,  direction cos of x,y direction
       call eprpst(comp, 7, 16, 1, 6)
!
!           check some values
       ra = Volat( comp%vol +  ira )
       rb = Volat( comp%vol +  irb )
       rap = Volat( comp%vol +  irap )
       h = Volat( comp%vol +  ih )
       xt = Volat( comp%vol +  ixt )
       sa = Volat( comp%vol +  isa )
       ea = Volat( comp%vol +  iea )
       if(ra  .lt. 0. .or. h .le. 0. .or. 
     *      rb  .lt. 0.  .or. rap .lt. 0.  .or. (
     *      ra .eq. 0. and. rap .eq. 0. )  ) then        
          write(msg, *) comp%cn, '-th component: Ra=', ra,
     *         ' Rb=',rb,  ' h=', h,  ' for ciecone;  invalid'
          call cerrorMsg(msg, 0)
       endif
       if( sa .lt. 0. .or. sa .gt. 360.  .or.
     *      ea .lt. 0. .or. ea .gt. 360. ) then
          write(msg, *) comp%cn, '-th compoent: sa=',sa,
     *         ' or ea=',ea, ' invalid'
          call cerrorMsg(msg, 0)
       endif
!             compute const for later use.
       Volat( comp%vol +  ieps ) = rb/ra
       Volat( comp%vol +  ib ) = (ra-rap)/h
       Volat( comp%vol +  ia ) = xt/h
!           center to the sa point length
       if(ra .gt. rap) then
!               lower ellipse
          s = 1.0
          shift = 0.
       else
!               upper ellipse
          s = rap/ra
          shift = xt
       endif
       r = 1.d0/sqrt((cos(sa*Torad)/(ra*s))**2 +
     *          (sin(sa*Torad)/(rb*s))**2)
       Volat( comp%vol +  ixs) = r*cos(sa*Torad) + shift
       Volat( comp%vol +  iys) = r*sin(sa*Torad)
       r = 1.d0/sqrt((cos(ea*Torad)/(ra*s))**2 
     *   + (sin(ea*Torad)/(rb*s))**2)
       Volat( comp%vol +  ixe) = r*cos(ea*Torad) + shift
       Volat( comp%vol +  iye) = r*sin(ea*Torad)


       end
!   ***************************************
      subroutine epbciecone(comp, pos, dir, length, icon)
       implicit none
#include "Zglobalc.h"
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"
       include "Zciecone.h"
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

       integer where(5), which
       common /Zciecone/  where, which

       real*8 ra, rb, h, sa, ea
       real*8  r0, x0, c0, c1h, c2, dq, xp, error
       integer lc
       real*8  la(5)


       integer  jcon, kcon
!


       real*8   l, l1, l2
       integer  i
       complex*16  zz0, zz1, zz2, expia 
       real*8 x, y, z, p, a, b, eps, sint
       real*8 xs, ys, xe, ye
       real*8  ellip, f
       real*8  right, inout, nsmall, psmall

       f(x,y) = (ye-ys)*(x-xs) - (xe-xs)*(y-ys)
!       if  ellip <0;inside ellipse.
       ellip(x,y,z)  = (x-a*z)**2 + (y/eps)**2 - (ra-b*z)**2

       data  error/1.d-10/, nsmall/-1.0d-10/ psmall/1.0d-10/
       save error, nsmall, psmall
!
!         let  the crossing point  be
!              x = pos.x +  l* dir.x
!              y = pos.y +  l* dir.y
!              z = pos.z +  l* dir.z
!         and obtain l.

!        inside or outside
       call epsciecone(comp, pos, kcon)
!///////////
!       write(0,*) '--------------------------'
!       write(0,*)
!       write(0,*) 'boundary search of ciecone'
!       write(0,*) ' comp no=', comp.cn
!       write(0,*) ' pos is in side if next is 0: kcon=',kcon
!       write(0,*) ' dir is =', dir.x, dir.y, dir.z
!//////////////
       ra = Volat( comp%vol +  ira)
       h = Volat( comp%vol +  ih)
       a = Volat( comp%vol + ia)
       b = Volat(comp%vol+ib)
       eps = Volat(comp%vol + ieps)
       xs =  Volat(comp%vol + ixs)
       ys =  Volat(comp%vol + iys)
       xe =  Volat(comp%vol + ixe)
       ye =  Volat(comp%vol + iye)
!//////////////
!       write(0,*) ' ra=',ra,' h=',h, ' a,b=',a,b
!       write(0,*) 'eps=', eps, ' xs,ys=',xs,ys
!       write(0,*) ' xe,ye=',xe,ye
!/////////////
       lc = 0
       which = 0
!            see top/bottom  ellipse
       if(dir%z .ne. 0.) then
          l = (h - pos%z)/dir%z
!/////////////
!          write(0,*)'top   l=',l
!/////////////
          if(l .gt. 0.) then
             x = pos%x + l*dir%x
             y = pos%y + l*dir%y
             z = h
             inout = ellip(x,y,z)
             right = f(x,y) 
!///////////////
!             write(0,*) 'top inout =', inout,  ' right=',right
!////////////

             if(pos%z .gt. h .and. inout .lt. psmall  .and. 
     *         right  .gt. nsmall) then
!                 this should be the crossig point
                length  = l
                icon = 1
                return   ! ********************
             elseif( kcon  .eq. 0  .and. right .gt. nsmall
     *              .and. inout .lt. psmall) then
!                 this should be the crossig point
                length  = l
                icon = 0
                return   ! ********************
             elseif(inout .lt. psmall .and. right .gt. nsmall) then
                lc = lc + 1
                la(lc)= l
                where(lc) = 1
             endif
          endif
          l = - pos%z/dir%z
!/////////////
!          write(0,*) 'bottom  l=',l
!/////////////
          if(l .gt. 0.) then
             x = pos%x + l*dir%x
             y = pos%y + l*dir%y
             z = 0.
             inout = ellip(x,y,z)
             right = f(x,y) 
!///////////////////
!             write(0,*) ' xs, ys=',xs,ys, ' xe, ye=', xe, ye
!             write(0,*) ' a=',a, ' eps=',eps,' ra=',ra, ' b=',b
!             write(0,*) ' x,y,z=',x,y,z
!            write(0,*) ' bottom inout=', inout, ' right=', right,
!     *          'l =',l, ' kcon=',kcon
!/////////////////

             if(pos%z .lt. 0. .and. inout .lt. psmall .and. 
     *         right  .gt. nsmall) then
!                 this should be the crossig point
                length  = l
                icon = 1
                return   ! ********************
             elseif(kcon .eq. 0   .and.  inout .lt. psmall
     *         .and.   right .gt. nsmall) then
!                 this should be the crossig point
                length  = l
                icon = 0
                return   ! ********************
             elseif(inout .le. psmall .and. right .ge. nsmall) then
                lc = lc + 1
                la(lc)= l
                where(lc)=2
             endif
          endif
       endif
!          
!           see cone part
!
!         crossing point with (x-x')^2 + (y/eps)^2 = r^2
!
!          coeff. of l^2, l^1, l^0
       r0 = ra - b*pos%z
       x0 = pos%x -a*pos%z

       c2 = (dir%x - a*dir%z)**2 + (dir%y/eps)**2 
     *      - (b*dir%z)**2
       c1h =  x0*(dir%x - a*dir%z) +
     *        dir%y*pos%y/eps**2 + r0*dir%z*b 
       c0 =   x0**2 + (pos%y/eps)**2 - r0**2

       dq =  c1h**2 - c0*c2
!////////
!       write(0,*) 'corn part:  dq=',dq
!///////
       if(dq .ge. 0.) then
          dq = sqrt(dq)
          l1 = (-c1h - dq)/c2
          l2 = (-c1h + dq)/c2 
!//////////
!          write(0,*) 'corn l1=',l1,' l2=',l2
!/////////
          if(l1 .gt. 0. ) then
             z = pos%z + l1*dir%z
!/////
!             x = pos.x + l1*dir.x
!             y = pos.y + l1*dir.y
!
!             write(0,*)  'l1=',l1, ' ellip=',ellip(x,y,z)
!//////
             if( z .gt. 0 .and. z .lt. h) then
                x = pos%x + l1*dir%x
                y = pos%y + l1*dir%y
                if( f(x,y) .gt. 0.) then
                   lc = lc + 1
                   la(lc) =l1
                   where(lc) = 3
                endif
             endif
          endif
          if(l2 .gt. 0.) then
             z = pos%z + l2*dir%z
!/////
!             x = pos.x + l2*dir.x
!             y = pos.y + l2*dir.y
!             write(0,*)  'corn l2=',l2, ' ellip=',ellip(x,y,z)
!             write(0,*)  ' f=',f(x,y)
!             write(0,*)  ' x,y,z=',x,y,z
!///////////
             if( z .gt. 0 .and. z .lt. h) then
                x = pos%x + l2*dir%x
                y = pos%y + l2*dir%y
                if( f(x,y) .gt. nsmall) then
                   lc = lc + 1
                   la(lc)= l2
                   where(lc) = 4
                endif
             endif
          endif
       endif
!            cut plane; 
       if( abs(dir%z) .lt. 1.0d0) then
          zz0 = cmplx(pos%x, pos%y,8)
          sint =sqrt(1.d0-dir%z**2) 
          expia = cmplx(dir%x, dir%y,8)/sint
          zz1 = cmplx(xs, ys,8)
          zz2 = cmplx(xe, ye,8)
          call kxplsl(zz0, expia, zz1,zz2, error, p, l, jcon)
!//////////
!          write(0,*) 'cut plane  kxplsl; jcon=',jcon, ' l=',l, ' p=',p
!////////
!ccc          if( (jcon .eq. 0 .or. jcon .eq. 2) .and. l .gt. 0.) then
          if( (jcon .eq. 0 .or. jcon .eq. 2) .and. l .ge. 0.) then
!              jcon =2---> outside of the segment, but upper cut part 
!             may wider than the one at the bottom so we have to examine it.
!             l = sqrt( ( (xs+p*(xe-xs) - pos.x)**2 +
!     *            (ys+p*(ye-ys)-pos.y)**2)/(1.0d0-dir.z**2))

             l = l/sint
!               check if outside of ellipse region
             x = pos%x + l*dir%x
             y = pos%y + l*dir%y
             z = pos%z + l*dir%z
!///////
!             if(abs(x) .gt. 1.d-4) then
!                write(0,*) ' cut pl', x, zz0, expia, zz1, zz2, p,l,jcon
!             endif
!              write(0,*) 'cut pl x,y,z =', x, y, z
!             write(0,*) ' // =', xs+p*(xe-xs),
!     *            ys+p*(ye-ys)
!             write(0,*) 'cut pl.  ellip, f=', ellip(x,y,z),
!     *        f(x,y)
!////////
             if(  ellip(x,y,z) .lt. psmall  .and.
     *               z .gt. 0 .and. z .lt. h ) then
                lc = lc + 1
                la(lc) = l
                where(lc)=5
             endif
          endif
       endif
       if(lc .gt. 0) then
          length = la(1)
          which = where(1)
          do i = 2, lc
             if(length .gt. la(i) ) then
                length = la(i)
                which=where(i)
             endif
          enddo
          icon = kcon
!////////////
!          write(0,*) ' **********found icon=',icon, '********' 
!/////////////
       else
          if(kcon .eq. 0) then
             write(0,*)
     *         'there should be crossing point with ciecone'
             write(0,*) 'local coordinate'
             write(0,*) 'pos%x,y,z=',pos%x, pos%y, pos%z
             write(0,*) 'dir%x,y,z=',dir%x, dir%y, dir%z
             call epl2w(comp%cn, pos, pos)
             call epl2wd(comp%cn, dir, dir)
             write(0,*) 'world coordinate'
             write(0,*) 'pos%x,y,z=',pos%x, pos%y, pos%z
             write(0,*) 'dir%x,y,z=',dir%x, dir%y, dir%z
             stop 6789
          endif
          icon = -1
       endif
       end          
       subroutine getwhich(i)
       integer where(5), which
       common /Zciecone/  where, which
       integer i
       i = which
       end
!      **********************************
      subroutine epsciecone(comp, pos, icon)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       include "Zciecone.h"
!
!           judge if a given 'pos' is inside 'comp'
!         
       type(Component)::  comp !input component
       type(epPos)::  pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside
!

       real*8 xp, xs, ys, xe, ye, r
       real*8 f, x, y
       real*8 ra

       f(x,y) = (ye-ys)*(x-xs) - (xe-xs)*(y-ys)

       if( pos%z .gt. Volat( comp%vol +  ih) ) then
          icon =1
       elseif(pos%z .lt. 0.) then
          icon = 1
       else
          ra = Volat( comp%vol + ira)
!//////////
!          write(0,*) ' pos=',pos.x, pos,y, pos.z
!          write(0,*) ' a,b,ra=',Volat( comp.vol + ia ) ,
!     *     Volat( comp.vol + ib ), ra 
!///////
          xp = Volat( comp%vol + ia ) * pos%z
          r = ra- Volat( comp%vol + ib)*pos%z
!////////////
!          write(0,*) ' xp=',xp,' r=',r
!          write(0,*) ' inout=',
!     *        (pos.x-xp)**2 + 
!     *        (pos.y/Volat(comp.vol+ieps))**2 - r**2
!///////
          if( (pos%x-xp)**2 + 
     *        (pos%y/Volat(comp%vol+ieps))**2 .gt. r**2) then
             icon = 1
          else
             xs = Volat(comp%vol+ ixs)
             ys = Volat(comp%vol+ iys)
             xe = Volat(comp%vol+ ixe)
             ye = Volat(comp%vol+ iye)
!///////
!             write(0,*) ' f=',  f(pos.x, pos.y)
!////////
             if( f(pos%x, pos%y) .lt. 0.) then
                icon = 1
             else
                icon = 0
             endif
          endif
       endif
       end

!     **************************************
      subroutine epenvlpciecone(comp, org, abc)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       include "Zciecone.h"
!
!        give the envloping box of the component
!
       type(Component)::  comp  ! input.   component.
       type(epPos)::  org       ! output.  origin of the enveloping box
                               !          in local coord. 
       type(ep3Vec)::  abc      ! output.  a,b,c of the box



       real*8 x, y, z,  xmx, xmn, ymx, ymn
       real*8  xs, ys, xe, ye, r, xt
       real*8  ra, eps, a, b, h, l
       real*8  wx, wy
       real*8  c0, c1, c2, d
       real*8 f
       real*8  nsmall
       
       f(x,y) = (ye-ys)*(x-xs) - (xe-xs)*(y-ys)
       data  nsmall/-1.d-6/
!       data  nsmall/-1.d-8/
       save  nsmall

       ra = Volat( comp%vol +  ira)
       eps = Volat( comp%vol +  ieps)
       a = Volat( comp%vol +  ia)
       b = Volat( comp%vol +  ib)
       h = Volat( comp%vol +  ih)
       xt =  Volat( comp%vol + ixt )
       xs =  Volat(comp%vol+ixs)
       ys =  Volat(comp%vol+iys)
       xe =  Volat(comp%vol+ixe)
       ye =  Volat(comp%vol+iye)

       xmx = max(xs, xe)
       ymx = max(ys, ye)
       xmn = min(xs, xe)
       ymn = min(ys, ye)

       if(eps .gt. 1.) then
!         there is possibilty that cut part at the top ellipse  becomes largest
!          get the coordinate
!         x =lwx  y = lwy          
!  solve:        (x-xt)**2 + (y/eps)**2 = r**2  
!          (lwx-xt)**2 + (lwy/eps)**2 -r**2
!         (wx**2 + (wy/eps)**2)l**2 - 2wx xt*l + xt**2 -r**2 = 0
!           bottom l

          l = sqrt( xs**2 + ys**2 )
          wx = xs/l
          wy = ys/l
          r = ra - b*h
          c2 = wx**2 + (wy/eps)**2
          c1 = -2*wx*xt
          c0 = xt**2 - r**2
          d = c2**2 - r*c2*c0
          if(d .gt. 0.) then
             d = sqrt(d)
             l = (-c1 + d )/2/c2
             xmx = max(xmx, l*wx)
             ymx = max(ymx, l*wy)
             xmn = min(xmn, l*wx)
             ymn = min(ymn, l*wy)
          endif

          l = sqrt( xe**2 + ye**2 )
          wx = xe/l
          wy = ye/l
          c2 = wx**2 + (wy/eps)**2
          c1 = -2*wx*xt
          c0 = xt**2 - r**2
          d = c2**2 - r*c2*c0
          if(d .gt. 0.) then
             d = sqrt(d)
             l = (-c1 + d )/2/c2
             xmx = max(xmx, l*wx)
             ymx = max(ymx, l*wy)
             xmn = min(xmn, l*wx)
             ymn = min(ymn, l*wy)
          endif
       endif
!          see right point
       x = ra
       y = 0
       if(f(x,y) .ge. nsmall ) then
          xmx = max(xmx, x)
          ymx = max(ymx, y)
          xmn = min(xmn, x)
          ymn = min(ymn, y)
       endif
!         see top part
       x = 0.
       y = Volat(comp%vol+irb)
       if( f(x,y) .gt. nsmall) then
          xmx = max(xmx, x)
          ymx = max(ymx, y)
          xmn = min(xmn, x)
          ymn = min(ymn, y)
       endif          
!         see left part
       x = -ra
       y = 0.
       if( f(x,y) .gt. nsmall) then
          xmx = max(xmx, x)
          ymx = max(ymx, y)
          xmn = min(xmn, x)
          ymn = min(ymn, y)
       endif          
!         see bottom part
       x = 0
       y = -Volat(comp%vol+irb)
       if( f(x,y) .gt. nsmall) then
          xmx = max(xmx, x)
          ymx = max(ymx, y)
          xmn = min(xmn, x)
          ymn = min(ymn, y)
       endif          
       org%x = xmn
       org%y = ymn
       org%z = 0.
       abc%x = xmx - xmn
       abc%y = ymx - ymn
       abc%z = h
       NVTX = 0
      end
!     *************************************
      subroutine epatlocciecone(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
       include "Zciecone.h"
       type(Component)::  comp ! input.
      integer loc(7)
      integer i

      do i = 1, 7
         loc(i) = i
      enddo
      end
