      module honeycomb
      integer,save::Compnum=-1000
      real(8),save::LL, LI, LS, h, W, Xs, a, Ys, b
      real(8),save:: sint
      real(8),save:: cost
      real(8),save:: tant

      real(8),save:: HSx  ! = Ls+ LI*cost  
      real(8),save:: Sx  ! = 2HSx
      real(8),save:: HSy  ! = LI*sint +W
      real(8),save:: Sy  ! = 2*HSy
      real(8),save:: offx 
      real(8),save:: HW
      real(8),save:: HWsec
      real(8),save:: LIsint, LIcost
      real(8),save:: HLS  ! = LS/2
      real(8),save:: Yu, Yl  ! horizontal edge upper lower y
      real(8),save:: Xl, Xl0  ! f1i(Xl) =Yl f1i(Xl0)=0
      real(8),save:: Xu, Xu0  ! f1o(Xl) =Yu f1o(Xu0)=0
      real(8),save:: Xuyw     ! 
      real(8),save:: grad  !  gradient at connection part of 1 to 2 
                           ! of 1 to neighbor 3
      real(8),parameter::eps= 1.d-8 
      real(8),parameter::eps2= 1.d-4  ! 1 micron m
      real(8),save:: Xc, Yc
      integer,save:: inpos      ! when pos is in honeycom, it's position
                      ! is classified by 1,..6. from right wing
                      ! with anti-clockwise.  inpos=6  is 
                   !         2 
                   !         --
                   !     3 /    \ 1
                   !     4 \    / 6 
                   !         --5
                   !

      contains

      function f1i(x) result(ans)  ! inner slant wall
      implicit none
      real(8),intent(in)::x
      real(8)::ans
      ans=-tant*(x-HLS) + Yl + HW - HWsec
      end      function f1i
      function f1o(x) result(ans)  !outer slat wall
      implicit none
      real(8),intent(in)::x
      real(8)::ans
      ans=-tant*(x-HLS) + Yl +  HW + HWsec
      end   function f1o

      function f3i(x) result(ans)
      implicit none
      real(8),intent(in)::x
      real(8)::ans
      ans = tant*(x+HLS) + Yl +  HW - HWsec
      end   function f3i

      function f3o(x) result(ans)
      implicit none
      real(8),intent(in)::x
      real(8)::ans
      ans = tant*(x+HLS) + Yl +  HW + HWsec
      end   function f3o

      function f4i(x) result(ans)
      implicit none
      real(8),intent(in)::x
      real(8)::ans
      ans=-tant*(x+HLS) - Yl - HW + HWsec      
      end   function f4i

      function f4o(x) result(ans)
      implicit none
      real(8),intent(in)::x
      real(8)::ans
      ans=-tant*(x+HLS) - Yl - HW - HWsec
      end function f4o

      function f6i(x) result(ans)
      implicit none
      real(8),intent(in)::x
      real(8)::ans
      ans = tant*(x-HLS) - Yl - HW + HWsec
      end function f6i

      function f6o(x) result(ans)
      implicit none
      real(8),intent(in)::x
      real(8)::ans
      ans = tant*(x-HLS)- Yl- HW  - HWsec
      end       function f6o

      function g12(x) result(ans)
!           line (Xl,Yl) --(Xu, Yu)
!           grad = (Yu-Yl)/(Xu-Xl)  = sin/(1-cos)
      implicit none
      real(8),intent(in)::x
      real(8)::ans
      ans = grad*(x-HLS) +  Yl + HW
      end      function g12

      function g1b(x) result(ans)  ! not used
!           line (Xl0,0) --(Xuyw, W)
!           grad =  sin/(1-cos)
      implicit none
      real(8),intent(in)::x
      real(8)::ans
      ans = grad*(x-Xl0)
      end   function g1b

      function nodal_y(x, ny) result(ans)
      implicit none
      real(8),intent(in):: x  ! nodal point
      integer,intent(in):: ny ! y band region.
      real(8)::ans(2)  !  ans(1)=x; ans(2) y at x.   x may be
                       ! modifed to gsep or other point

      real(8)::xp, yp, xo, yo
      integer:: n
      n = (x+eps2)/Sx
      xp = x - n*Sx
      yp = ny*HSy
      xo = x
      if( mod(ny,2) == 0 ) then
         if(xp < HLS) then
            yo = HW + yp
         elseif( xp < offx ) then
            yo = f6i(xp) - HWsec + yp + HSy
         elseif( xp < offx +LS ) then
            yo = HSy + yp - HW
         else
            xp = xp - LIcost-LS
            yo = f1i(xp)  + HWsec + yp 
         endif
      else
         if(xp < HLS ) then
            yo =   yp + HSy - HW
         elseif( xp < offx ) then
            yo = f1i(xp) + HWsec + yp 
         elseif( xp < offx+LS ) then
            yo =  yp + HW
         else
            xp = xp - LIcost-LS
            yo = f6i(xp) - HWsec + HSy+  yp
         endif
      endif
      ans(1) = x
      ans(2) = yo
      end  function nodal_y


!                
!               
      end module honeycomb
!  
!
!   Data format in config is:
!       ox oy oz  LL, LI, LS, h, W,  Xs, a, Ys, b
!    
!           LL: if negative, |LL| is angle in deg. LL>0 length
!           W>0, h>0,  0<=Xs< width of honeycomb unit
!  a>0.   0<=Ys< unit for y
!  b>0 
!
!      
      subroutine eprhoneycomb(comp)
      use honeycomb
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!
       type(Component)::  comp  ! output. to recieve the config data.


!           9 volume attributes and the direction cosines
!           (1~6)
!
!             next is mandatory
        call eprpst(comp, 9, 9, 1, 6)
!
        call  ep_honeycmCnst(comp)

        if(LI  <= 0. .or. LS <= 0 .or.  h <= 0. .or.
     *     W<=0. .or. Xs >= Sx+offx
     *    .or.  Xs < 0. .or. Ys < 0   ) then
           write(0, *)
     *      comp%cn, '-th component: some of LL=', LL,
     *    ' LI=',LI, ' LS=',LS, ' h=',h,
     *    ' W=', W,  ' Xs=',Xs, ' Ys=',Ys,
     *    ' for honeycomb  invalid'
           stop 98765
        endif
        if( a <= 0 .or. b <= 0. ) then
           write(0, *)
     *      comp%cn, '-th component: some of a=', a,
     *    ' b=', b, ' for honeycomb  invalid'
           stop 98766
        endif
       end subroutine eprhoneycomb
!   ***************************************
      subroutine epbhoneycomb(comp, pos, dir, length, icon)
      use honeycomb
      implicit none
#include "Zglobalc.h"
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"

!
!    find length to the boundary of a honeycomb specfied by 
!    'comp' from 'pos' with direction cos 'dir'
!    'pos' and 'dir' are given in this 'comp' local coordinate.
!     pos  may or may not be inside comp.
! if pos is outside of B.B, and partilce can cross B.B, the
! length is to the XP. So to get length to the honeycomb,
! this may be called again.  (then, no xp may happen)
!
       type(Component):: comp    ! input.  honeycomb

       type(epPos)::  pos        ! input.  position.
       type(epDirec)::  dir      ! input. direction cosinse

      real(8),intent(out):: length    !  length cm from pos to the boundary
      integer,intent(out):: icon  
          !  0: length obtained. pos  is inside
          !  1:  //                     outside ( xp is to the
          !        B.B not to the honeycomb)
          ! -1: the line dose not cross the volume


       type(epPos)::  cpos        
       type(epDirec)::  cdir      

      call epv2c_honeycomb(comp, pos, cpos) !  local to canonical
      call epv2cd_honeycomb(comp, dir, cdir)  ! both pos and dir
!        see if there is possibilty of xrossing.
      call epbhoneyBB(comp, cpos, cdir, length,icon)
      if( icon == 0 ) then  ! yes; so examine in detail
         call epbhoneycm0(comp, cpos, cdir, length, icon)
      else
!          no xsing /or xing with B.B
      endif
      end  subroutine epbhoneycomb

      subroutine epbhoneyBB(comp, pos, dir, length, icon)
!       see if the line starting from pos with dir can
!       cross the bounding box of the honeycomb
      use honeycomb
      implicit none
#include "Zglobalc.h"
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"

       type(Component):: comp    ! input.  honeycomb
     
      real(8),intent(in):: pos(3)  !  position.  canonical
      real(8),intent(in):: dir(3) !   direction cosinse. canonical
      real(8),intent(out):: length ! length to the XP, if exists 
      integer,intent(out):: icon  !   0 cross with the bb.
                      !  -1   no  cross. so
                      !   no possibility xsing honeycomb

      logical:: kInsideBox  ! external func.
      real(8)::  org(3), abc(3), spos(3)
      call epenvlphoneycomb(comp, org, abc)
! if config is:
! honeycomb_zx  SCIN    0  0  0 /  7   5  6  -60 1 0.6  10  0.05 0\
! .5 60  0.3 30
!      org becomes: 0.5 0.3 0  (i.e, Xs, Ys, 0)
      call epv2c_honeycomb(comp, abc, abc)
!         abc becomes 60 30 10
!        to see x.p existence,
!        bounding box origin should be 0; so if Xs Ys is non zero
!        we must shit particle position
!      spos(1) = pos(1) - Xs
!      spos(2) = pos(2) - Ys
!      spos(3) = pos(3)  
      spos(:) = pos(:) - org(:)   ! same as abvoe
      if( kInsideBox(spos(1), spos(2), spos(3),
     *   abc(1), abc(2), abc(3) ) )  then
         icon = 0
         length = 0.  ! dummy
      else
!           search xsing point with the b.b
         call kxplbx(spos(1), spos(2), spos(3), dir(1), dir(2), 
     *     dir(3), abc(1), abc(2), abc(3), length, icon)
      endif
      end  subroutine epbhoneyBB


      subroutine epbhoneycm0(comp, pos, dir, length, icon)
      use honeycomb
      implicit none
#include "Zglobalc.h"
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"

!
!    find length to the boundary of a honeycomb specfied by 
!    'comp' from 'pos' with direction cos 'dir'
!    'pos' and 'dir' are given in this 'comp' local coordinate.
!     pos  may or may not be inside comp.

       type(Component):: comp    ! input.  honeycomb

       type(epPos)::  pos        ! input.  position.
       type(epDirec)::  dir      ! input. direction cosinse

      real(8),intent(out):: length    !  length cm from pos to the boundary
      integer,intent(out):: icon      ! t 0: length obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume
      real(8)::xb,yb,zb
      real(8)::leng(4), minlen
      integer:: lengc, i

      call epbhoneycm1(comp, pos, dir, length, icon)
      if(icon ==  0 ) then
!            the boundary may be out side of the square
!         region formed by (Xs,Ys) (Xs+s, Ys+b)
!         (since search is done without such boundary)
!         if (xb, yb) below is outside, we must reset
!         the boundary
         xb = pos%x + length*dir%x
         yb = pos%y + length*dir%y
         lengc = 0
         if(xb < Xs ) then
            lengc = lengc + 1
            leng(lengc) = (Xs-pos%x)/dir%x  ! dir.x should not be 0
         elseif(xb > Xs+a) then
            lengc = lengc + 1
            leng(lengc) = ((Xs+a)-pos%x)/dir%x
         endif
         if(yb < Ys ) then
            lengc = lengc + 1
            leng(lengc) = (Ys-pos%y)/dir%y
         elseif(yb > Ys+b) then
            lengc = lengc + 1
            leng(lengc) = (Ys+b-pos%y)/dir%y
         endif
         if(lengc > 0 ) then
            minlen = 1.d20
            do i = 1, lengc
               if(leng(i) >= 0.) then
                  if(minlen > leng(i)) then
                     minlen = leng(i)
                  endif
               endif
            enddo
            length = minlen
         endif
      endif
      end  subroutine epbhoneycm0

!   ***************************************
      subroutine epbhoneycm1(comp, pos, dir, length, icon)
      use honeycomb
      implicit none
#include "Zglobalc.h"
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"

!
!    find length to the boundary of a honeycomb specfied by 
!    'comp' from 'pos' with direction cos 'dir'
!    'pos' and 'dir' are given in this 'comp' local coordinate.
!     pos  may or may not be inside comp.

       type(Component):: comp    ! input.  honeycomb

       type(epPos)::  pos        ! input.  position.
       type(epDirec)::  dir      ! input. direction cosinse

      real(8),intent(out):: length    !  length cm from pos to the boundary
      integer,intent(out):: icon      !  0: length obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume

      real(8)::  el
       type(epPos)::  p1, p2, p3, p4
       type(epPos)::  OO, AA, BB, CC
       type(epPos)::  OOp, AAp, BBp, CCp
       type(epPos)::  xyz
       type(epPos)::  normal
       type(epDirec)::  dcos

      real(8):: x, y

      real(8):: xb, yb, zb
      real(8):: xa(4), ya(4)
      integer:: kcon, jcon
      real(8):: dd
      integer i

!          see if pos is inside or out side of honeycomb
!          inpos is also set if inside
      call epshoneycm0(comp, pos, icon)   ! this should be cm0. 
!          move to unit honeycomb coordinate.
      x = pos%x - Xc
      y = pos%y - Yc

      if(icon == 1 ) then
         ! pos is outside; see if it crosses top or bottom
         ! (prob. large)
         !      z =  pos.z + el*dir.z  = 0.
         if( dir%z >0. ) then
               !  we can assume pos.z >=0 and pos.z<=h
            el = (h-pos%z)/dir%z
            xb = x + el*dir%x
            yb = y + el*dir%y
           ! see if xb,yb is in inner hexagon
            call ep_honeyin2(xb, yb,  kcon) 
            if( kcon == 0 ) then
               ! inner hexagon. so does not cross 
               icon = -1
               return
            endif
         elseif(dir%z < 0.) then
            el = -pos%z/dir%z
            xb = x + el*dir%x
            yb = y + el*dir%y
           ! see if xb,yb is in inner hexagon
            call ep_honeyin2(xb, yb,  kcon) 
            if( kcon == 0 ) then
               ! inner hexagon. so does not cross 
               icon = -1
               return
            endif
         else
            !  horizontal beam; considered below
         endif
         ! see the crossing points with hexagon wall
         !  No. 2 wall
         normal = epPos(0.d0, 1.d0, 0.d0)
         dd = Yl
         xyz = epPos(x, y, pos%z)
!           inner vertex  
         p1 = epPos(Xl, Yl, 0.d0)
         p2 = epPos(Xl, Yl, h)
         p3 = epPos(-Xl, Yl, h)
         p4 = epPos(-Xl, Yl, 0.d0)
         call epxpLand4vp2(p1, p2, p3, p4,
     *       xyz, dir, normal, dd, el, kcon, jcon)
!///////////////////
!            write(0,'(a,i3,1p,4G12.4)') 'w2 kcon, el', kcon, el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!///////////////////
         if(kcon .le. 4 .and. el >  0.d0 ) then
            length = el
            return
         endif

         !       wall 5
         normal = epPos(0.d0, -1.d0, 0.d0)
         dd = Yl
!           inner vertex
         p1 = epPos(Xl, -Yl, 0.d0)
         p2 = epPos(Xl, -Yl, h)
         p3 = epPos(-Xl, -Yl, h)
         p4 = epPos(-Xl, -Yl, 0.d0)
         call epxpLand4vp2(p1, p2, p3, p4,
     *       xyz, dir, normal, dd, el, kcon, jcon)
!///////////////////
!            write(0,'(a,i3,1p,4G12.4)') 'w5 kcon el',kcon, el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!///////////////////
         if(kcon .le. 4 .and. el > 0.d0) then
            length = el
            return
         endif

!          no 1 wall
         normal = epPos(sint, cost, 0.d0)
         dd = Xl*sint+ Yl*cost
!           wall inner vertex
!           on x axis; f1b(x)=0
         p1 = epPos(Xl0, 0.d0, 0.d0)
         p2 = epPos(Xl0, 0.d0, h)
         p3 = epPos(Xl, Yl, h)
         p4 = epPos(Xl, Yl, 0.d0)
         call epxpLand4vp2(p1, p2, p3, p4,
     *       xyz, dir, normal, dd, el, kcon, jcon)
!///////////////
!            write(0,'(a,i2,1p,4G12.4)') 'w1 kcon, el',kcon, el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!///////////////////

         if(kcon .le. 4 .and. el > 0.d0 ) then
            length = el
            return
         endif

         
!          no 4 wall
         normal = epPos(-sint, -cost, 0.d0)
         dd = Xl*sint+ Yl*cost
!           wall inner vertex
         p1 = epPos(-Xl0, 0.d0, 0.d0)
         p2 = epPos(-Xl0, 0.d0, h)
         p3 = epPos(-Xl, -Yl, h)
         p4 = epPos(-Xl, -Yl, 0.d0)
         call epxpLand4vp2(p1, p2, p3, p4,
     *       xyz, dir, normal, dd, el, kcon, jcon)
!///////////////
!            write(0,'(a,i3,1p,4G12.4)') 'w4 kocn, el=',kcon,el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!///////////////////
         if(kcon .le. 4 .and. el >  0.d0) then
            length = el
            return
         endif

         
         
!          no 3 wall;  cos(90+t), sin
         normal = epPos(-sint, cost, 0.d0)
         dd = Xl*sint+ Yl*cost  ! (-Xl,Yl,0)*(-sint,cost,0)
!           wall inner vertex
!           on x axis; f1b(x)=0
         p1 = epPos(-Xl, Yl, 0.d0)
         p2 = epPos(-Xl, Yl, h)
         p3 = epPos(-Xl0, 0.d0, h)
         p4 = epPos(-Xl0, 0.d0, 0.d0)
         call epxpLand4vp2(p1, p2, p3, p4,
     *       xyz, dir, normal, dd, el, kcon, jcon)
!///////////////
!            write(0,'(a,i3,1p,4G12.4)') 'w3,kcon,el=',kcon, el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!///////////////////
         if(kcon .le. 4 .and. el > 0.d0) then
            length = el
            return
         endif

!          no 6 wall;  
         normal = epPos(sint, -cost, 0.d0)
         dd = Xl*sint+ Yl*cost  ! (Xl,-Yl,0)*(sint,-cost,0)
!           wall inner vertex
!           on x axis; f1b(x)=0
         p1 = epPos(Xl, -Yl, 0.d0)
         p2 = epPos(Xl, -Yl, h)
         p3 = epPos(Xl0, 0.d0, h)
         p4 = epPos(Xl0, 0.d0, 0.d0)
         call epxpLand4vp2(p1, p2, p3, p4,
     *       xyz, dir, normal, dd, el, kcon, jcon)
!///////////////
!            write(0,'(a,i3,1p,4G12.4)') 'w6 kcon,el= ',kcon, el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!///////////////////
         if(kcon .le. 4 .and. el > 0.d0) then
            length = el
            return
         endif
         write(0,*) 'should not come here: inside epbhoneycm1'
         stop
      else
         icon = 0 
!         pos is inside of honeycomb walls
         OOp = epPos(Xl, Yl, 0.d0)
         AAp = epPos(Xu, Yu, 0.d0)
         BBp = epPos(Xu, Yu, h)
         CCp = epPos(Xl, Yl, h)
         xyz = epPos(x, y, pos%z)
         dcos = dir
         if(inpos == 2 .or. inpos == 5) then
            if(inpos == 5) then
               xyz%y = -xyz%y
               dcos%y = -dcos%y
            endif
            OO=epPos(-Xl, Yl, 0.d0)
            AA=epPos(-Xl, Yl, h)
            BB=epPos(-Xu, Yu, h)
            CC=epPos(-Xu, Yu, 0.d0)
!              inner horizontal wall
            normal = epPos(0.d0, 1.d0, 0.d0)
            dd = Yl
            call epxpLand4vp2(OOp, CCp, AA, OO,
     *       xyz, dcos, normal, dd, el, kcon, jcon)
            if(kcon .le. 4 .and. el > 0.d0 ) then
               length = el
               return
            endif
!              outer horizontal wall
            dd = Yu
            call epxpLand4vp2(AAp, BBp, BB, CC,
     *       xyz, dcos, normal, dd, el, kcon, jcon)
!////////////////////////
!               write(0, '(a,i3,1p,3g12.4)') 
!     *          'ohw ',inpos, 
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!////////////////////////
            if(kcon .le. 4 .and. el > 0.d0 ) then
               length = el
               return
            endif
!              2 --> 1 / wall
            normal =epPos(-sint/(1.d0-cost)*sqrt((1.0+cost)/2),
     *              sqrt( (1.0+cost)/2.), 0.d0)
            dd = Xl*normal%x + Yl*normal%y

            call epxpLand4vp2(AAp, BBp, CCp, OOp,
     *       xyz, dcos, normal, dd, el, kcon, jcon)
!////////////////////////
!               write(0, '(a,i3,1p,3g12.4)') 
!     *          'hw21 ',inpos, 
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!////////////////////////
            if(kcon .le. 4 .and. el > 0.d0) then
               length = el
               return
            endif
!            2 --> 3 / wall
            normal =epPos(sint/(1.d0-cost)*sqrt((1.0+cost)/2),
     *              sqrt( (1.0+cost)/2.), 0.d0)
            dd =-Xl*normal%x + Yl*normal%y

            call epxpLand4vp2(CC, BB, AA, OO,
     *       xyz, dcos, normal, dd, el, kcon, jcon)
!////////////////////////
!               write(0, '(a,i3,1p,3g12.4)') 
!     *          'hw23 ',inpos, 
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!////////////////////////
            if(kcon .le. 4 .and. el > 0.d0) then
               length = el
               return
            endif
!              horizontal bottom wall
            normal =epPos(0.d0, 0.d0, -1.d0)
            dd = 0.

            call epxpLand4vp2(AAp, CC, OO, OOp,
     *       xyz, dcos, normal, dd, el, kcon, jcon)
!////////////////////////
!               write(0, '(a,i3,1p,3g12.4)') 
!     *          'hb ',inpos, 
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!////////////////////////
            if(kcon .le. 4 .and. el > 0.d0) then
               length = el
               return
            endif

!            horizontal top wall
            normal =epPos(0.d0, 0.d0, 1.d0)
            dd = h
            call epxpLand4vp2(BBp, BB, AA, CCp,
     *       xyz, dcos, normal, dd, el, kcon, jcon)
!////////////////////////
!               write(0, '(a,i3,1p,3g12.4)') 
!     *          'ht ',inpos, 
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!////////////////////////
            if(kcon .le. 4 .and. el > 0.d0 ) then
               length = el
               return
            endif
         else
         ! see  1 3 4 6, use wall 1
         ! actually 4 and 6 will never come
         !  they are treated as 1 or 3 of another
         !  cell 
            OO = epPos(Xl0, 0.0d0, 0.0d0)
            AA = epPos(Xu0, 0.d0, 0.0d0)
            BB = epPos(Xu0, 0.d0, h)
            CC = epPos(Xl0, 0.0d0, h)

            if(inpos == 1) then
            elseif( inpos == 3 ) then
!               mirror left right
               xyz%x = - xyz%x
               dcos%x  =  -dcos%x
            elseif( inpos == 4 ) then
               xyz%x = - xyz%x
               xyz%y = - xyz%y
               dcos%x = - dcos%x
               dcos%y = - dcos%y
            elseif( inpos == 6 ) then
               xyz%y = - xyz%y
               dcos%y = -dcos%y
            endif

            !    inner  wall by f1i
            normal = epPos(sint, cost, 0.d0)
            dd = Xl*sint+ Yl*cost
            call epxpLand4vp2(OO, CC, CCp, OOp,
     *       xyz, dcos, normal, dd, el, kcon, jcon)
!////////////////////////
!               write(0, '(a,i3,1p,4g12.4)') 
!     *          'iw ',inpos, el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!////////////////////////
            if(kcon .le. 4 .and. el > 0.d0) then
               length = el
               return
            endif
            !    outer  wall by f1o
            dd = Xu*sint+ Yu*cost
            call epxpLand4vp2(AA, BB, BBp, AAp,
     *       xyz, dcos, normal, dd, el, kcon, jcon)
!////////////////////////
!               write(0, '(a,i3,1p,4g12.4)') 
!     *          'ow ',inpos, el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!////////////////////////
            if(kcon .le. 4 .and. el > 0.d0) then
               length = el
               return
            endif
!             side wall on X axis
            normal = epPos(0.d0, 1.d0, 0.d0)
            dd = 0.
            call epxpLand4vp(OO, AA, BB, CC, xyz, dcos,
     *        el, kcon, jcon)
!////////////////////////
!               write(0, '(a,i3,1p,4g12.4)') 
!     *          'sw0 ',inpos, el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!////////////////////////

            if(kcon .le. 4 .and. el > 0 ) then
               length = el
               return
            endif
!            side wall between 1 and 2
            normal =epPos(-sint/(1.d0-cost)*sqrt((1.0+cost)/2),
     *              sqrt( (1.0+cost)/2.), 0.d0)
            dd = Xl*normal%x + Yl*normal%y
            
            call epxpLand4vp2(AAp, BBp,  CCp, OOp,
     *       xyz, dcos, normal, dd, el, kcon, jcon)
!////////////////////////
!               write(0, '(a,i3,1p,4g12.4)') 
!     *          'sw12 ',inpos, el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!////////////////////////
            if(kcon <= 4  .and. el > 0.d0) then
               length = el
               return
            endif

!            bottom wall
            normal = epPos(0.d0, 0.d0, -1.d0) 
            dd = 0.
            call epxpLand4vp2(OO, AA, AAp, OOp,
     *       xyz, dcos, normal, dd, el, kcon, jcon)
!////////////////////////
!               write(0, '(a,i3,1p,4g12.4)') 
!     *          'bot ',inpos, el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!////////////////////////
            if(kcon <= 4 .and. el > 0.d0) then
               length = el
               return
            endif
!            top wall
            normal = epPos(0.d0, 0.d0, 1.d0)
            dd = h
            call epxpLand4vp2(CC, BB, BBp, CCp,
     *       xyz, dcos, normal, dd, el, kcon, jcon)
!////////////////////////
!               write(0, '(a,i3,1p,4g12.4)') 
!     *          'top ',inpos, el,
!     *            el*dir.x+x+Xc,
!     *            el*dir.y+y+Yc,
!     *            el*dir.z+pos.z
!////////////////////////
            if(kcon <= 4  .and. el > 0.d0) then
               length = el
               return
            endif
         endif
         icon = 1
      endif
      end subroutine epbhoneycm1

      subroutine epshoneycomb(comp, pos, icon)
      use honeycomb
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "ZepDirec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
!           judge if a given 'pos' is inside 'comp'
!         
       type(Component)::  comp !input honeycomb component 
       type(epPos)::  pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside

       type(epPos)::  cpos   


      call epv2c_honeycomb(comp, pos, cpos)
      call epshoneycm0(comp, cpos, icon)
      end      subroutine epshoneycomb


      subroutine epshoneycm0(comp, pos, icon)
      use honeycomb
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "ZepDirec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
!           judge if a given 'pos' is inside 'comp'
!         
       type(Component)::  comp !input honeycomb component 
       type(epPos)::  pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside


      real(8)::x, y
      integer jcon

      call ep_honeycmCnst(comp)

      inpos = 0
      if( pos%z < 0. .or.  pos%z > h) then
         icon = 1
      elseif( pos%x < Xs .or. pos%x > Xs + a ) then
         icon = 1
      elseif( pos%y < Ys .or. pos%y > Ys + b ) then
         icon = 1
      else
         call ephoneyxy2unit(pos%x, pos%y,  jcon)
!                 pos.x,y is in a unit  honeycomb (honeycomb may not be
!                    closed one) 
         x = pos%x - Xc
         y = pos%y - Yc
         call ep_honeyin(x, y,  icon)
      endif
      end subroutine epshoneycm0

      subroutine ep_honeyin(x, y,  icon) 
      use honeycomb
      implicit none
      real(8),intent(in):: x
      real(8),intent(in):: y
      integer,intent(out):: icon
 
      icon = 1

      if(  abs(x) < Xl ) then
         if( abs(y) < Yl ) then
            return              ! in central inner square
         else
!            if( abs(y)>=Yl abs(y)<=Yu ) OK
            icon = 0
            if( y > 0 ) then
               inpos = 2
            else
               inpos = 5
            endif
         endif
      else  ! abs(x) >= Xl  
         if( abs(y) < f1i(abs(x)) ) then
!            in one of 4  inner triangles
            return
         elseif( abs(y) <= f1o(abs(x))) then
            icon = 0
            if( abs(y) <= g12(abs(x)) ) then
               if(x >= 0 ) then
                  if( y >= 0 ) then
                     inpos = 1
                  else
                     inpos = 6
                  endif
               else
                  if( y >= 0 ) then
                     inpos = 3
                  else
                     inpos = 4
                  endif
               endif
            else
!              small area near corner of horizontal bar 
               if(x > 0 ) then
                  if(y> 0 ) then
                     inpos = 2
                  else
                     inpos = 5
                  endif
               else
                  if(y> 0 ) then
                     inpos = 2
                  else
                     inpos = 5
                  endif
               endif
            endif
            return
         endif
      endif
      end subroutine  ep_honeyin

      subroutine ep_honeyin2(x, y,  icon) 
!      see if x,y is inner hexgon;  (x,y) in 
!          unit hexagon 
      use honeycomb
      implicit none
      real(8),intent(in):: x
      real(8),intent(in):: y
      integer,intent(out):: icon  ! 0--> inner hexagon else 1
      icon = 0
      if(  abs(x) < Xl ) then
         if( abs(y) < Yl ) then
            return              ! in central inner square
         endif
      else  ! abs(x) >= Xl  
         if( abs(y) < f1i(abs(x)) ) then
!            in one of 4  inner triangles
            return
         endif
      endif
      icon = 1
      end subroutine  ep_honeyin2

!     **************************************
      subroutine epenvlphoneycomb(comp, org, abc)
      use honeycomb
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


      call epenvlphoneycm0(comp, org, abc)
      call epc2v_honeycomb(comp, abc, abc)

      end    subroutine epenvlphoneycomb


      subroutine epenvlphoneycm0(comp, org, abc)
      use honeycomb
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



      real*8  zmin, zmax
      integer n
!
      call ep_honeycmCnst(comp)

      org%x = Xs
      org%y = Ys
      org%z = 0.
      abc%x = a
      abc%y = b
      abc%z = h
      NVTX = 0
      end  subroutine epenvlphoneycm0
!     _____________________________ 
      subroutine ep_honeycmCnst(comp)
!      set consts 
      use honeycomb
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

       type(Component)::  comp  ! input.   component.
      if(CompNum /= comp%vol ) then
         LL = Volat( comp%vol + 1)
         LI = Volat( comp%vol + 2)
         Ls = Volat( comp%vol + 3)
         h = Volat( comp%vol + 4)
         W = Volat( comp%vol + 5)
         Xs = Volat( comp%vol + 6)
         a =  Volat( comp%vol + 7)
         Ys = Volat( comp%vol + 8)
         b = Volat( comp%vol + 9)
         if(LL < 0 ) then
            if(abs(LL) == 60.) then
               cost=0.5d0
               sint=sqrt(1.d0-cost**2)
               tant=sint/cost
            else
               cost=cos(abs(LL)*Torad)
               sint=sin(abs(LL)*Torad)
               tant=sint/cost
            endif
         else
            cost = (LL-LS)/2/LI
            if(cost <= 0. .or. cost>=1.0d0) then
               write(0,*) ' LL,LI,LS for honeycm relataion invalid'
               write(0,*) ' they are =', LL, LI, LS
               write(0,*) ' cos angle is ',cost
               stop
            endif
            sint=sqrt(1.d0-cost**2)
            tant=sint/cost
         endif 
         LIsint = LI*sint
         LIcost = LI*cost
         HSy = LIsint + W  
         Sy = 2*HSy              ! y step
         HSx = LIcost + LS
         Sx = 2*HSx             ! xstep
         HLS = LS/2         
         offx  = LIcost + HLS
         HW = W/2

         HWsec = HW/cost
         Yu = LIsint + W
         Yl = LIsint 
!         x at which y=Yl for f1i:; i.e; f1i(x)=Yl
!       -tant*(Xl-HLS) + Yl + HW - HWsec=Yl
         Xl =HW*(cost-1.)/sint + HLS
         Xl0 =-(HWsec-HW)/tant + offx ! at y=0
!            = HW(cost-1)/sint + offx
!         x at which y=Yu for f1o, i.e, f1o(x) = Yu
!        -tant*(x-HLS)+ Yl+ HW + HWsec = Yu
         Xu=HW*(1.d0-cost)/sint + HLS
!         x at which y=0 for f1o i.e, f1o(x) = 0
         Xu0= (HWsec+HW )/tant   + offx
!           x at which y=W for f1o f1(o(x) = W
!         Xuyw = HW*(1.0 - cost)/sint + offx =
         Xuyw = Xu + LIcost
!          gradient of line (Xl,Yl)--(Xu,Yu);  (Yu-Yl)/(Xu-Xl)=
!          W/{(Yl+HWsec-HW-Yu+HWsec+Hw)/tant}
!          W/{-W +2HWsec}/tant
!          W/{W(sec-1)/tant} = (1/cost-1)cost/sint
!          (1-cost)/sint (if t=90, 1/1 = 1; ok)   
         grad = sint/(1.0-cost)   ! common to (xl0,0)-(xuyw,w) grad

!c//////////////
!         write(0,*) ' honeycom cnst'
!         write(0,*) ' LL, LI, LS, h, W=', LL,LI,LS,h, W
!         write(0,*) ' Xs, a=',Xs,a
!         write(0,*) ' Ys, b=',Ys,b
!         write(0,*) ' offx=',offx
!         write(0,*) ' Sx, Sy=',Sx, Sy
!         write(0,*) ' LIcost=',LIcost, ' LIsint=',LIsint
!         write(0,*) ' Yu, Yl, Xl=',Yu, Yl, Xl
!         write(0,*) ' Xu, Xu0=', Xu, Xu0
!c//////////////
!         Compnum = comp.cn
         Compnum = comp%vol
      endif

      end subroutine ep_honeycmCnst
!      program main
!          to test ephoneyxy2unit
!      implicit none
!      real(8):: x, y
!      real(8):: Xc, Yc
!      integer:: nx, ny
!      integer:: icon
!      
!
!      do while(.true.)
!         write(0,*) "Eneter x,y values (>=0)"
!         read(*,*) x, y
!         if(x< 0.) exit
!         write(0,*) ' x y =', x, y
!         call ephoneyxy2unit(x, y, nx, ny, Xc, Yc, icon)
!         write(*,'(1p,4g13.3,2i5, i2)') x, y, Xc, Yc, nx, ny, icon
!      enddo
!      end


      subroutine ephoneyxy2unit(x, y, icon)
!        find unit honeycomb figure in which  given (x,y) 
!     point is contained.
!     The center of the unit is (Xc, Yc).       
!     This routine does not consider Xs, Ys, a, b
      use honeycomb
      implicit none
      real(8),intent(in):: x, y  ! given input point in local coord.
      integer,intent(out):: icon  ! =0 ok. else  x,y outside
                              !  Xc, Yc.  are alos output
      integer nx, ny
      integer n

      real(8):: xp, yp
      
      ny = y/HSy
      yp = y - ny*HSy 
      n = x/HSx
      xp = x-n*HSx
      if( mod(ny,2) == 0 ) then
         if(mod(n,2) /= 0 ) then
            if(yp <= f1o(xp)) then
               nx=n
               Xc = nx*HSx
               Yc = ny*HSy
            else
               nx=n+1
               Xc=nx*HSx
               Yc=(ny+1)*HSy
            endif
         else
            yp = yp-HSy
            if(yp >=f6i(xp)) then
               nx=n
               Xc=nx*HSx
               Yc = (ny+1)*HSy
            else
               nx=n+1
               Xc=nx*HSx
               Yc=ny*HSy
            endif
         endif
      else
         if( mod(n,2) /= 0) then
            yp = yp - HSy
            if( yp >= f6i(xp) ) then
               nx= n
               Xc=nx*HSx
               Yc=(ny+1)*HSy
            else
               nx=n+1
               Xc=nx*HSx
               Yc=ny*HSy
            endif
         else
            if(yp<=f1o(xp)) then
               nx=n
               Xc=nx*HSx
               Yc=ny*HSy
            else
               nx=n+1
               Xc=nx*HSx
               Yc=(ny+1)*HSy
            endif
         endif
      endif
!/////////////
!      write(0, *) ' nx, ny=', nx, ny
!      write(0, *) ' Xc, Yc=', Xc, Yc
!///////////
      icon  = 0

      end  subroutine ephoneyxy2unit
      subroutine epatlochoneycomb(comp,loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
       type(Component)::  comp   ! input.
      integer loc(*)
      integer i
      do i = 1, 9
         loc(i) = i
      enddo
      end subroutine epatlochoneycomb

!=============================
      subroutine epv2c_honeycomb(comp, posv, posc)
      use honeycomb
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp  ! input.   component. 
       type(epPos)::  posv     ! input   
       type(epPos)::  posc     ! output    can be posv

       type(epPos)::  temp

      call ep_honeycmCnst(comp)
      if( comp%struc  == 'honeycomb' .or.
     *    comp%struc  == 'honeycomb_w'  ) then
         posc = posv
      elseif( comp%struc(1:12)  == 'honeycomb_xy' ) then
         posc = posv
      elseif( comp%struc(1:12) == "honeycomb_yz" ) then
         temp =epPos(posv%y, posv%z, posv%x)
         posc = temp
      elseif( comp%struc(1:12) == "honeycomb_zx" ) then
         temp = epPos( posv%z, posv%x, posv%y)
         posc = temp
!      elseif( comp.struc == "honeycomb_x" ) then
!         posc = posv
!      elseif( comp.struc == "honeycomb_y" ) then
!         posc = epPos( posv.y, b-posv.x, posv.z)
!      elseif( comp.struc == "honeycomb_z" ) then
!         posc = epPos( posv.z, b-posv.y, posv.x)
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epv2c_honeycomb is invalid'
         stop
      endif
      end   subroutine epv2c_honeycomb

      subroutine epc2v_honeycomb(comp, posc, posv)
      use honeycomb
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp  ! input.   component.
       type(epPos)::  posc     ! input 
       type(epPos)::  posv     ! output can be posc

       type(epPos)::  temp

      call ep_honeycmCnst(comp)
      if( comp%struc  == 'honeycomb' .or.
     *    comp%struc  == 'honeycomb_w'  ) then
         posv = posc
      elseif( comp%struc(1:12)  == 'honeycomb_xy' ) then
         posv = posc
      elseif( comp%struc(1:12) == "honeycomb_yz" ) then
!                yzx:  y<->x   z<->y  x<->z
         temp = epPos( posc%z, posc%x, posc%y)
         posv = temp
      elseif( comp%struc(1:12) == "honeycomb_zx" ) then
!            zxy 
         temp= epPos(posc%y, posc%z, posc%x)
         posv = temp
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epc2v_honeycomb is invalid'
         stop
      endif
      end   subroutine epc2v_honeycomb

      subroutine epv2cd_honeycomb(comp, dirv, dirc)
      use honeycomb
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
       type(Component)::  comp  ! input.   component.
       type(epDirec)::  dirv     ! input
       type(epDirec)::  dirc     ! output can be dirv

       type(epDirec)::  temp

      call ep_honeycmCnst(comp)
      if( comp%struc  == 'honeycomb' .or.
     *    comp%struc  == 'honeycomb_w'  ) then
         dirc = dirv
      elseif( comp%struc(1:12)  == 'honeycomb_xy' ) then
         dirc = dirv
      elseif( comp%struc(1:12) == "honeycomb_yz" ) then
!          yzx
         temp = epDirec(dirv%y, dirv%z, dirv%x)
         dirc = temp
      elseif( comp%struc(1:12) == "honeycomb_zx" ) then
         temp = epDirec(dirv%z, dirv%x, dirv%y)
         dirc = temp
!      elseif( comp.struc == "honeycomb_x" ) then
!         dirc = dirv
!      elseif( comp.struc == "honeycomb_y" ) then
!         dirc = epDirec( dirv.y, -dirv.x, dirv.z)
!      elseif( comp.struc == "honeycomb_z" ) then
!         dirc = epDirec( dirv.z, -dirv.y, dirv.x)
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epv2cd_honeycomb is invalid'
         stop
      endif
      end subroutine epv2cd_honeycomb

      subroutine epc2vd_honeycomb(comp, dirc, dirv)
      use honeycomb
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
       type(Component)::  comp  ! input.   component.
       type(epDirec)::  dirc     ! input
       type(epDirec)::  dirv     ! output can be dirc

       type(epDirec)::  temp
      call ep_honeycmCnst(comp)
!     (X,Y,Z)<---(x,y,z)
      if( comp%struc  == 'honeycomb' .or.
     *    comp%struc  == 'honeycomb_w' ) then
         dirv = dirc
      elseif( comp%struc(1:12)  == 'honeycomb_xy' ) then
         dirv = dirc
      elseif( comp%struc(1:12) == "honeycomb_yz" ) then
!         cdir = epDirec(dir.y, dir.z, dir.x)  X:y  Y:z  Z:x
!                 yzx
         temp = epDirec(dirc%z, dirc%x, dirc%y)
         dirv = temp
      elseif( comp%struc(1:12)  == 'honeycomb_zx' ) then
!         cdir = epDirec(dir.z, dir.x, dir.y)  X:z Y:x Z:y 
         temp = epDirec(dirc%y, dirc%z, dirc%x)
         dirv = temp 
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epc2vd_honeycomb is invalid'
         stop
      endif
      end   subroutine epc2vd_honeycomb




