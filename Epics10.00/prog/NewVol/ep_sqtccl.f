!          square to circle.
!    This uses approximation using triangles to express the wall from 
!    ccl to  sq.  All triagles'  vertex points are aligned counter-clock wise
!    and this leads to that even if the shape is completely symetric as to
!    the axis, there appears some assymetry.  (order of 0.1 mm for few cm
!    scale figures).
!     
#include "Zcondc.h"
      module sqtccl
!           to be stored as vol attrib.
      real(8),save:: a, b, h, r
      real(8),save:: ri   ! inscribing  ccl raius when the circle is approximated
                         ! by a polygon.
      real(8),save::Fai0           ! azim. angle to corner of the square (rad)
      real(8),save::Dfai1, Dfai2   !  unit opeining angle in 1st and 2nd region (rad)

!      --------------------------- computed every time
      integer,save::Nfai1, Nfai2   !  number of edges in first and 2nd region
                               !   (fai=0~90: deg)
      logical,save:: cylinderInside, boxInside
! ------------------------------------------------------
      integer,save::Compnum=-1000 
      real(8),parameter::eps= 1.d-8
      real(8),parameter::eps2= 1.d-5
      integer,save:: xconst, yconst
#if defined (MATHLOUSY)
      real(8),parameter:: hpi = 1.57079632679489661923d0
#else
      real(8),parameter:: hpi = asin(1.0d0)
#endif
      real(8),parameter:: pi =  hpi*2
      real(8),parameter:: twopi =  pi*2
      real(8),parameter:: Torad =  pi/180.d0
      integer,parameter::Npolygon=60   ! circle is approximated by npolygon
                                       ! must be multiple of 4
      integer,save::check=0
      integer,save::topbot

      end module sqtccl



!
!     see Fig/NewVol2.pdf 
!      
      subroutine eprsqtccl(comp)
      implicit none
#include "ZepManager.h"
#include "Zep3Vec.h"
#include "ZepPos.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "sqtccl"
!
       type(Component)::  comp   ! output. to recieve the config data.

!             next is mandatory
      call eprpst(comp, 4, 8, 1, 6)
      call ep_sqtcclCnst0(comp)
      if(MsgLevel >= 2 ) then
         call ep_sqtcclPrintCnst
      endif

      end subroutine eprsqtccl

!   ***************************************
      subroutine epbsqtccl(comp, pos, dir, length, icon)
!           find length to the boundary
      use sqtccl
      implicit none
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"
#include "Zepdebug.h"

!
!    find length to the boundary of a sqtccl specfied by 
!    'comp' from 'pos' with direction cos 'dir'
!    'pos' and 'dir' are given in this 'comp' local coordinate.
!     pos is may or may not be inside comp.

       type(Component):: comp    ! input.  sqtccl

       type(epPos)::  pos        ! input.  position.
       type(epDirec)::  dir      ! input. direction cosinse

      real(8),intent(out):: length    !  length cm from pos to the boundary
      integer,intent(out):: icon      !  0: length obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume

      real(8)::el
       type(epPos)::  cpos
       type(epDirec)::  cdir 

!   
      call epv2c_sqtccl(comp, pos, cpos)    ! for _y _x, get canonical pos
      call epv2cd_sqtccl(comp, dir, cdir)     ! and dir

      call epbsqtccl0(comp, cpos, cdir, length, icon)

      end subroutine epbsqtccl


      subroutine epbsqtccl0(comp, pos, dir, length, icon)
!           find length to the boundary
      use sqtccl
      implicit none
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "ZepPos.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
#include "Zepdebug.h"

!
!    find length to the boundary of a sqtccl specfied by 
!    'comp' from 'pos' with direction cos 'dir'
!    'pos' and 'dir' are given in canonical  'comp'  coordinate.
!     pos is may or may not be inside comp.

       type(Component):: comp    ! input.  sqtccl

       type(epPos)::  pos        ! input.  position.
       type(epDirec)::  dir      ! input. direction cosinse

      real(8),intent(out):: length    !  length cm from pos to the boundary
      integer,intent(out):: icon      !  0: length obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume

      real(8)::el, temp

       type(epPos)::  org, abc
      
      call epenvlpsqtccl0(comp, org, abc)
!           see bounding box (el is length to xing point, 
!           if any, else icon =-1
      call kxplbx( pos%x-org%x, pos%y-org%y, pos%z-org%z, 
     *               dir%x, dir%y, dir%z, 
     *               abc%x, abc%y, abc%z, el, icon)

      if(icon == -1)  return  ! no crossing 
      if(icon == 1)  then
!            pos is outside but may cross
         if( pos%z  < 0.) then
            if( abs(pos%z+ el*dir%z) < eps ) then
               temp =pos%x + el*dir%x 
               if( temp > -a .and. temp < a ) then
                  temp =pos%y + el*dir%y
                  if(temp > -b .and. temp< b) then
!                    xp is at bottom
                     length = el ! xing with the bottom
                     return
                  endif
               endif
            endif
         endif
      endif

!          icon = 1 may come here
      if( icon == 0 ) then   ! inside box so
!          see if pos is inside. next is exact (polygon approx. of ccl
!          is consideed)
         call epssqtccl0(comp, pos, icon)

         topbot = 0  !  used delicate case where wall may intercept 
                     !  inner point and ceil or floor.
         if(icon ==  0 ) then
!           inside.  so there should   be x-ing point;  so icon =0 
            call epbsqtcclin(comp, pos, dir, length)
         else
             !  pos is  outside
            call epbsqtcclot(comp, pos, dir, length, icon)
         endif
      else
!            pos is outside
         call epbsqtcclot(comp, pos, dir, length, icon)
      endif
      end       subroutine epbsqtccl0 

      subroutine epbsqtcclin(comp, pos, dir, el)
!         get length to the boundary when pos is  inside
!       of sqtccl
      use sqtccl
      implicit none
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component):: comp    ! input.  sc2ccl

       type(epPos)::  pos
       type(epDirec)::  dir
      real(8),intent(out):: el  ! length to the boundary

      real(8)::x, y
      integer::icon, jcon, kcon
       type(epPos)::  posh
       type(epPos)::  p1, p2, p3, p4
      integer::q
      real(8):: el1, el2, es1, es2, temp, elsave
      real(8):: fai1, fai2, faiin, fai, Dfai
      logical::delicate, posIsInCcl, posIsInSq
      integer inout


      posIsInSq  = -a <= pos%x .and. pos%x<=a
     *       .and. -b <= pos%y .and. pos%y<=b .and. boxInside
      inout = -1
      delicate = .false.
      if(dir%z /= 0.0) then
!          see bottom
         el = -pos%z/dir%z
         if(el > 0.) then
            x = pos%x + el*dir%x
            y = pos%y + el*dir%y
            if(x >= -a .and. x<= a) then
               if( y >= -b .and. y<= b) then
                  if( posIsInSq ) then
                     return
                  endif
!                      before reaching bottom, wall may intercept. so must see
!                      
                  topbot = 2
                  elsave = el
                  inout = 3
               endif
            endif
         endif
!           see top
         posIsInCcl = pos%x**2 + pos%y**2 <= ri**2 .and. cylinderInside

         el = (h-pos%z)/dir%z
         if( el > 0.)  then
            elsave = el
            x = pos%x + el*dir%x
            y = pos%y + el*dir%y
            if( x**2 + y**2  <= ri**2) then  ! inside inner circle
               if( posIsInCcl) then
                  return
               endif
               topbot=1 !  must see if wall may intercept
               inout = 4
            elseif(x**2 + y**2  <= r**2) then

               !  we should see if inside the tangential polygon
               faiin = atan2(y, x)
               call epsqtcclStepFai(faiin, fai, Dfai)
               call epsqtcclGet4p(fai, Dfai, p1, p2, p3, p4, q)
               call k3inout(x, y,
     *              0.d0, 0.d0, p1%x, p1%y,  p4%x, p4%y, inout)
!               if(icon == 0 ) then
!                  return
!               endif
               if(inout == 0 ) then
                  if( posIsInCcl ) then
                     return
                  endif
                  topbot = 1
               endif
               delicate = posIsInCcl
!               if(inout == 1 .and. icon == 1) then
!                  write(0,*)' pos and triagle'
!                  write(0,*) ' x,y=', pos.x, pos.y
!                  write(0,*) ' 1st p=', p1.x, p1.y
!                  write(0,*) ' 2nd p=', p4.x, p4.y
!               endif
            else
               inout = 2
            endif
         endif
      else
         inout = 5
      endif
!         x-ing point should be  on the side wall
!         get fai region to inspect
!            get xp with ccl  on the projection
      call kxplccl(pos%x, pos%y, dir%x, dir%y, r, el1, el2, jcon)
          !     jcon = 0 pos is inside ccl. el1 is usable
          !            1 pos is outside el1 and el2 are xp
          !           -1 no xp 
          ! see xp with square
      call kxplsq(pos%x+a, pos%y+b, dir%x, dir%y,
     *      2.d0*a, 2.d0*b, es1, es2, kcon)
          !  kcon: same meaning as jcon above.

      fai1 = atan2(pos%y, pos%x)
      if( jcon ==  0 .and. kcon == 0 ) then
         ! xpoint should be between  pos ~ max(es1,el1) 
         if( es1 > el1 ) then
            fai2 = atan2(pos%y+es1*dir%y, pos%x + es1*dir%x)
         else
            fai2 = atan2(pos%y+el1*dir%y, pos%x + el1*dir%x) 
         endif
       elseif( jcon == 1 .and.  kcon == 0 ) then
          if( es1 >  el2 ) then
             fai2 = atan2(pos%y+es1*dir%y, pos%x + es1*dir%x)
          else
             fai2 = atan2(pos%y+el2*dir%y, pos%x + el2*dir%x) 
          endif
       elseif( jcon == 0 .and. kcon == 1 ) then
          ! ccl is larger than square; es2 ~ el1
          if( es2 > el1 ) then
             fai2 = atan2(pos%y+es2*dir%y, pos%x + es2*dir%x)
          else
             fai2 = atan2(pos%y+el1*dir%y, pos%x + el1*dir%x) 
          endif
       elseif( jcon == -1 .and. kcon == 0 ) then
          !  pos~ es1
          fai2 = atan2(pos%y+es1*dir%y, pos%x + es1*dir%x)
       elseif( jcon == 0 .and. kcon == -1) then
          ! pos ~ el1
          fai2 = atan2(pos%y+el1*dir%y, pos%x + el1*dir%x)
       elseif( jcon == 1 .and. kcon == -1 ) then
          fai2 = atan2(pos%y+el2*dir%y, pos%x + el2*dir%x) 
       elseif( jcon == -1 .and. kcon == 1 ) then
          fai2 = atan2(pos%y+es2*dir%y, pos%x + es2*dir%x)
       elseif( jcon == 1 .and. kcon == 1 ) then  
          if(es2 > el2) then
             fai2 = 
     *            atan2(pos%y+es2*dir%y, pos%x + es2*dir%x)
          else
             fai2 = 
     *            atan2(pos%y+el2*dir%y, pos%x + el2*dir%x) 
          endif
       else
!     can happen 
          call epsqtcclStepFai(fai1, fai, Dfai)
          fai2 = fai1 - Dfai
          fai1 = fai1 + Dfai
!          ! impossilbe
!          write(0,*) ' strange: inside but outside ?'
!          write(0,*) 'error in ep_sqtccl Cn=',comp.cn
!          write(0,*) ' jcon,kcon =',jcon, kcon
!          write(0,*)' cpos=',pos
!          posh = epPos(pos.x+dir.x*15.d0, pos.y+dir.y*15.d0, 
!     *         pos.z + dir.z*15.d0)
!          write(0,*) ' cpos2=',posh
!          write(0,*) 'c dir=',dir
!
!         call epc2v_sqtccl(comp, pos, pos)
!         call epl2w(comp.cn,  pos, pos)
!         call epc2v_sqtccl(comp, posh, posh)
!         call epl2w(comp.cn,  posh, posh)
!         call epc2vd_sqtccl(comp, dir, dir)
!         call epl2wd(comp.cn, dir, dir)
!         write(0,*) ' world pos  ='
!         write(0,*)  pos
!         write(0,*) ' world pos2 ='
!         write(0,*)  posh
!         write(0,*) ' world dir ='
!         write(0,*)  dir
!
!         write(0,*) ' a, b=', a, b
!         write(0,*) ' es1,2 =', es1, es2
!         stop
       endif

       if( delicate ) then
          fai1 = atan2(y, x)
          fai2 = fai1
       endif

       if(fai1 < 0.) fai1 = twopi + fai1
       if(fai2 < 0.) fai2 = twopi + fai2
!         now fai is 0 to twopi.
       temp = fai1
       if( fai1 > fai2 ) then
          fai1 = fai2
          fai2 = temp
       endif
       if(  (fai2 - fai1) > pi ) then
          temp = fai1
          fai1 = fai2
          fai2 = temp
       endif
!       examine from fai1 to fai2 
!          fai1 may be > fai2 (say fai1 = 270*Torad
!          fai2 = 10*Torad. )
          ! put some margin
!       fai1 = fai1 - Dfai1       ! Dfai1 or 2 may be ok
!       fai2 = fai2 + Dfai1
!                 crossing point should be in fai1 to fai2
      call epxpsqtcclW(pos, dir, 0, fai1, fai2, el, icon)
      if(icon /= 0) then
         if( topbot /=  0 ) then
            el = elsave   ! ceil or floor xp is correct
            return
         endif

         write(0,*)
     *    ' No crossing is strange in ep_sqtccl;cn= ', comp%cn 
         write(0,*) ' struc=',comp%struc
         write(0,*) 'icon=',icon,' jcon=', jcon, ' kcon =',kcon
         write(0,*) ' inout =', inout, ' topbot=',topbot

         write(0,'(a, 1p,3g16.8)') ' canonical pos  =', pos
         posh =
     *    epPos(pos%x+15*dir%x, pos%y+15*dir%y,  pos%z+15*dir%z)
         write(0,'(a, 1p,3g16.8)') ' canonical pos2 =', posh
         write(0,'(a,1p,3g16.7)') ' dir =', dir
         write(0, *)  ' fai(deg)=', atan2(dir%y, dir%x)/Torad
         call epc2v_sqtccl(comp, pos, pos)
         call epl2w(comp%cn,  pos, pos)
         call epc2v_sqtccl(comp, posh, posh)
         call epl2w(comp%cn,  posh, posh)
         call epc2vd_sqtccl(comp, dir, dir)
         call epl2wd(comp%cn, dir, dir)
         write(0,*) ' world pos  ='
         write(0,*)  pos
         write(0,*) ' world pos2 ='
         write(0,*)  posh
         write(0,*) ' world dir ='
         write(0,*)  dir
         write(0,*)  ' fai=', atan2(dir%y, dir%x)/Torad
         write(0,*) ' canonical fai region fai1 =', fai1/Torad,
     *    ' fai2 =',fai2/Torad
!         write(0,*) '-------------------'
!         call epfordebug('from sqTcll')
!         write(0,*) '-------------------'
         stop
      endif
      end subroutine epbsqtcclin


      subroutine epbsqtcclot(comp, pos, dir, length, icon)
!              pos.z > 0  and pos is outside of
!              comp
      use sqtccl
      implicit none
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component):: comp    ! input.  sc2ccl
       type(epPos)::  pos
       type(epDirec)::  dir
      real(8),intent(out)::length
      integer,intent(out)::icon  ! 1 if cross -1 if no

       type(epPos)::  posh
      real(8)::x, y
      integer:: face
      real(8)::el, fai1, fai2      
      real(8)::ela(2)      
      integer::icona(2)
      integer:: xpc, jcon

      xpc = 0
!         first check if x with bounding horse
      posh = epPos( pos%x + a, pos%y + b, pos%z)
      call kxplhorse(2.0d0*a, 2.0d0*b, h, a-r, b-r, 
     *  2.d0*r, 2.d0*r, posh, dir,  el, icon, face)
      if(icon == -1) then
         return                 !************  no crossing
      endif
!             from top side.
      if(face == 6 .and.  dir%z <  0.0 .and. pos%z > h ) then
         el =  (h-pos%z)/dir%z
         x = pos%x + el*dir%x
         y = pos%y + el*dir%y
         if(x**2 + y**2 < ri**2 ) then
            icon = 1
            length = el
            return   ! ******************
         elseif(x**2 + y**2 < r**2) then
              !  within top ccl. but  we should see if inside the polygon
              !  use epssqtccl0  (prob. small)
            posh = epPos(x, y, h)
            call epssqtccl0(comp, posh, icon)
            if(icon == 0 ) then
               icon = 1
               length = el
               return
            endif
         endif
      endif
!          see details with side walls
      if( face == 2 ) then
         fai1 = 0.
         fai2 = hpi
         call epxpsqtcclW(pos, dir, 1, fai1, fai2, length, icon)      
         if( icon == 0 ) then
            icon = 1
            return  !************
         endif
         fai1 = 3*hpi
         fai2 = twopi
         call epxpsqtcclW(pos, dir, 1, fai1, fai2, length, icon)      
         if( icon == 0 ) then
            icon = 1
            return
         endif
      elseif( face == 4 ) then
         fai1 = 0.
         fai2 = pi
         call epxpsqtcclW(pos,dir, 1, fai1, fai2, length, icon)      
         if( icon == 0 ) then
            icon = 1
            return
         endif
      elseif( face == 5 ) then
         fai1 = hpi
         fai2 = 3*hpi
         call epxpsqtcclW(pos,dir, 1, fai1, fai2, length, icon)      
         if( icon == 0 ) then
            icon = 1
            return
         endif
      elseif( face == 3 ) then
         fai1 = pi
         fai2 = twopi
         call epxpsqtcclW(pos,dir, 1, fai1, fai2, length, icon)      
         if( icon == 0 ) then
            icon = 1
            return
         endif
      endif
      icon = -1
      end    subroutine epbsqtcclot


!   ***************************************
      subroutine epsqtcclStepFai(faiin, fai, Dfai)
      use sqtccl
      real(8),intent(in):: faiin ! may be < 0, >twopi
      real(8),intent(out):: fai  ! faiin belongs
      real(8),intent(out):: dfai !     (fai, fai+dfai)  region

      real(8):: faix
      integer:: j
      faix = faiin
      if(faix < 0.d0 ) faix = faix+twopi

      faix = mod(faix, twopi)

      if(faix < Fai0 ) then
         j = faix/Dfai1
         fai =j*Dfai1
         dfai = Dfai1
      elseif( faix < pi-Fai0 ) then
         j = ( faix-Fai0)/Dfai2
         fai = j*Dfai2 + Fai0
         dfai = Dfai2
      elseif( faix < pi+Fai0) then
         j = ( faix-(pi-Fai0))/Dfai1
         fai = j*Dfai1 + (pi-Fai0)
         dfai = Dfai1
      elseif( faix < twopi-Fai0) then
         j = ( faix-(pi+Fai0) )/Dfai2
         fai = j*Dfai2 + (pi+Fai0)
         dfai = Dfai2
      else
         j = ( faix-(twopi-Fai0))/Dfai1
         fai = j*Dfai1 + (twopi-Fai0)
         dfai = Dfai1
      endif
      end    subroutine epsqtcclStepFai

      subroutine epsqtcclGet4p(faiin, Dfai, p1, p2, p3, p4, q)
      use sqtccl
      implicit none
#include "ZepPos.h"
!       for fai region (faiin, faiin+Dfai),  get 4 points
!       so that p1,p2,p3 be a triangle with p1 on the
!       ccl, and p2, p3 on the square. It is closer
!       to the x axis than the other triangle which is
!       formed by p1,p3,p4 where p4 is on the ccl.
!       the triangle is 1st or 3rd quadrant,
!       p1,p2,p3 and  p1,p3,p4 are counter clock-wise 
!       when seen   outside.  Otherwise, clock-wise
 
      real(8),intent(in):: faiin  ! >=0  must be stepped one
      real(8),intent(in):: Dfai
       type(epPos):: p1, p2, p3, p4
      integer,intent(out)::q  ! on of 1, 2,3,4 which quadrant  
      
      real(8):: fai, temp
      fai = mod(faiin, twopi)
      p1 = epPos(r*cos(fai), r*sin(fai), h)
      p4 = epPos(r*cos(fai+Dfai), r*sin(fai+Dfai), h)
      if( fai < Fai0-eps2 ) then
         q=1
         p2 = epPos(a, a*tan(fai), 0.d0)
         p3 = epPos(a, a*tan(fai+Dfai), 0.d0)
      elseif( fai < pi-Fai0-eps2 ) then
         q =2
         p2 = epPos(b*tan(hpi-fai), b, 0.d0)
         p3 = epPos(b*tan(hpi-fai-Dfai), b, 0.d0)
         if( fai < hpi-eps2) q = 1
      elseif( fai < pi+Fai0 -eps2 ) then
         q =3
         p2 = epPos(-a,  -a*tan(fai), 0.d0)
         p3 = epPos(-a,  -a*tan(fai+Dfai), 0.d0)
         if(fai < pi-Fai0-eps2) q = 2
      elseif( fai < twopi-Fai0-eps2 ) then
         q =4
         p2 = epPos(-b*tan(hpi-fai), -b, 0.d0)
         p3 = epPOs(-b*tan(hpi-fai-Dfai), -b, 0.d0)
         if( fai < 3*hpi-eps2) q = 3
      else
         q = 4
         p2 = epPos(a, a*tan(fai), 0.d0)
         p3 = epPos(a, a*tan(fai+Dfai), 0.d0)
      endif
      end subroutine epsqtcclGet4p

      function epsqtcclDfai(faiin) result(Dfai)
!          fai is classified  by Dfai1 or Dfai2
!        is to be used.  Also, xconst and yconst
!        is set to be -1,1 or 0 dending on
!        whether square edge is const in x or y
!        in that fai region.  0 means not const.
!         
      use sqtccl
      implicit none
      real(8),intent(in):: faiin
      real(8):: Dfai, fai

      if( faiin < 0.) then
         fai=twopi + faiin
      else
         fai = faiin
      endif
      fai = mod(fai, twopi)
      
      if( fai < Fai0 - eps2 )  then
         Dfai = Dfai1
         yconst = 0
         xconst = 1
      elseif( fai < pi-Fai0-eps2 ) then
         Dfai = Dfai2
         yconst = 1
         xconst = 0
      elseif( fai < pi+Fai0-eps2 ) then
         Dfai = Dfai1
         yconst = 0
         xconst = -1
      elseif( fai < twopi-Fai0-eps2 ) then
         Dfai = Dfai2
         yconst = -1
         xconst = 0
      else
         Dfai = Dfai1
         yconst = 0
         xconst = 1
      endif
      end function epsqtcclDfai


      subroutine epxpsqtcclW(pos, dir, io, ifai1, ifai2, length, icon)
!       exmine crossing point within fai1 and fai2
!
      use sqtccl
      implicit none
#include "ZepPos.h"
#include "ZepDirec.h"

       type(epPos)::  pos
       type(epDirec)::  dir
      real(8),intent(in)::ifai1, ifai2
      integer,intent(in):: io !  0 pos is inside 1 pos is outside
      real(8),intent(out)::length
      integer,intent(out)::icon  ! 0 if crossing point found
                                 ! -1 no crossing
 
      real(8):: Dfai, fai
      real(8):: fai1, fai2
      real(8)::  epsqtcclDfai

      integer:: xpc,  jcon
      real(8):: ela(2), temp, x, el
       type(epPos)::  p1, p2, p3, p4
      real(8):: towhich
      integer:: q


      fai1 = ifai1
      fai2 = ifai2
      xpc = 0
      Dfai = epsqtcclDfai(fai1)
      fai1 = fai1 - Dfai
      
      if(fai1 < 0.) then
         fai1 = twopi + fai1
      endif
      fai2 = fai2 + epsqtcclDfai(fai2)
      if( fai1 > fai2) then
         fai2 = fai2 + twopi
      endif
!        find starting  fai which is integral
!        of Dfai1 and/or Dfai2 and closest
!        to fai1 but <=fai1  
      call  epsqtcclStepFai(fai1, fai, Dfai)


      
      do while( fai < fai2+eps)
         call epsqtcclGet4p(fai, Dfai, p1, p2, p3, p4,q)
         call  epxpLand3vpB(
     *    p1, p2, p3, pos, dir, el,  towhich, icon)
         if(icon <= 3  .and. el > 0 ) then
            icon = 0
            xpc = xpc + 1
            ela(xpc) = el
            if(xpc == 2 .or. io == 0) exit
!            if( io == 1 ) then   
            if( towhich < 0.d0) exit   !!!
               ! should from out side to  inside
               ! and nearer xp than possible the ohter one
               ! so we need not wait xpc=2 and this is xp
         endif
         call epxpLand3vpB(
     *   p1, p3, p4, pos, dir, el, towhich, icon)
         if( icon <= 3 .and. el> 0. ) then
            icon = 0
            xpc = xpc + 1
            ela(xpc) = el
            if(xpc == 2 .or. io == 0 ) exit
            if( towhich < 0.d0 ) exit   !//////
         endif
         Dfai =  epsqtcclDfai(fai)
         fai = fai + Dfai
         Dfai =  epsqtcclDfai(fai)
      enddo
      if( xpc == 2 ) then
         length = min(ela(1), ela(2))
         icon =0
      elseif(xpc == 1) then
         length = ela(1)
         icon =0  
      else
         icon = -1
      endif
      end subroutine epxpsqtcclW

      subroutine epssqtccl(comp, pos, icon)
!
!           judge if a given 'pos' is inside 'comp'
!         
      use sqtccl
      implicit none

#include "Zep3Vec.h"
#include "ZepDirec.h"
#include "ZepPos.h"
#include "Zcnfig.h"

       type(Component)::  comp !input sqtccl component 
       type(epPos)::  pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside

       type(epPos)::  cpos
      call epv2c_sqtccl(comp, pos, cpos)
      call epssqtccl0(comp, cpos, icon)
      end      subroutine epssqtccl

      subroutine epssqtccl0(comp, pos, icon)
!
!           judge if a given 'pos' is inside 'comp'
!         
      use sqtccl
      implicit none

#include "Zep3Vec.h"
#include "ZepDirec.h"
#include "ZepPos.h"
#include "Zcnfig.h"

       type(Component)::  comp !input sqtccl component 
       type(epPos)::  pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside

       type(epPos)::  org, abc

      if(pos%z > h ) then
         icon = 1 
         return
      elseif(pos%z < 0.d0) then
         icon = 1
         return
      endif

      call epenvlpsqtccl0(comp, org, abc)
      call kseeinbox(org%x, org%x + abc%x,
     *               org%y, org%y + abc%y,
     *               org%z, org%z + abc%z,
     *     pos%x, pos%y, pos%z, icon)
      !  if not in the bounding box, outside
      if(icon == 1 ) return
      !  further check
      call  epssqtccl1(comp, pos, icon)
      end subroutine epssqtccl0


      subroutine epssqtccl1(comp, pos, icon)
!       see actually inside
      use sqtccl
      implicit none

#include "Zep3Vec.h"
#include "ZepDirec.h"
#include "ZepPos.h"
#include "Zcnfig.h"

!
!           judge if a given 'pos' is inside 'comp'
!         
       type(Component)::  comp !input sqtccl component 
       type(epPos)::  pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside


       type(epPos)::  p1, p2, p3, p4
       type(epPos)::  p14, p32, vprod
       type(ep3Vec)::  nm 
      real(8):: k,  faiij, fai, temp, fais, Dfai
      real(8):: diff1,  diff2
      integer:: j, q


      if( cylinderInside ) then
         call kseeincyl(0.d0, 0.d0, 0.d0, ri, h, 
     *        pos%x, pos%y, pos%z, icon)
        ! if inside of inner inscribing cyl, then inside
         if(icon == 0 ) return
      elseif( boxInside) then
         call kseeinbox(-a, a, -b, b, 0.d0, h,
     *      pos%x, pos%y, pos%z, icon)
        ! if inside of inner box, then inside
         if(icon == 0 ) return
      endif

      fai = atan2(pos%y, pos%x)
      call epsqtcclStepFai(fai,  fais, Dfai)
      call epsqtcclGet4p(fais, Dfai, p1, p2, p3, p4,q)
      call epwhichside0(p1, p2, p3, pos, nm, k, diff1)
      call epwhichside0(p1, p3, p4, pos, nm, k, diff2)
      if(diff1  < 0.0 .and. diff2 < 0. ) then
!                inside
         icon = 0
         return
      elseif( diff1 > 0.0 .and. diff2 > 0.) then
         icon = 1
         return
      endif
!       very delicate point  which is very near to the surface
      p14%x = p4%x - p1%x
      p14%y = p4%y - p1%y
      p14%z = 0.d0
      p32%x = p2%x - p3%x
      p32%y = p2%y - p3%y
      p32%z = 0.d0
      call epvectorProd(p14, p32, vprod)
      if( vprod%z < 0.) then 
         icon = 0
         return
      else
         icon = 1
         return
      endif

      end subroutine epssqtccl1


!     **************************************
      subroutine epenvlpsqtccl(comp, org, abc)
      use sqtccl
      implicit none

#include "Zep3Vec.h"
#include "ZepPos.h"
#include "Zcnfig.h"


!
!        give the envloping box of the component
!
       type(Component)::  comp  ! input.   component.
       type(epPos)::  org       ! output.  origin of the enveloping box
                               !          in local coord. 
       type(ep3Vec)::  abc      ! output.  a,b,c of the box


      call ep_sqtcclCnst(comp)
      call epenvlpsqtccl0(comp, org, abc)
      call  epenvlpsqtccl1(comp, org, abc)
      end      subroutine epenvlpsqtccl

!     *******************************
      subroutine epenvlpsqtccl0(comp, org, abc)
      use sqtccl
      implicit none

#include "Zep3Vec.h"
#include "ZepPos.h"
#include "Zcnfig.h"
!             note: ..ccl, ..ccl0 ..ccl1  are
!                ..ccl0 : canonical
!                ..ccl1 : to local
!          different from epbsqtccl, ccl0, ccl1
!                   ccl0:  local   ccl1: canonical
!        give the envloping box of the component
!
       type(Component)::  comp  ! input.   component.
       type(epPos)::  org       ! output.  origin of the enveloping box
                               !          in canonical coord. 
       type(ep3Vec)::  abc      ! output.  a,b,c of the box

!

      org%x =min( -a, -r)
      org%y =min( -b, -r)
      org%z = 0.
      abc%x = abs(org%x*2)
      abc%y = abs(org%y*2)
      abc%z = h
      NVTX = 0

      end  subroutine epenvlpsqtccl0

      subroutine epenvlpsqtccl1(comp, org, abc)
      use sqtccl
      implicit none

#include "Zep3Vec.h"
#include "ZepPos.h"
#include "Zcnfig.h"
!
!        give the envloping box of the component
!
       type(Component)::  comp  ! input.   component.
       type(epPos)::  org       ! in/output.  origin of the enveloping box
                               !          in local coord. 
       type(ep3Vec)::  abc      ! in/output.  a,b,c of the box
      real(8):: temp


      if(comp%struc == 'sqtccl_y') then
         !  y<->z
         temp = abc%y
         abc = ep3Vec(abc%x, abc%z, temp)
         org =epPos(org%x, 0.d0, -max(a,r))
      elseif(comp%struc == 'sqtccl_x') then
         !  z<-->x
         temp = abc%x
         abc = ep3Vec(abc%z, abc%y, temp)
         org =epPos(0.d0, org%y, -max(a,r))
      endif
!                This dose not work.  do as above
!      call epc2v_sqtccl(comp, org, org)

      end  subroutine epenvlpsqtccl1


!     _____________________________ 
      subroutine ep_sqtcclCnst0(comp)
!      set consts 
      use sqtccl
      implicit none

#include "Zep3Vec.h"
#include "ZepPos.h"
#include "Zcnfig.h"
      
       type(Component)::  comp  ! input.   component.

      real(8):: temp
!           check
         a = Volat( comp%vol + 1)
         b = Volat( comp%vol + 2)
         h = Volat( comp%vol + 3) 
         r = Volat( comp%vol + 4) 
         Fai0 = atan2(b,a)
         Nfai1 = Npolygon/4*b/(a+b) + 0.5d0
         Nfai2 = Npolygon/4- Nfai1

         Dfai1 = Fai0/Nfai1
         Dfai2 = (90.0*Torad-fai0)/Nfai2
!           radius of the inscribing circle (smallest one)
         temp = Dfai1
         if(temp < Dfai2) then
            temp = Dfai2
         endif
         ri = r *cos(temp/2.d0)    ! inscribing ccl radius
         Volat(comp%vol + 5)  = ri
         Volat(comp%vol + 6)  = Fai0
         Volat(comp%vol + 7)  = Dfai1
         Volat(comp%vol + 8)  = Dfai2
      end  subroutine ep_sqtcclCnst0
      subroutine ep_sqtcclCnst(comp)
!      set consts 
      use sqtccl
      implicit none

#include "Zep3Vec.h"
#include "ZepPos.h"
#include "Zcnfig.h"

       type(Component)::  comp  ! input.   component.

!      if(comp.cn /= Compnum ) then  !  honeycomb
      if(comp%vol /= Compnum ) then
         a = Volat( comp%vol + 1)
         b = Volat( comp%vol + 2)
         h = Volat( comp%vol + 3) 
         r = Volat( comp%vol + 4) 
         ri = Volat( comp%vol + 5) 
         Fai0 =Volat(comp%vol + 6)  
         Dfai1 = Volat(comp%vol + 7) 
         Dfai2 = Volat(comp%vol + 8)  
         if( r <= min(a,b)) then
            cylinderInside = .true.
            boxInside = .false.
         elseif( r**2 >= a**2 + b**2 ) then
            boxInside = .true.
            cylinderInside = .false.
         else
            boxInside = .false.
            cylinderInside = .false.
         endif
!         if( a<=0 .or. b<=0 .or. r > a .or. r > b ) then
!            write(0,*) ' a, b, r or h of sqtccln are invalid'
!            write(0,*) a, b, r, h
!            stop
!         endif
         Nfai1 = Npolygon/4*a/(a+b) + 0.5d0
         Nfai2 = Npolygon/4- Nfai1
!         Compnum = comp.cn
         Compnum = comp%vol
      endif
      end  subroutine ep_sqtcclCnst

      subroutine ep_sqtcclPrintCnst
      use sqtccl
      implicit none

#include "Zep3Vec.h"
#include "ZepPos.h"
#include "Zcnfig.h"
      write(0,*) ' sqtccl cnst'
      write(0,*) ' a, b, r, h=', a, b, r, h
      write(0,*) ' ri, Fai0=', ri, Fai0/Torad
      write(0,*) ' Dfai1, Dfai2 =',Dfai1/Torad, Dfai2/Torad
      end subroutine ep_sqtcclPrintCnst


      subroutine epatlocsqtccl(comp,loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
       type(Component)::  comp   ! input.
      integer loc(*)
      integer i
      do i = 1, 4
         loc(i) = i
      enddo
      end subroutine epatlocsqtccl


      subroutine epv2c_sqtccl(comp, posv, posc)
      use sqtccl
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp  ! input.   component.
       type(epPos)::  posv     ! input
       type(epPos)::  posc     ! output

      call ep_sqtcclCnst(comp)
      if( comp%struc  == 'sqtccl' .or.
     *    comp%struc  == 'sqtccl_w') then
         posc = posv
      elseif( comp%struc(1:8)  == 'sqtccl_z' ) then
         posc = posv
      elseif( comp%struc(1:8)  == 'sqtccl_y' ) then
         posc = epPos( posv%x, -posv%z, posv%y)
      elseif( comp%struc(1:8)  == 'sqtccl_x' ) then
         posc = epPos( -posv%z, posv%y, posv%x)
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbsqtccl is invalid'
         stop
      endif
      end   subroutine epv2c_sqtccl

      subroutine epc2v_sqtccl(comp, posc, posv)
      use sqtccl
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp  ! input.   component.
       type(epPos)::  posc     ! input
       type(epPos)::  posv     ! output

      call ep_sqtcclCnst(comp)
!           (X,Y,Z)--->(x,y,z)
      if( comp%struc  == 'sqtccl' .or.
     *    comp%struc  == 'sqtccl_w' ) then
         posv = posc
      elseif( comp%struc(1:8)  == 'sqtccl_z' ) then
         posv = posc
      elseif( comp%struc(1:8)  == 'sqtccl_y' ) then
         posv = epPos( posc%x, posc%z,  -posc%y)
      elseif( comp%struc(1:8)  == 'sqtccl_x' ) then
         posv = epPos( posc%z, posc%y, -posc%x)
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbsqtccl is invalid'
         stop
      endif
      end   subroutine epc2v_sqtccl

      subroutine epv2cd_sqtccl(comp, dirv, dirc)
      use sqtccl
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
       type(Component)::  comp  ! input.   component.
       type(epDirec)::  dirv     ! input
       type(epDirec)::  dirc     ! output

      call ep_sqtcclCnst(comp)

      if( comp%struc  == 'sqtccl' .or. 
     *    comp%struc  == 'sqtccl_w' ) then
         dirc = dirv
      elseif( comp%struc(1:8)  == 'sqtccl_z' ) then
         dirc = dirv
      elseif( comp%struc(1:8)  == 'sqtccl_y' ) then
         dirc = epDirec( dirv%x, -dirv%z, dirv%y)
      elseif( comp%struc(1:8)  == 'sqtccl_x' ) then
         dirc = epDirec( -dirv%z,dirv%y, dirv%x)
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbsqtccl is invalid'
         stop
      endif
      end      subroutine epv2cd_sqtccl

      subroutine epc2vd_sqtccl(comp, dirc, dirv)
      use sqtccl
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
       type(Component)::  comp  ! input.   component.
       type(epDirec)::  dirc     ! input
       type(epDirec)::  dirv     ! output

      call ep_sqtcclCnst(comp)
!           (X,Y,Z)<---(x,y,z)
      if( comp%struc  == 'sqtccl' .or.
     *    comp%struc  == 'sqtccl_w' ) then
         dirv = dirc
      elseif( comp%struc(1:8)  == 'sqtccl_z' ) then
         dirv = dirc
      elseif( comp%struc(1:8)  == 'sqtccl_y' ) then
         dirv = epDirec( dirc%x, dirc%z,  -dirc%y)
      elseif( comp%struc(1:8)  == 'sqtccl_x' ) then
         dirv = epDirec( dirc%z, dirc%y, -dirc%x)
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbsqtccl is invalid'
         stop
      endif
      end   subroutine epc2vd_sqtccl


