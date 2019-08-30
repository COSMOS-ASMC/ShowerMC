!         fpolygon:  flat polygon.
!      
      module fpolygon
      real(8),save:: height
      real(8),save:: xmin, xmax, ymin, ymax
      real(8),parameter::bigdist=1.0d20
      real(8),save:: elmin
      integer,save::Compnum=-1000
      integer,save:: npoly
      integer,save:: posinside
      real(8),parameter:: eps= 1.d-9
!///////////
!     logical,save:: mydebug=.false.
!/////////////

      end module fpolygon

      subroutine eprfpolygon(comp)
      use fpolygon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
!         interface to read configuration data for "fpolygon"
!
       type(Component)::  comp   ! output. to recieve the config data.
!           read fpolygon data as 'new-*'
!           fpolygon  has height and npoly vertex  attributes and the direction cosines
!           of the 'x' and 'y'==> (1-6)
!
      integer i,  j
      character(len=MAX_STRUCCHR):: basename 
!             next is mandatory
      call eprpst(comp, 0, 0, 1, 6)
      call epfpolygonCnst(comp)  ! this call, xmin etc  in Volat
                              ! are undefined.
      j = comp%vol + 3
      xmin=  Volat( j )
      xmax = Volat( j )
      ymin = Volat( j+1 )
      ymax = Volat( j+1)

      j = comp%vol + 3
      do i = 2, npoly*2, 2
         if(xmin > Volat( j+i) ) xmin = Volat(j+i)
         if(xmax < Volat( j+i) ) xmax = Volat(j+i)
         if(ymin > Volat( j+i+1 ) ) ymin = Volat(j+i+1)
         if(ymax < Volat( j+i+1 ) ) ymax = Volat(j+i+1)
      enddo
      j = j + npoly*2 
      Volat( j ) = xmin
      Volat( j+1 ) = ymin
      Volat( j+2 ) = xmax
      Volat( j+3 ) = ymax
      call epGetBaseStrucName(comp%struc, basename)
      if( basename /= 'fpolygon' ) then
         write(0,*) 'structure=',comp%struc, ' not usable'
         stop
      endif
      end   subroutine eprfpolygon

      subroutine epfpolygonCnst(comp)
      use fpolygon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "fpolygon"
!
       type(Component)::  comp   ! output. to recieve the config data.
      integer::j
!           check some values
!      if( Compnum /= comp.cn ) then  ! this is dangerous.
!                        since diff. comp. may get same cn.
!                        at the time of config reading
      if( Compnum /= comp%vol ) then   ! while comp.vol is uniq
         npoly = Volat( comp%vol + 1)
         height = Volat( comp%vol + 2)
!         Compnum = comp.cn
         Compnum = comp%vol
         if( npoly<=2 .or. height<=0 ) then
            write(0,*) ' npoly=',npoly, ' for fpolygon '
            write(0,*) ' height=', height
            write(0,*) ' some of above invalid'
            stop
         endif
         j = comp%vol+ 3  + npoly*2 
         xmin = Volat( j ) 
         ymin = Volat( j+1 )
         xmax = Volat( j+2 ) 
         ymax = Volat( j+3 ) 
      endif
      end subroutine epfpolygonCnst
      subroutine epbfpolygon(comp, posl, dirl, el, icon)
      use fpolygon
      implicit none
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
 

       type(Component):: comp    ! input. you can extract volume parameters
                          !            by Volat( comp.vol + 1), etc
       type(epPos)::  posl       ! input.  position.
       type(epDirec)::  dirl      ! input. direction cosinse

      real(8):: el                !  output length cm from pos to the boundary
      integer:: icon              ! output 0: el obtained. pos    is inside
                     !        1:  //                       outside
                     !       -1: the line dose not cross the volume

!      record /epPos/ cposl
!      record /epDirec/ cdirl
      real(8):: cposl(3), cdirl(3)

      call epv2c_fpolygon(comp, posl, cposl)
      call epv2cd_fpolygon(comp, dirl, cdirl)
!//////////
!      if( abs(cposl(1) - 22.6299991148469d0) < 1.d-8 .and.
!     *    abs(cdirl(1) -0.716320230980483d0 ) < 1.d-8 ) then
!         mydebug= .true.
!      endif
!////////////////
      call epbfpolygon0(comp, cposl, cdirl, el, icon)
      end subroutine epbfpolygon

!      **********************************
      subroutine epsfpolygon(comp, pos, icon)
      use fpolygon
      implicit none
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
       type(epPos)::  cpos


      call epv2c_fpolygon(comp, pos, cpos)
      call epsfpolygon0(comp, cpos, icon)
      end subroutine epsfpolygon
!      **********************************
      subroutine epsfpolygon0(comp, pos, icon)
      use fpolygon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
!           judge if a given 'pos' is inside 'comp'
!         
       type(Component)::  comp !input component
!      record /epPos/ pos  ! input. position; canonical
      real(8),intent(in):: pos(3)
      integer,intent(out):: icon  ! output. 0--> pos is inside
                    !         1-->        outside
      integer j, inside
      call epfpolygonCnst(comp)
      if( pos(3) < 0.d0 .or. pos(3) > height ) then
         icon = 1
         return
      elseif(pos(1) < xmin .or. pos(1) > xmax) then
         icon = 1
         return
      elseif(pos(2) < ymin .or. pos(2) > ymax) then
         icon = 1
         return
      else
         j = comp%vol + 3
        call kinout2(Volat(j),2, Volat(j+1),2,
     *    npoly, pos(1), pos(2), inside)
!         call kinout(Volat(j),2, Volat(j+1),2,
!     *    npoly, pos(1), pos(2), eps, "any", inside) 
!                   inside=1  pos is outside
!                          -1 inside 0 on the edge
!                          -2 invalid polygon

         if(inside == 1 ) then
            icon = 1
         elseif(inside == 0 .or. inside == -1) then 
            icon = 0
         else
            write(0,*) 'in epsfpolygon0 '
            write(0,*) ' strange; npoly=', npoly
            stop
         endif
      endif
      end subroutine epsfpolygon0

!     **************************************
      subroutine epenvlpfpolygon(comp, org, abc)
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

      real(8):: temp, temp2
      call epenvlpfpolygon0(comp, org, abc)
      call epc2v_fpolygon(comp, abc, abc)
      end subroutine epenvlpfpolygon

!    *************************************
      subroutine epatlocfpolygon(comp, loc)
      use fpolygon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(*)
 
      integer i
      call epfpolygonCnst(comp)
      do i = 1, npoly*2 + 2
         loc(i) = i
      enddo
      end subroutine epatlocfpolygon



!*****************************
      subroutine epbfpolygon0(comp, posl, dirl, el, icon)
      use fpolygon
      implicit none
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
 

       type(Component):: comp    ! input. you can extract volume parameters
                          !            by Volat( comp.vol + 1), etc
      real(8),intent(in):: posl(3)     ! input.  position.
      real(8),intent(in):: dirl(3)     ! input. direction cosinse

      real(8):: el                !  output length cm from pos to the boundary
      integer:: icon              ! output 0: el obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume
 
      integer jcon
      real(8):: x, y, z


      integer::i, inside, j
      elmin = bigdist
      posinside  = -10
      if( posl(3) > height ) then
         !  point is above top flat plane;  get point at height
         !  height = posl(3) + el *dirl(3)
         if( dirl(3) >= 0. ) then
            icon  = -1   ! case 1 
            return   !***********
         endif
         el = ( height-posl(3) ) /dirl(3)
!          candidate xing pos         
         x = el * dirl(1) + posl(1)
         y = el * dirl(2) + posl(2)
         ! see if it is inside upper flat polygon
         call kinout2(Volat(comp%vol+3),2,  Volat(comp%vol+4),2,
     *   npoly, x, y, inside)
!         j = comp.vol+3
!         call kinout(Volat(j),2, Volat(j+1),2,
!     *    npoly, x, y, eps, "any", inside) 
         if(inside == -1 .or. inside == 0 ) then
            icon = 1  ! case 3
            return  !   x point obtained .  posl is outside
         endif
         !     still possible to cross at the side wall
         posinside =1
         call ep_fpolygon_side(comp, posl, dirl, el, icon)
         return
      elseif( posl(3) < 0. ) then
         ! posl is below flat polygon
         if( dirl(3) <= 0. ) then
            icon  = -1  ! case  1
            return   !***********
         endif
         el =  -posl(3)  /dirl(3)
!          candidate xing pos         
         x = el * dirl(1) + posl(1)
         y = el * dirl(2) + posl(2)
         ! see if it is inside lower flat polygon
         call kinout2(Volat(comp%vol+3),2,  Volat(comp%vol+4),2,
     *   npoly, x, y, inside)
!         j = comp.vol+3
!         call kinout(Volat(j),2, Volat(j+1),2,
!     *    npoly, x, y, eps, "any", inside) 

         if(inside == -1 .or. inside == 0 ) then
            icon = 1  ! case 3
            return  !   x point obtained .  posl is outside
         endif
                  ! still  possible to cross at the side wall
         posinside = 1
         call ep_fpolygon_side(comp, posl, dirl, el, icon)
         return
      elseif( posl(1) <= xmin  .and. dirl(1) <= 0.d0 ) then
         icon = -1  ! no crossing
         return
      elseif( posl(1) >= xmax  .and. dirl(1) >= 0.d0 ) then
         icon = -1  ! no crossing
         return
      elseif( posl(2) <= ymin  .and. dirl(2) <= 0.d0 ) then
         icon = -1  ! no crossing
         return
      elseif( posl(2) >= ymax  .and. dirl(2) >= 0.d0 ) then
         icon = -1  ! no crossin
         return
      elseif( posl(3) == height) then
         if( dirl(3) >= 0. ) then
            icon = -1           ! case 4
            return              ! ************
         endif
      elseif( posl(3) == 0. ) then
         if( dirl(3) <= 0. ) then
            icon = -1           ! case 4
            return              ! ************
         endif
      endif
        !   posl(3) >=0. .and. posl(3) <= height
      if( posinside == -10 ) then
         call epsfpolygon0(comp, posl, posinside)
      endif
!          posinside = 0.  posl is inside (inc on surface)
!                 1.  outside
           ! posl is inside of fpolygon
      if( posinside == 0 ) then
         if( dirl(3) > 0. ) then
            el = ( height-posl(3) ) /dirl(3)         
            x = el*dirl(1) + posl(1)
            y = el*dirl(2) + posl(2)
            call kinout2(Volat(comp%vol+3),2,  Volat(comp%vol+4),2,
     *           npoly, x, y, inside)
!            j = comp.vol+3
!            call kinout(Volat(j),2, Volat(j+1),2,
!     *           npoly, x, y, eps, "any", inside) 
            if(inside <= 0 ) then
               if(elmin > el ) then
                  elmin = el
               endif
            endif
         elseif( dirl(3) < 0. ) then
            el = -posl(3)/dirl(3)
            x = el*dirl(1) + posl(1)
            y = el*dirl(2) + posl(2)
            call kinout2(Volat(comp%vol+3),2, Volat(comp%vol+4),2,
     *           npoly, x, y, inside)

!            j = comp.vol+3
!            call kinout(Volat(j),2, Volat(j+1),2,
!     *           npoly, x, y, eps, "any", inside) 
            if(inside <= 0 ) then
               if(elmin > el ) then
                  elmin = el
               endif
            endif
         endif
      endif
              ! examine wall
      call ep_fpolygon_side(comp, posl, dirl, el, icon)
      end subroutine epbfpolygon0

      subroutine ep_fpolygon_side(comp, posl, dirl, el, icon)
      use fpolygon
      implicit none
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"

!          seee side wall
 
       type(Component):: comp    ! input. you can extract volume parameters
                          !            by Volat( comp.vol + 1), etc
      real(8),intent(in):: posl(3)     ! input.  position.
      real(8),intent(in):: dirl(3)     ! input. direction cosinse

      real(8):: el                !  output length cm from pos to the boundary
      integer:: icon              ! output 0: el obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume
 
      real(8)::p1(3), p2(3)
      real(8)::p3(3), p4(3)
      integer i, j, jcon, kcon
      real(8),parameter::eps2=2.d-5

      j = comp%vol + 2
      p2(1) = Volat(j+1)
      p2(2) = Volat(j+2)
      p2(3) = 0.d0
      do i = 1, npoly*2, 2
         p1(:) = p2(:) 
!!!  for epxpLand4vp
         p4(1:2) = p1(1:2)
         p4(3) = height
!!!
         if(i == npoly*2-1) then
            p2(1) = Volat(j+1)
            p2(2) = Volat(j+2)
            p2(3) = 0.
!!!  for epxpLand4vp
            p3(1:2)= p2(1:2)
            p3(3) = height
!!!
         else
            p2(1) = Volat(j+i+2)
            p2(2) = Volat(j+i+3)
            p2(3) = 0.

!!!  for epxpLand4vp
            p3(1:2) = p2(1:2)
            p3(3) = height
!!!
         endif

         call epxpLand4vp(p1, p2, p3, p4, posl,dirl,
     *         el, kcon, jcon)

!!!          near the boundary, next sometimes  fials
!!!       call kxplvsq(posl, dirl, p1, p2, height, el, jcon)
!//////////////
!         if(mydebug) then
!            write(0,*) 'after kxplvsq : H=', height, i
!            write(0,*) 'posl=',posl(:)
!            write(0,*) 'dirl=',dirl(:)
!            write(0,*) 'p1=', p1(:)
!            write(0,*) 'p2=', p2(:)
!            write(0,*) ' jcon, el=',jcon, el
!         endif
!///////////////
!!         if(jcon == 0 ) then
         if( kcon <= 4 .and. el > 0.d0) then
            if(elmin > el ) then
               elmin = el
            endif
         endif
      enddo
      el = elmin
      if( el /= bigdist ) then
         if( posinside == 0   ) then
            icon = 0
         else
            icon = 1
         endif
      else
         if(posinside == 0 ) then
            ! strange. xpoint should exist

            ! always very close to the surface and goes almost parallel
            ! to the surface. So we may regards safely that
            ! the particle is outsied and no X-ing point
            ! When the incident is vertical, this apts to happen.
            !
!!!            posinside = 1

!           if(mydebug) then
!             write(0,*) 'posl =',posl(:)
!             write(0,*) 'posl+dr=',posl(:)+1*dirl(:)
!             write(0,*) 'dirl =',dirl(:)
!             write(0,*) 'comp# =',comp.cn
!             write(0,*) ' name=',comp.struc
!             write(0,*) ' posl is inside but no Xpoint'
!             j = comp.vol + 1
!             do i = 1, npoly*2, 2
!                write(0,*) ' vertex ' 
!                write(0,*)  i, Volat(j+i+1), Volat(j+i+2)
!             enddo
!             stop
!           endif
         endif
         icon = -1
      endif
      end subroutine ep_fpolygon_side
!     **************************************
      subroutine epenvlpfpolygon0(comp, org, abc)
      use fpolygon
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


      call epfpolygonCnst(comp)  ! for safety 


      org = epPos(xmin, ymin, 0.d0)
      abc = ep3Vec( xmax - xmin,  ymax - ymin, height)

      NVTX = 0
      end subroutine epenvlpfpolygon0

      subroutine epv2c_fpolygon(comp, posv, posc)
      use fpolygon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp  ! input.   component.
       type(epPos)::  posv     ! input
       type(epPos)::  posc     ! output can be posv

       type(epPos)::  temp

      call epfpolygonCnst(comp)
      if( comp%struc  == 'fpolygon' .or. 
     *    comp%struc  == 'fpolygon_w' ) then
         posc = posv
      elseif( comp%struc(1:11)  == 'fpolygon_xy' ) then
         posc = posv
      elseif( comp%struc(1:11) == "fpolygon_yz" ) then
         temp =epPos(posv%y, posv%z, posv%x)
         posc = temp
      elseif( comp%struc(1:11) == "fpolygon_zx" ) then
         temp = epPos( posv%z, posv%x, posv%y)
         posc = temp 
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbfpolygon is invalid'
         stop
      endif



      end   subroutine epv2c_fpolygon

      subroutine epc2v_fpolygon(comp, posc, posv)
      use fpolygon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp  ! input.   component.
       type(epPos)::  posc     ! input
       type(epPos)::  posv     ! output (can be posc)

       type(epPos)::  temp    !

      call epfpolygonCnst(comp)
      if( comp%struc  == 'fpolygon' .or.
     *    comp%struc  == 'fpolygon_w'  ) then
         posv = posc
      elseif( comp%struc(1:11)  == 'fpolygon_xy' ) then
         posv = posc
      elseif( comp%struc(1:11) == "fpolygon_yz" ) then
!                X:y   Y:z    Z:x
         temp = epPos( posc%z, posc%x, posc%y)
!!         posc = temp !!!!!!!!!!!! aho
         posv = temp
      elseif( comp%struc(1:11) == "fpolygon_zx" ) then
!         cpos = epPos( pos.z, pos.x, pos.y)
!                  X:z    Y:x          Z:y
         temp = epPos(posc%y, posc%z, posc%x)
!!         posc = temp  !!!!!!!!!! aho
         posv = temp  
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbfpolygon is invalid'
         stop
      endif
      end   subroutine epc2v_fpolygon

      subroutine epv2cd_fpolygon(comp, dirv, dirc)
      use fpolygon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
       type(Component)::  comp  ! input.   component.
       type(epDirec)::  dirv     ! input
       type(epDirec)::  dirc     ! output  can be dirv
 
       type(epDirec)::  temp
     
      call epfpolygonCnst(comp)

      if( comp%struc  == 'fpolygon' .or.
     *    comp%struc  == 'fpolygon_w' ) then
         dirc = dirv
      elseif( comp%struc(1:11)  == 'fpolygon_xy' ) then
         dirc = dirv
      elseif( comp%struc(1:11) == "fpolygon_yz" ) then
         temp = epDirec(dirv%y, dirv%z, dirv%x)
         dirc = temp
      elseif( comp%struc(1:11) == "fpolygon_zx" ) then
         temp = epDirec(dirv%z, dirv%x, dirv%y)
         dirc = temp
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbfpolygon is invalid'
         stop
      endif

      end subroutine epv2cd_fpolygon

      subroutine epc2vd_fpolygon(comp, dirc, dirv)
      use fpolygon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
       type(Component)::  comp  ! input.   component.
       type(epDirec)::  dirc     ! input
       type(epDirec)::  dirv     ! output can be dirc

       type(epDirec)::  temp

      call epfpolygonCnst(comp)
!           (X,Y,Z)<---(x,y,z)
      if( comp%struc  == 'fpolygon' .or.
     *    comp%struc  == 'fpolygon_w') then
         dirv = dirc
      elseif( comp%struc(1:11)  == 'fpolygon_xy' ) then
         dirv = dirc
      elseif( comp%struc(1:11) == "fpolygon_yz" ) then
!         cdir = epDirec(dir.y, dir.z, dir.x)  X:y  Y:z  Z:x
         temp = epDirec(dirc%z, dirc%x, dirc%y) 
         dirv = temp   
      elseif( comp%struc(1:11)  == 'fpolygon_zx' ) then
!         cdir = epDirec(dir.z, dir.x, dir.y)  X:z Y:x Z:y
         temp = epDirec(dirc%y, dirc%z, dirc%x) 
         dirv = temp
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbfpolygon is invalid'
         stop
      endif

      end   subroutine epc2vd_fpolygon
