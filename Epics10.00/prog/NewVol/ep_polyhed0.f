!         polyhedra0: 
!      
      module modpolyhed0
      integer,save::Compnum=-1000
      real(8),save:: xmin, xmax, ymin, ymax, zmin, zmax
      integer,save:: npoly
      integer,save:: convex
      integer,save:: p1st, q1st
      end module modpolyhed0

      subroutine eprpolyhed0(comp)
      use modpolyhed0
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp   ! output. to recieve the config data.
!         polyhed0 has npoly*2 3D vertexes  attributes and optional
!         direction cosines
      integer i,  j, jmax
      character(len=MAX_STRUCCHR):: basename 
!             next is mandatory
      call eprpst(comp, 0, 6, 1, 6)
!            data should be 
!              N convex  x y z x y z ... \
!                 x y z x y z ... \
!     N is the total # of items convex, x, y, z.. and N itself.
!     The # of vertexes for p is n=(N-2)/2. First n
!     items are for p. which must be followed by n 
!     items for q.  6 is for (xmin,ymin,zmin), (xmax,ymax,zmax)
!     convex: 1 --> p,q forms convex polygon.  0--> p,q may form
!     concave or convec  polygon. 
      call eppolyhed0Cnst(comp)  ! at this moment,, xmin etc  in Volat
                  ! are undefined. so compute them here and put
                  ! them in Volat
      if( npoly<=2 ) then
         write(0,*) ' npoly=',npoly, ' for polyhed0 invalid'
         stop
      endif


      jmax = (2*npoly-1)*3
      xmin = minval( Volat(p1st:p1st+jmax:3) )
      xmax = maxval( Volat(p1st:p1st+jmax:3) )

      ymin = minval( Volat(p1st+1:p1st+jmax+1:3) )
      ymax = maxval( Volat(p1st+1:p1st+jmax+1:3) )

      zmin = minval( Volat(p1st+2:p1st+jmax+2:3) )
      zmax = maxval( Volat(p1st+2:p1st+jmax+2:3) )

      j = p1st + jmax + 3
      Volat( j ) = xmin
      Volat( j+1 ) = ymin
      Volat( j+2 ) = zmin
      Volat( j+3 ) = xmax
      Volat( j+4 ) = ymax
      Volat( j+5 ) = zmax
      call epGetBaseStrucName(comp%struc, basename)
      if( basename /= 'polyhed0' ) then
         write(0,*) 'structure=',comp%struc, ' not usable'
         stop
      endif
      end   subroutine eprpolyhed0
      subroutine eppolyhed0Cnst(comp)
      use modpolyhed0
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!    interface to read configuration data for "polyhed0"
!
       type(Component)::  comp   ! output. to recieve the config data.
      integer::j
!           check some values
!      if( Compnum /= comp.cn ) then  ! this is dangerous.
!                        since diff. comp. may get same cn.
!                        at the time of config reading
      if( Compnum /= comp%vol ) then   ! while comp.vol is uniq
         npoly = (Volat( comp%vol + 1)-2)/3/2
         p1st= comp%vol + 3
         q1st = p1st+ npoly*3
!         Compnum = comp.cn
         Compnum = comp%vol
         j = comp%vol + 3  + npoly*2 *3
         xmin = Volat( j ) 
         ymin = Volat( j+1 )
         zmin = volat( j+2 )
         xmax = Volat( j+3 ) 
         ymax = Volat( j+4 ) 
         zmax = Volat( j+5 ) 
      endif
      end subroutine eppolyhed0Cnst

      subroutine epbpolyhed0(comp, posl, dirl, el, icon)
      use modpolyhed0
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
      real(8),intent(in):: posl(3)
      real(8),intent(in):: dirl(3)

      real(8),intent(out):: el !  output length cm from pos to the boundary
      integer,intent(out):: icon  ! output 0: el obtained. 
                          !     pos  is inside
                     !        1:  //                       outside
                     !       -1: the line dose not cross the volume

      real(8):: cposl(3), cdirl(3)

!      call epv2c_polyhed0(comp, posl, cposl)
!      call epv2cd_polyhed0(comp, dirl, cdirl)
!      call epbpolyhed00(comp, cposl, cdirl, el, icon)
      call eppolyhed0Cnst(comp)

      call epbpolyhed00(comp, posl, dirl, el, icon)

      end subroutine epbpolyhed0

!      **********************************
      subroutine epspolyhed0(comp, pos, icon)
      use modpolyhed0
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


!      call epv2c_fpolygon(comp, pos, cpos) ! not needed ?
      call eppolyhed0Cnst(comp)
      call epspolyhed00(comp, pos, icon)
      end subroutine epspolyhed0

!      **********************************
      subroutine epspolyhed00(comp, pos, icon)
      use modpolyhed0
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
      if( pos(1) < xmin .or. pos(1) > xmax ) then
         icon = 1
      elseif(pos(2) < ymin .or. pos(2) > ymax) then
         icon = 1
      elseif(pos(3) < zmin .or. pos(3) > zmax) then
         icon = 1
      else
         call kisInPolyhed0(Volat(p1st), Volat(q1st), npoly,
     *  pos, icon)
!!!!!!!!!!
      write(0,*) ' in poly by kisIn  icon=',icon
!!!!!!!!!!
      endif
!!!!!!!!!!
      write(0,*) ' in poly? icon=',icon
!!!!!!!!!!
      end subroutine epspolyhed00

!     **************************************
      subroutine epenvlppolyhed0(comp, org, abc)
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
      call epenvlppolyhed00(comp, org, abc)
!      call epc2v_polyhed0(comp, abc, abc)
      end subroutine epenvlppolyhed0

!    *************************************
      subroutine epatlocpolyhed0(comp, loc)
      use modpolyhed0
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(*)
 
      integer i
      call eppolyhed0Cnst(comp)
      do i = 1, npoly*2*3 + 2
         loc(i) = i
      enddo
      end subroutine epatlocpolyhed0




!*****************************
      subroutine epbpolyhed00(comp, posl, dirl, el, icon)
      use modpolyhed0
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
 
      integer::icon1,icon2
      
      call kxplPolyhed0(Volat(p1st), Volat(q1st), npoly,
     *      convex, posl,  dirl, el, icon1, icon2)
      if(icon2 == -1 ) then
         icon = icon2
      else
         icon = icon1
      endif 
      
      end subroutine epbpolyhed00

!     **************************************
      subroutine epenvlppolyhed00(comp, org, abc)
      use modpolyhed0
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

      call eppolyhed0Cnst(comp)  ! for safety 

      org = epPos(xmin, ymin, zmin)
      abc = ep3Vec( xmax - xmin,  ymax - ymin, zmax-zmin)
      NVTX = 0
      end subroutine epenvlppolyhed00

      subroutine epv2c_polyhed0(comp, posv, posc)
      use modpolyhed0
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp  ! input.   component.
       type(epPos)::  posv     ! input
       type(epPos)::  posc     ! output can be posv

       type(epPos)::  temp

!      call epfpolygonCnst(comp)
      if( comp%struc  == 'polyhed0' .or. 
     *    comp%struc  == 'polyhed0_w' ) then
         posc = posv
      elseif( comp%struc(1:11)  == 'polyhed0_xy' ) then
         posc = posv
      elseif( comp%struc(1:11) == "polyhed0_yz" ) then
         temp =epPos(posv%y, posv%z, posv%x)
         posc = temp
      elseif( comp%struc(1:11) == "polyhed0_zx" ) then
         temp = epPos( posv%z, posv%x, posv%y)
         posc = temp 
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbpolyhed0 is invalid'
         stop
      endif

      end   subroutine epv2c_polyhed0

      subroutine epc2v_polyhed0(comp, posc, posv)
      use modpolyhed0
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp  ! input.   component.
       type(epPos)::  posc     ! input
       type(epPos)::  posv     ! output (can be posc)

       type(epPos)::  temp    !

!      call eppolyhed0Cnst(comp)
      if( comp%struc  == 'polyhed0' .or.
     *    comp%struc  == 'polyhed0_w'  ) then
         posv = posc
      elseif( comp%struc(1:11)  == 'polyhed0_xy' ) then
         posv = posc
      elseif( comp%struc(1:11) == "polyhed0_yz" ) then
!                X:y   Y:z    Z:x
         temp = epPos( posc%z, posc%x, posc%y)
!!         posc = temp !!!!!!!!!!!! aho
         posv = temp
      elseif( comp%struc(1:11) == "polyhed0_zx" ) then
!         cpos = epPos( pos.z, pos.x, pos.y)
!                  X:z    Y:x          Z:y
         temp = epPos(posc%y, posc%z, posc%x)
!!         posc = temp  !!!!!!!!!! aho
         posv = temp  
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbpolyhed0 is invalid'
         stop
      endif
      end   subroutine epc2v_polyhed0

      subroutine epv2cd_polyhed0(comp, dirv, dirc)
      use modpolyhed0
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
       type(Component)::  comp  ! input.   component.
       type(epDirec)::  dirv     ! input
       type(epDirec)::  dirc     ! output  can be dirv
 
       type(epDirec)::  temp
     
!      call eppolyhed0Cnst(comp)

      if( comp%struc  == 'polyhed0' .or.
     *    comp%struc  == 'polyhed0_w' ) then
         dirc = dirv
      elseif( comp%struc(1:11)  == 'polyhed0_xy' ) then
         dirc = dirv
      elseif( comp%struc(1:11) == "polyhed0_yz" ) then
         temp = epDirec(dirv%y, dirv%z, dirv%x)
         dirc = temp
      elseif( comp%struc(1:11) == "polyhed0_zx" ) then
         temp = epDirec(dirv%z, dirv%x, dirv%y)
         dirc = temp
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbpolyhed0 is invalid'
         stop
      endif

      end subroutine epv2cd_polyhed0

      subroutine epc2vd_polyhed0(comp, dirc, dirv)
      use modpolyhed0
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
       type(Component)::  comp  ! input.   component.
       type(epDirec)::  dirc     ! input
       type(epDirec)::  dirv     ! output can be dirc

       type(epDirec)::  temp

!      call eppolyhedCnst(comp)
!           (X,Y,Z)<---(x,y,z)
      if( comp%struc  == 'polyhed0' .or.
     *    comp%struc  == 'polyhed0_w') then
         dirv = dirc
      elseif( comp%struc(1:11)  == 'polyhed0_xy' ) then
         dirv = dirc
      elseif( comp%struc(1:11) == "polyhed0_yz" ) then
!         cdir = epDirec(dir.y, dir.z, dir.x)  X:y  Y:z  Z:x
         temp = epDirec(dirc%z, dirc%x, dirc%y) 
         dirv = temp   
      elseif( comp%struc(1:11)  == 'polyhed0_zx' ) then
!         cdir = epDirec(dir.z, dir.x, dir.y)  X:z Y:x Z:y
         temp = epDirec(dirc%y, dirc%z, dirc%x) 
         dirv = temp
      else
         write(0,*) ' comp%struc=', comp%struc,
     *      ' to epbpolyhed0 is invalid'
         stop
      endif

      end subroutine epc2vd_polyhed0
