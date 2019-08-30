      module prism
      integer,save::Compnum=-1000
      real(8),save:: a, b, c, h    ! 
      end module prism

      subroutine  eprprs(comp)
      use prism
      use modNegativeC
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
       type(Component)::  comp   ! output. to recieve the config data.

      character(len=MAX_STRUCCHR):: basename
      
      call eprpst(comp, 4, 4, 1, 6)
      call epGetBaseStrucName(comp%struc, basename)
      if( basename /= 'prism' ) then
         write(0,*) ' at prism program'
         write(0,*) ' error structure=', comp%struc
         stop
      endif
      call epprismCnst(comp)
      if( b < 0. .or. c < 0. .or. h < 0. .or. a < 0.) then
         negativePrism = negativePrism + 1
      endif
      if( a < 0.) then
         comp%orgy = comp%orgy + a
         a = -a
         Volat( comp%vol + 1 ) = a

      endif
      if( b < 0.) then
         comp%orgy = comp%orgy + b
         b = -b
         Volat( comp%vol +2 ) = b
      endif
      end      subroutine  eprprs

      subroutine epbprs(comp, pos, dir, length, icon)
      implicit none
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"
       type(Component):: comp    ! input.  prism
       type(epPos)::  pos        ! input.  position.
       type(epDirec)::  dir      ! input. direction cosinse                                          
      real(8),intent(out):: length !  length cm from pos to the boundary
      integer,intent(out):: icon      !  0: length obtained. pos is inside
                             !        1:  //                        outside
                             !       -1: the line dose not cross the volume                                 
       type(epPos)::  cpos
       type(epDirec)::  cdir

      call epv2c_prism(comp, pos, cpos)
      call epv2cd_prism(comp,  dir, cdir)
      call epbprs0(comp, cpos, cdir, length, icon)
      end subroutine epbprs

      subroutine epv2c_prism(comp, pos, cpos)
      use prism
      implicit none
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
!            local coord of variant to local coord 
!            of canonical prism
       type(Component):: comp    ! input.  honeycomb                                                 
       type(epPos)::  pos        ! input.  position.
       type(epPos)::  cpos       ! output  pos-->canonical one

      real(8):: temp
      call epprismCnst(comp)
      if( comp%struc == "prism" .or.
     *    comp%struc == "prism_w" ) then
         cpos = pos
      elseif(  comp%struc(1:8) == "prism_xy" ) then
         cpos = pos
      elseif(  comp%struc(1:7) == "prism_y" .or. 
     *       comp%struc(1:8) == "prism_yx"  ) then
         cpos =epPos(pos%y, b-pos%x, pos%z)
      elseif(  comp%struc(1:8) == "prism_xz" ) then
         temp =  b - pos%z
         cpos = epPos( pos%x, temp, pos%y)
      elseif(   comp%struc(1:8) == "prism_zx" ) then
         cpos = epPos( pos%z, pos%x, pos%y)
      elseif(   comp%struc(1:8) == "prism_yz" ) then
         cpos =epPos(pos%y, pos%z, pos%x)
      elseif(   comp%struc(1:7) == "prism_z" .or.
     *          comp%struc(1:8) == "prism_zy" ) then
         temp =  b - pos%y
         cpos =epPos(pos%z, temp, pos%x)
      else
         write(0,*) ' struc=',comp%struc, ' invalid'
         stop
      endif
      end subroutine epv2c_prism

      subroutine epv2cd_prism(comp, dir, cdir)
      use prism
      implicit none
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
!            local coord of variant to local coord 
!            of canonical prism
       type(Component):: comp    ! input.  honeycomb                                                 
       type(epDirec)::  dir      ! input. direction cosinse                                          
       type(epDirec)::  cdir   ! output    dir--> //

      call epprismCnst(comp)
      if( comp%struc == "prism" .or.
     *    comp%struc == "prism_w" ) then
         cdir = dir
      elseif(  comp%struc(1:8) == "prism_xy" ) then
         cdir = dir
      elseif(  comp%struc(1:7) == "prism_y" .or. 
     *       comp%struc(1:8) == "prism_yx"  ) then
         cdir =epDirec(dir%y, -dir%x, dir%z)
      elseif(  comp%struc(1:8) == "prism_xz" ) then
         cdir = epDirec(dir%x, -dir%z, dir%y)
      elseif(   comp%struc(1:8) == "prism_zx" ) then
         cdir = epDirec(dir%z, dir%x, dir%y)
      elseif(   comp%struc(1:8) == "prism_yz" ) then
         cdir = epDirec(dir%y, dir%z, dir%x)
      elseif(   comp%struc(1:7) == "prism_z" .or.
     *          comp%struc(1:7) == "prism_zy" ) then
         cdir = epDirec(dir%z, -dir%y, dir%x)
      else
         write(0,*) ' struc=',comp%struc, ' invalid'
         stop
      endif
      end subroutine epv2cd_prism



      subroutine epc2v_prism(comp, cpos, pos)
      use prism
      implicit none
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
!             canonical coord to variant local
       type(Component):: comp    ! input.  honeycomb                                                 
       type(epPos)::  cpos        ! input.  canonical position.
       type(epPos)::  pos        ! output  pos-->varant local one

      real(8):: temp
      call epprismCnst(comp)
      if( comp%struc == "prism" .or.
     *    comp%struc == "prism_w"  ) then
         pos = cpos
      elseif(  comp%struc(1:8) == "prism_xy" ) then
         pos = cpos
      elseif(  comp%struc(1:7) == "prism_y" .or. 
     *       comp%struc(1:8) == "prism_yx"  ) then
!         cpos =epPos(pos.y, b-pos.x, pos.z)  V-->C
!                     X : y    Y : b-x   Z :z
         pos =epPos(b-cpos%y, cpos%x, cpos%z)  ! C-->V
      elseif(  comp%struc(1:8) == "prism_xz" ) then
!         temp =  b - pos.z
!         cpos = epPos( pos.x, temp, pos.y)
!              X: x     Y: b-z   Z: y
         pos = epPos(cpos%x, cpos%z,  b-cpos%y) 
      elseif(   comp%struc(1:8) == "prism_zx" ) then
!         cpos = epPos( pos.z, pos.x, pos.y)
!                  X:z    Y:x          Z:y
         pos = epPos(cpos%y,  cpos%z, cpos%x)
      elseif(   comp%struc(1:8) == "prism_yz" ) then
!         cpos =epPos(pos.y, pos.z, pos.x)
!                X:y   Y:z    Z:x
         pos= epPos( cpos%z, cpos%x, cpos%y)
      elseif(   comp%struc(1:7) == "prism_z" .or.
     *          comp%struc(1:8) == "prism_zy" ) then
!         temp =  b - pos.y
!         cpos =epPos(pos.z, temp, pos.x)
!              X: z      Y: b-y   Z:x
         pos = epPos(cpos%z, b-cpos%y, cpos%x)
      else
         write(0,*) ' struc=',comp%struc, ' invalid'
         stop
      endif
      end subroutine epc2v_prism

      subroutine epc2vd_prism(comp, cdir, dir)
      use prism
      implicit none
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
!            local coord of variant to local coord 
!            of canonical prism
       type(Component):: comp    ! input.  honeycomb                                                 
       type(epDirec)::  cdir      ! input. direction cosinse in canonical coord.
                                         
       type(epDirec)::  dir   ! output.  local coord.
      real(8):: temp

      call epprismCnst(comp)
      if( comp%struc == "prism" .or.
     *    comp%struc == "prism_w" ) then
         dir = cdir
      elseif(  comp%struc(1:8) == "prism_xy" ) then
         dir = cdir
      elseif(  comp%struc(1:8) == "prism_y" .or. 
     *       comp%struc(1:8) == "prism_yx"  ) then
!         cdir =epDirec(dir.y, -dir.x, dir.z)   X:y   Y:-x  Z:z
         dir =epDirec(-cdir%y, cdir%x, cdir%z)
      elseif(  comp%struc(1:8) == "prism_xz" ) then
!         temp =  b - pos.z
!         cdir = epDirec(dir.x, -dir.z, dir.y)  X:x  Y:-z   Z:y
         dir = epDirec(cdir%x, cdir%z, -cdir%y) 
      elseif(   comp%struc(1:8) == "prism_zx" ) then
!         cdir = epDirec(dir.z, dir.x, dir.y)  X:z Y:x Z:y
         dir = epDirec(cdir%y, cdir%z, cdir%x) 
      elseif(   comp%struc(1:8) == "prism_yz" ) then
!         cdir = epDirec(dir.y, dir.z, dir.x)  X:y  Y:z  Z:x
         dir = epDirec(cdir%z, cdir%x, cdir%y) 
      elseif(   comp%struc(1:7) == "prism_z" .or.
     *          comp%struc(1:8) == "prism_zy" ) then
!         temp =  b - pos.y
!         cdir = epDirec(dir.z, -dir.y, dir.x) X:z  Y;-y  Z:x
         dir = epDirec(cdir%z, -cdir%y, cdir%x)
      else
         write(0,*) ' struc=',comp%struc, ' invalid'
         stop
      endif
      end subroutine epc2vd_prism

      subroutine epprismCnst(comp)
      use prism
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
       type(Component)::  comp   ! output. to recieve the config data.

      if( Compnum /= comp%cn ) then
         a = Volat( comp%vol + 1)
         b = Volat( comp%vol + 2)
         c = Volat( comp%vol + 3)
         h = Volat( comp%vol + 4)
         Compnum = comp%cn
      endif
      end    subroutine epprismCnst


      subroutine epsprism(ncx, pos, icon)
!          search specified prism
!         Is  pos in ncx-th comp.  ? if yes, icon=0 else
!               icon = 1
!         pos is assumed to be Local coord.
       implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epPos)::  pos
           integer ncx,  icon
       type(epPos)::  cpos
!              to canonical coord.              
      call epv2c_prism(Det%cmp(ncx),  pos, cpos)
      call epsprism0(ncx, cpos, icon)
      end    subroutine epsprism

      subroutine epenvlpPrism(comp, org, abc)
      use prism
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
       type(Component)::  comp
       type(epPos)::  org
       type(ep3Vec)::  abc
      real(8)::temp, temp2
       type(epPos)::  vtx, tvtx
      integer i
   
      call epenvlpPrism0(comp, org, abc)
   
      if( comp%struc == "prism" .or.
     *    comp%struc == "prism_w" ) then
      elseif(  comp%struc(1:8) == "prism_xy" ) then
      elseif(  comp%struc(1:7) == "prism_y" .or. 
     *       comp%struc(1:8) == "prism_yx"  ) then
         temp = org%x
         org = epPos( org%y, temp, org%z)
         temp = abc%x
         abc = ep3Vec( abc%y, temp, abc%z)
      elseif(  comp%struc(1:8) == "prism_xz" ) then
         temp = org%y
         org = epPos( org%x, org%z, temp)
         temp = abc%y
         abc = ep3Vec( abc%x, abc%z, temp)
      elseif(   comp%struc(1:8) == "prism_zx" ) then
!                z->x  x->y
         temp = org%x
         temp2 = org%y
         org = epPos( temp2, org%z, temp)
         temp = abc%x
         temp2 = abc%y
         abc = ep3Vec( temp2, abc%z, temp)
      elseif(   comp%struc(1:8) == "prism_yz" ) then
!         temp =  Volat(comp.vol+2) - org.y
         temp =  org%y
         temp2 = org%x
         org = epPos(org%z, temp2, temp)
         temp = abc%y
         temp2= abc%x
         abc = ep3Vec( abc%z, temp2, temp)
      elseif(   comp%struc(1:7) == "prism_z" .or.
     *          comp%struc(1:8) == "prism_zy" ) then
!         temp =  Volat(comp.vol+prismb) - org.y
         temp =  org%y
         temp2 = org%x
         org = epPos(org%z, temp, temp2)
         temp = abc%y
         temp2 = abc%x
         abc = ep3Vec( abc%z, temp, temp2)
      else
         write(0,*) ' struc=',comp%struc, ' invalid'
         write(0,*) ' detected at epenvlpPrism'
         stop
      endif


      end  subroutine epenvlpPrism

      subroutine epatlocprism(comp,loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
       type(Component)::  comp   ! input.
      integer loc(*)
      integer i
      do i = 1, 4
         loc(i) = i
      enddo
      end subroutine epatlocprism

      subroutine epbprs0(comp, posl, dirl, length, icon)
      use prism
      implicit none
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"
       type(Component):: comp    ! input.  honeycomb
       type(epPos)::  posl        ! input.  position
       type(epDirec)::  dirl      ! input. direction cosinse
      real(8),intent(out):: length    !  length cm from pos to the boundary
      integer,intent(out):: icon ! t 0: length obtained. pos    is inside
                      !        1:  //                        outside                                                      !       -1: the line dose not cross the volume
      call kxplPrism(posl%x, posl%y, posl%z,
     *     dirl%x, dirl%y, dirl%z,
     *     a, b, c, h, length, icon)

      end      subroutine epbprs0

      subroutine epenvlpPrism0(comp, org, abc)
      use prism
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

       type(Component)::  comp
       type(epPos)::  org
       type(ep3Vec)::  abc

      call epprismCnst(comp)

!      org.x = min(Volat(base+prisma), Volat(base + prismc), 0.d0)
!      org.y = min(Volat(base+prismb), 0.d0)
!      org.z = min(Volat(base+prismh), 0.d0)
      org = epPos(  min( 0.d0, c), 0.d0,  min( h, 0.d0))
      abc = ep3Vec( max(a, c ) - org%x, b,
     *              max(h, 0.d0)  - org%z )
      NVTX = 0
      end  subroutine epenvlpPrism0
      subroutine epsprism0(ncx, pos, icon)
      use prism
!          search specified prism
!         Is  pos in ncx-th comp.  ? if yes, icon=0 else
!               icon = 1
!         pos is assumed to be Local coord.
       implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epPos)::  pos
           integer ncx,  icon

        call kioPrism(pos%x, pos%y, pos%z, 
     *        a, b, c, h, icon)
       end    subroutine  epsprism0

