!
!  torus
!                                        
! canonical one.   center circle  of the torus body is at 
!      (x-y) plane wit radius R.  the torus body has radius r.
!      
!
!
!   Data format in config is:
!       ox oy oz  R  r
!
!      where (ox,oy,oz) is the origin in the world coord.
!      
      subroutine eprtorus(comp)
      implicit none
#include "Zglobalc.h"
#include "ZepTrackv.h"
!   #include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepManager.h"

!
!         interface to read configuration data for "torus"
!
       type(Component)::  comp   ! output. to recieve the config data.
      character*80 msg

      integer iR, isr
      parameter( iR = 1,  isr = 2)


      real*8 R, sr
      character(len=MAX_STRUCCHR):: basename
!
      call eprpst(comp, 2, 2, 7, 9) ! 7 9 is Z axis direction
      call epGetBaseStrucName(comp%struc, basename)
      if( basename /= "torus" ) then
         write(msg,*) ' error structure=', comp%struc
         call cerrorMsg(msg, 0)
      endif
!        why confdata has been  used so far ?????
!      if(index(confdata, 'torus_z') .ne. 0 .or.
!     *        trim(Field(2)) .eq. 'torus' ) then
!c               this is canonical form; nothing to do                 
!      elseif(index(confdata, 'torus_x') .ne. 0) then
!      elseif(index(confdata,'torus_y') .ne. 0) then
!      else
!         write(msg,*) ' error structure=', confdata
!         call cerrorMsg(msg, 0)
!      endif

!           check some values
      R = Volat( comp%vol +  iR)
      sr = Volat( comp%vol +  isr)
      if(R  .le. 0. .or. sr .le. 0. .or. sr > R) then
         write(msg, *) comp%cn, '-th component=', comp%struc
         call cerrorMsg(msg, 1)
         write(msg, *) ' R, r=', R, sr, ' invalid'
         call cerrorMsg(msg, 0)
      endif
!             compute const for later use. nothing yet
      end
!      ***********
      subroutine epbtorus(comp, posli, dirli,  el, jcon)
!      ***********
!        find length to the boundary of comp component
!        from posli with dire cos dirli in local coord.

      implicit none
#include "Zglobalc.h"
#include "ZepPos.h"
#include "Zep3Vec.h"
#include "ZepDirec.h"
#include "Zcnfig.h"
!
       type(Component):: comp
      integer jcon
      real*8 el
!

       type(epPos)::   posli
       type(epDirec)::  dirli
      integer base
       type(epPos)::   posl
       type(epDirec)::  dirl

      real(8):: k(4)
      integer::i, nx, icon, inside

      call epv2c_torus(comp, posli, posl)
      call epv2cd_torus(comp,dirli, dirl)


      call kxpLineAndTorus( posl, dirl, Volat(comp%vol+1), 
     *                     Volat(comp%vol+2),  k, nx)

      call kisInsideTorus( posl, 
     *     Volat(comp%vol+1),  Volat(comp%vol+2), inside)
      jcon = -1
      do i = 1, nx
         if( k(i) >= 0.d0 ) then
            if( inside <= 1 ) then   ! inside or surface
               jcon = 0
            else    !  inside=2 outside
               jcon = 1
            endif   
            el = k(i)
            return   ! *****************
         endif
      enddo
      if(jcon == -1 .and. inside <= 1 ) then
         ! contradict.  pos is inside but 
         ! no crossing or crossing pt is backside.
         !   numerical error. (of order of 10^-6)
         ! inside is more reliable so we adjust 
         if( nx > 0 ) then
            el = abs(k(nx))
         else
            el = 1.d-4    ! 1 micron meter.
         endif
      endif
      jcon = 0
      end
      subroutine epstorus(comp, posi, icon)
        implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
!           judge if a given 'pos' is inside 'comp'
!         
       type(Component)::  comp !input component
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside
       type(epPos)::  posi
       type(epPos)::  pos

!          pos in local-c. 
!  &&&&         convert to canonial one
        call epv2c_torus(comp, posi, pos)
        call kisInsideTorus( pos,  Volat( comp%vol+1), 
     *  Volat( comp%vol+2), icon)
        if( icon <= 1) then
           icon = 0
        else
           icon = 1
        endif
        end

      subroutine epv2c_torus(comp, posv, posc)
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  posv ! input
       type(epPos)::  posc ! output. 
       type(Component)::  comp
      if( comp%struc(1:7) == "torus_x") then
!         posc = epPos( posv.z, posv.y, posv.x)
         posc = epPos( posv%y, posv%z, posv%x)
      elseif(comp%struc(1:7) == "torus_y") then
!         posc = epPos( posv.x, posv.z, posv.y)
         posc = epPos( posv%z, posv%x, posv%y)
      else
         posc = posv
      endif
      end
      subroutine epc2v_torus(comp, posc, posv)
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  posc ! input
       type(epPos)::  posv ! output. 
       type(Component)::  comp
      if( comp%struc(1:7) == "torus_x") then
!         posv = epPos( posc.z, posc.y, posc.x)
         posv = epPos( posc%z, posc%x, posc%y)
      elseif(comp%struc(1:7) == "torus_y") then
!         posv = epPos( posc.x, posc.z, posc.y)
         posv = epPos( posc%y, posc%z, posc%x)
      else
         posc = posv
      endif
      end
      subroutine epv2cd_torus(comp, dirv, dirc)
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  dirv ! input
       type(epPos)::  dirc ! output. 
       type(Component)::  comp
      if( comp%struc(1:7) == "torus_x") then
!         dirc = epPos( dirv.z, dirv.y, dirv.x)
         dirc = epPos( dirv%y, dirv%z, dirv%x)
      elseif(comp%struc(1:7) == "torus_y") then
         dirc = epPos( dirv%z, dirv%x, dirv%y)
      else
         dirc = dirv
      endif
      end
      subroutine epc2vd_torus(comp, dirc, dirv)
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  dirc ! input
       type(epPos)::  dirv ! output. 
       type(Component)::  comp
      if( comp%struc(1:7) == "torus_x") then
!         dirv = epPos( dirc.z, dirc.y, dirc.x)
         dirv = epPos( dirc%z, dirc%x, dirc%y)
      elseif(comp%struc(1:7) == "torus_y") then
         dirv = epPos( dirc%y, dirc%z, dirc%x)
      else
         dirv = dirc
      endif
      end
      subroutine epenvlpTorus(comp, org, abc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"


       type(Component)::  comp
       type(epPos)::  org
       type(ep3Vec)::  abc
      real(8):: temp

!                     R                  r
      org%x = -(Volat( comp%vol+ 1 ) + Volat(comp%vol+2))
      org%y = org%x
      org%z = -Volat(comp%vol+2)
      abc%x = -org%x*2
      abc%y = abc%x
      abc%z = Volat( comp%vol+2 ) *2

      NVTX = 0
      if(comp%struc(1:7) == "torus_y") then
         call epc2v_torus(comp, org, org)
!         temp = abc.y
!         abc = ep3Vec(abc.x, abc.z, temp)
!         abc = ep3Vec(temp, abc.z, abc.x)
         call epc2v_torus(comp, abc, abc)
      elseif(comp%struc(1:7) == "torus_x") then
         call epc2v_torus(comp, org, org)
!         temp=abc.x
!         abc = ep3Vec(abc.z, abc.y, temp)
!         abc = ep3Vec(abc.z, temp, abc.y)
         call epc2v_torus(comp, abc, abc)
      endif
      end

      subroutine epatlocTorus(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"


       type(Component)::  comp   ! input.                                
      integer loc(*)
      integer i

      do i = 1, 2
         loc(i) = i
      enddo
      end


