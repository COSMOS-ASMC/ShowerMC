!    This is to manage a new non-default  volume/shape which is not
!    included in Epics.

      subroutine epbNew(comp,  posl, dirl, length, icon)
      implicit none
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
       type(Component):: comp  ! input. given component
       type(epPos)::  posl  ! input. position given in local coord.
       type(epDirec)::  dirl ! input. direction cos. given in
                            !   local  coordinate
      real*8  length  ! output. minimum length to the boundary of 'comp'
!                          from posl with dirl
      integer icon    ! output.
!                       0-->length obtained. pos  is  inside 
!                       1-->length obtained. pos  is  outside
!                      -1-->not cross 
! 
      integer uscl
!           compute the minimum length to the boundary
      call epseeUnderScore(comp%struc, uscl)
      if(comp%struc(1:uscl) .eq. '%1') then
         call epb%1(comp, posl, dirl,  length, icon)
!---  elseif(comp%struc(1:uscl) .eq. '%2') then
!---     call epb%2(comp, posl, dirl,  length, icon)
      else
         call cerrorMsg(comp%struc(1:8), 1)
         call cerrorMsg(' is not supported: epbNew. Consider', 1)
         call cerrorMsg(
     *   "  usenewvol 'your-config-file'", 0)
      endif

      end
!     ***********************************
      subroutine epsNew(comp, pos, icon)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!             see if pos is inside the comp component

       type(Component)::   comp
       type(epPos):: pos  ! input. given point
      integer icon  !   output.  0-->inside. 1-->outside

      integer uscl

      call epseeUnderScore(comp%struc, uscl)
      if(comp%struc(1:uscl) .eq. '%1')  then
         call eps%1(comp, pos, icon)
!---  elseif(comp%struc(1:uscl) .eq. '%2')  then
!---     call eps%2(comp, pos, icon)
      else
         call cerrorMsg(trim(comp%struc), 1)
         call cerrorMsg(' is not supported: epsNew', 0)
      endif
         
      end
!     ***********************************
      subroutine eprNew(comp, shape)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
      
       type(Component)::  comp  ! input/output. where config data is put

      character*(*)  shape     ! input.  config data has this shape
!            read configuration data 
      integer uscl
      call epseeUnderScore(shape, uscl)
      if(shape(1:uscl)  .eq. '%1') then
         call epr%1(comp)
!---  elseif(shape(1:uscl)  .eq. '%2') then
!---     call epr%2(comp)
      else
         call cerrorMsg(trim(shape), 1)
         call cerrorMsg(' is not supported: eprNew ', 0)
      endif

      end
!     ***********************************
      subroutine epenvlpNew(comp,  orig, abc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

!         returns  box parameters which envelopes the
!         given component comp. (local cordinate).
!         
       type(Component)::  comp  ! input. given component
       type(epPos)::  orig      ! output. origin
       type(ep3Vec)::  abc      ! output. a,b,c of the box.
      integer uscl

      call epseeUnderScore(comp%struc, uscl)

      if(comp%struc(1:uscl) .eq. '%1') then
         call epenvlp%1(comp, orig, abc)
!---  elseif(comp%struc(1:uscl) .eq. '%2') then
!---     call epenvlp%2(comp, orig, abc)
      else
         call cerrorMsg(comp%struc(1:8), 1)
         call cerrorMsg(' is not supported: epenvlpNew', 0)
      endif
      end
!     ***********************************
      subroutine epatlocNew(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         
       type(Component)::  comp  ! input. given component
      integer loc(*)     !  output. loc(i) is the location of the 
                         !  attribute given in the i-th position in
                         !  the config file. It's value is in 
                         !  Volat( comp.vol + loc(i)).
                         !  loc size must be >= comp.Nattributes
      integer uscl

      call epseeUnderScore(comp%struc, uscl)

      if(comp%struc(1:uscl) .eq. '%1') then
         call epatloc%1(comp, loc)
!---  elseif(comp%struc(1:uscl) .eq. '%2') then
!---     call epatloc%2(comp, loc)
      else
         call cerrorMsg(comp%struc, 1)
         call cerrorMsg(' is not supported: epatlocNew', 0)
      endif
      end

