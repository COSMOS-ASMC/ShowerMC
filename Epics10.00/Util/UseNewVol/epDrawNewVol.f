!    This is to manage a new non-default  volume/shape which is not
!    included in Epics.
!   This is the interface to draw new volumes.
      subroutine epDrawNew(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component):: comp  ! input. given component to be drawn
      integer n   !  ouput.  number of segments given in p 
       type(epPos)::  p(*)  ! output.   position given in local coord.

!      
      integer uscl

      call epseeUnderScore(comp%struc, uscl)
      if(comp%struc(1:uscl) .eq. 'octagon') then
         call epDraw_octagon(comp, p, n)
      elseif(comp%struc(1:uscl) .eq. 'honeycomb') then
         call epDraw_honeycomb(comp, p, n)
      elseif(comp%struc(1:uscl) .eq. 'horse') then
         call epDraw_horse(comp, p, n)
      elseif(comp%struc(1:uscl) .eq. 'oecyl') then
         call epDraw_oecyl(comp, p, n)
      elseif(comp%struc(1:uscl) .eq. 'ocyl') then
         call epDraw_ocyl(comp, p, n)
      elseif(comp%struc(1:uscl) .eq. 'sqpipe') then
         call epDraw_sqpipe(comp, p, n)
      elseif(comp%struc(1:uscl) .eq. 'fpolygon') then
         call epDraw_fpolygon(comp, p, n)
      elseif(comp%struc(1:uscl) .eq. 'sqtccl') then
         call epDraw_sqtccl(comp, p, n)
!---  elseif(comp%struc(1:uscl) .eq. '%2') then
!---     call epDraw_%2(comp, p, n)
      else
         call cerrorMsg(comp%struc(1:8), 1)
         call cerrorMsg(' is not supported: epDrawNew; Maybe', 1)
         call cerrorMsg(
     *    "You havent done mkNewVolume 'your-config-file'",
     *     0)
      endif

      end
