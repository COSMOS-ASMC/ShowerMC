!       This is a dummy when no new- volume is used.
!    This is to manage a new non-default  volume/shape which is not
!    included in Epics.
!   This is the interface to draw new volumes.
      subroutine epDrawNew(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component):: comp  ! input. given component to be drawn
      integer n   !  ouput.  number of segments given in from and to pairs
       type(epPos)::  p(*)  ! output.   position given in local coord.

!      
      if(comp%struc(1:8) .eq. ' ') then
!        call epDraw_%1(comp, p, n)
!---  elseif(comp%struc(1:8) .eq. '%2') then
!---     call epDraw_%2(comp, p, n)
      else
         call cerrorMsg(comp%struc(1:8), 1)
         call cerrorMsg(' is not supported: epDrawNew; Maybe', 1)
         call cerrorMsg(
     *    "You havent done mkNewVolume 'your-config-file'",
     *     0)

      endif

      end

