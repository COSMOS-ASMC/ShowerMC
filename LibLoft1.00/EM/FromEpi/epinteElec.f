      subroutine epinteElec( pj)
!     electron interaction.   
      use modSetIntInf
      implicit none
#include "Zptcl.h"
      type(ptcl),intent(in):: pj

      character(120):: msg

      if(IntInfArray(ProcessNo)%process .eq. 'brem') then
         call epbrem(pj)
      elseif(IntInfArray(ProcessNo)%process .eq. 'knoc') then
         call epknoc(pj)
      elseif(IntInfArray(ProcessNo)%process .eq. 'hcs' ) then
!           it has been done in  epdoMixedMCS2; 
!            so simply puth current track in stack
!     call eppush(cTrack)
!           ????
      elseif(IntInfArray(ProcessNo)%process .eq. 'anih') then
         call epanih(pj)
      elseif(IntInfArray(ProcessNo)%process .eq. 'sync') then
         call epsync(pj)
      else
         write(msg,
     *        '("process=",a4, " for e is undef.")')
     *          IntInfArray(ProcessNo)%process
         call cerrorMsg(msg, 0)
      endif
      end
