      subroutine epinteG( pj)
!     g interaction.   if photo hadron -> call cfixTargee must have been
!        called      
      use modSetIntInf
      implicit none
#include "Zptcl.h"
      type(ptcl),intent(in):: pj

      character(120):: msg
      
      if(IntInfArray(ProcessNo)%process .eq. 'comp') then
         call epcmpt(pj)
      elseif(IntInfArray(ProcessNo)%process .eq. 'pair') then
         call eppair(pj)
      elseif(IntInfArray(ProcessNo)%process .eq. 'phot') then
         call epphot(pj)
      elseif(IntInfArray(ProcessNo)%process .eq. 'coh') then
         call epcoher(pj)
      elseif(IntInfArray(ProcessNo)%process .eq. 'photop') then
!                     need not now ?
!            call cfixTarget(xmedia(mediumNo)) ! 
!                       
         call cphotop( pj ) 
      elseif(IntInfArray(ProcessNo)%process .eq. 'mpair') then
         call epmpair(pj)
      else
         write(msg,
     *        '('' proccess='',a4,'' for gamma undefined'')')
     *        IntInfArray(ProcessNo)%process
         call cerrorMsg(msg,0)
      endif
      end
