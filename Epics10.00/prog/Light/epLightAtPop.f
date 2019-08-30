      subroutine epLightAtPop(aTrack, icon )
!       when a light is poped up, this examines something
      use modepLightPty
      implicit none
#include  "ZepTrackv.h"
#include  "Zcode.h"

       type(epTrack)::  aTrack  ! input track info. ( klight, kEdepp,
                               ! kchgPath track)
      integer,intent(out):: icon   ! 0--> need further trakcing 
                                   ! 1--> not need to further tracking
                                   !  so pop another track 

      if( aTrack%p%code == klight ) then
         icon = 0               ! need further tracking
      elseif( aTrack%p%code == kchgPath ) then
         call epLightCerenkov(aTrack)
         icon = 1
      elseif( aTrack%p%code == kEdepo ) then
         call epLightScintiFromCell( aTrack )
         icon  = 1
      else
         write(0,*) ' strange code =',aTrack%p%code, 
     *         ' in  epLightAtPop '
         stop
      endif
      end
