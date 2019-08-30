      subroutine epLightchkE( aLight, icon)
        use modepLightPty
        implicit none
#include "ZepTrackv.h"
      
       type(epTrack)::  aLight  ! input
       !      check wave length and if outside of the range
       !      make icon = 1
        integer,intent(out)::icon ! if outside of the range 1, else 0
      
!//////////
!        character*16 matter
!
!        if(cPtyNo <= 0 ) then
!           write(0,*) ' in chkE , cPtyNo<=0, Cn=',Cn
!           call epqmat(Cn, matter)
!           write(0,*) ' media=',matter
!           call epqstruc(Cn, matter)
!           write(0,*) ' struc=',matter
!        endif
!///////
        if( cPtyNo <= 0 ) then  !  'sp' case comes here
           icon =0    
!        elseif( aLight.p.subcode > 10 ) then
!           ! starter or stopper; so no check
!           icon = 0
        elseif( comInfo( cPtyNo)%refracIndex(1) == 0.) then
           !  this should be sensor so don't check w.l
           icon =0
        else
           if( aLight%wl < comInfo( cPtyNo)%minmaxWL(1)  .or. 
     *         aLight%wl > comInfo( cPtyNo )%minmaxWL(2)  ) then
              icon  = 1
           else
              icon = 0
           endif
        endif
      end subroutine epLightchkE
