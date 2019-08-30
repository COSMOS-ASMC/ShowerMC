      subroutine epSetEmin(cmpn)
!     Set Emin in ZepTrackv.h
!      From compnent #, get Emin for that compnent.
!      can be used after reading config.  typially called
!      from epLightNewcomp.f 
      use epModify
! >>>>>>>>>>>>>>>>>>>>>>>light
!      use modepLightPty
! <<<<<<<<<<<<<<<<<<<<<light
      implicit none

#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"


      integer:: cmpn      ! component#
      integer::modi ! modifier index
!         when a particle moves  to a new component
!         the following must be done

      EminGamma =  Det%cmp(cmpn)%EminG
      EminElec  =  Det%cmp(cmpn)%EminE
      RecoilKEmin = Det%cmp(cmpn)%RecoilE
      call epSetTcut(Media(MediaNo),  RecoilKEmin ) ! v9.154
      KEmin =  KEminsave
      EminH = Enminsave
      ImperativeEmin = .false.
      if( Det%cmp(cmpn)%modifier > 0 ) then
         modi =Det%cmp(cmpn)%modifier
         if( allocated ( modify ) ) then
            if( IBITS( modify(modi)%kind, bitEmin, 1) > 0) then
               KEmin = modify(modi)%Em%KEmin
               EminH = modify(modi)%Em%Enmin
               ImperativeEmin = modify(modi)%Em%imperative
!            else
!               KEmin = KEminsave
!               EminH = Enminsave
            endif
!         else
!            KEmin = KEminsave
!            EminH = Enminsave
         endif
      endif
      end

