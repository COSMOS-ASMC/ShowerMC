      subroutine cresetMom(aTrack)
!      Assum the case where the direction cos has been changed
!      and the directionof  vector p should be  changed.
!        (see cresetMombyEDir)
      implicit none

#include  "Ztrack.h"
      type(track)::aTrack  ! input/output.
!
      real*8 p
      call cpxyzp(aTrack%p%fm, p)
      aTrack%p%fm%p(1) = p * aTrack%vec%w%r(1)
      aTrack%p%fm%p(2) = p * aTrack%vec%w%r(2)
      aTrack%p%fm%p(3) = p * aTrack%vec%w%r(3)

      end       subroutine cresetMom

      subroutine cresetMombyEdir(aTrack)
!          Assume the case where the direction cos is correctly given
!     and E has beeen changed, Vec p should not be used.

!     If E< mass, E is reset to mass. and vec p becomes 0.
      implicit none
#include  "Ztrack.h"
      type(track)::aTrack  ! input/output.
      real(8)::pabs
      pabs = aTrack%p%fm%p(4)**2 - aTrack%p%mass**2
      if( pabs <= 0. ) then
         aTrack%p%fm%p(4) = aTrack%p%mass
         pabs = 0.
      else
         pabs = sqrt(pabs)
      endif
      aTrack%p%fm%p(1:3) = aTrack%vec%w%r(1:3)*pabs
      end      subroutine cresetMombyEdir
