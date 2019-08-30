!        call this when momentum has been changed to reset the
!     direction information.
!
      subroutine cresetDirec(aTrack)
      implicit none

#include  "Ztrack.h"
      type(track)::aTrack
!
      real*8 pabs

      call cpxyzp(aTrack%p%fm, pabs)
      if(pabs .gt.  0) then
         aTrack%vec%w%r(1) = aTrack%p%fm%p(1)/pabs
         aTrack%vec%w%r(2) = aTrack%p%fm%p(2)/pabs
         aTrack%vec%w%r(3) = aTrack%p%fm%p(3)/pabs
      else
         aTrack%vec%w%r(1) = 0.
         aTrack%vec%w%r(2) = 0.
         aTrack%vec%w%r(3) = 1.
      endif
      call cgetZenith(aTrack, aTrack%vec%coszenith)
      end


         
