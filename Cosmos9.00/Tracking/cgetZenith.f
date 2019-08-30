#include "ZsubstRec.h"
!    
      subroutine csetDirCos(dc, aTrack)
      implicit none

#include  "Ztrack.h"
      type(track)::aTrack
      type(coord)::dc
!
#ifdef SUBSTREC      
      aTrack%vec%w = dc
#else
      aTrack%vec%w%r = dc%r
#endif
      call cgetZenith(aTrack, aTrack%vec%coszenith)
      end
!         get cos of zenith angle from track information
!
      subroutine cgetZenith(aTrack, cosz)
      implicit none

#include  "Ztrack.h"
      type(track)::aTrack
      real*8 cosz

!
      cosz = - (aTrack%vec%w%r(1) * aTrack%pos%xyz%r(1) + 
     *          aTrack%vec%w%r(2) * aTrack%pos%xyz%r(2) + 
     *          aTrack%vec%w%r(3) * aTrack%pos%xyz%r(3))/
     *          aTrack%pos%radiallen
      if(cosz .gt. 1.d0) then
         cosz = 1.d0
      elseif(cosz .lt. -1.d0) then
         cosz = -1.d0
      endif
           ! or
!     *     sqrt(aTrack.pos.xyz.r(1)**2 +
!     *          aTrack.pos.xyz.r(2)**2 +
!     *          aTrack.pos.xyz.r(3)**2)
      end
