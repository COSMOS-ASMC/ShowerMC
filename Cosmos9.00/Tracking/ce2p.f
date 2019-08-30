!        ce2p:  compute px, py, pz from energy with a given
!              direction vector and mass.
!
      subroutine ce2p(aTrack)
      implicit none
!----      include 'Ztrack.h'
#include  "Ztrack.h"
      type(track)::aTrack
!
      real*8 p
      p  = aTrack%p%fm%p(4)**2 - aTrack%p%mass**2
      if(p .lt. 0.) then
!         write(*,*) ' p=',p,' code =',aTrack.p.code
!            since |p| is order of  10^-12        
         p = 0.
      endif
      p = sqrt(p)
!      p  = sqrt( aTrack.p.fm.p(4)**2 - aTrack.p.mass**2 )
      aTrack%p%fm%p(1) = p * aTrack%vec%w%r(1)
      aTrack%p%fm%p(2) = p * aTrack%vec%w%r(2)
      aTrack%p%fm%p(3) = p * aTrack%vec%w%r(3)
      end
