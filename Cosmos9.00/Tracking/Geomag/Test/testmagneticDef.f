!c       testing cmagneticDef
!
       implicit none
#include  "Zcode.h"
#include  "Ztrack.h"
 
       type(magfield):: B
       type(track):: aTrack
       type(coord):: dr, dd
       real*8 len, norm
       integer i
 
       call cmkptc(kmuon, 0, 1, aTrack%p)
       aTrack%p%fm%p(4) = 10.
       aTrack%vec%w%r(1) = sqrt(0.1)
       aTrack%vec%w%r(2) = 0.
       aTrack%vec%w%r(3) = sqrt(0.9)
       call ce2p(aTrack)
       
       B%x = 1.
       B%y = 0.
       B%z = 0.
 
       aTrack%pos%xyz%r(1) = 0.
       aTrack%pos%xyz%r(2) = 0. 
       aTrack%pos%xyz%r(3) = 0.
       
       len = .5       ! if there is no drift, len ~ 1/6 of radius is vergy good.  Even len =1/3 is torelable.
!                but if there is  drift, we should be careful about len size.  1/100 is needed.
!                    For cosmos application, however, len = radius/5 ~3 has no probolem.
!     This case is shown in eps files.
 
       do i = 1, 5000
          call cmagneticDef(aTrack, B, len, dr, dd)
       
          aTrack%pos%xyz%r(1) = aTrack%pos%xyz%r(1)  +
     *        aTrack%vec%w%r(1) *len   + dr%r(1)
          aTrack%pos%xyz%r(2) = aTrack%pos%xyz%r(2)  +
     *        aTrack%vec%w%r(2) *len   + dr%r(2)
          aTrack%pos%xyz%r(3) = aTrack%pos%xyz%r(3)  +
     *        aTrack%vec%w%r(3) *len   + dr%r(3)
 
          write(*, *) sngl(aTrack%pos%xyz%r(1)), 
     *             sngl(aTrack%pos%xyz%r(2)), 
     *             sngl(aTrack%pos%xyz%r(3))
 
          aTrack%vec%w%r(1) = aTrack%vec%w%r(1) + dd%r(1)
          aTrack%vec%w%r(2) = aTrack%vec%w%r(2) + dd%r(2)
          aTrack%vec%w%r(3) = aTrack%vec%w%r(3) + dd%r(3)
!           normalize      
          norm = sqrt(aTrack%vec%w%r(1)**2 +
     *            aTrack%vec%w%r(2)**2 +
     *            aTrack%vec%w%r(3)**2)
          aTrack%vec%w%r(1) = aTrack%vec%w%r(1)/norm
          aTrack%vec%w%r(2) = aTrack%vec%w%r(2)/norm
          aTrack%vec%w%r(3) = aTrack%vec%w%r(3)/norm
!              reset momentum
          call cresetMom(aTrack)
       enddo
       end
 
