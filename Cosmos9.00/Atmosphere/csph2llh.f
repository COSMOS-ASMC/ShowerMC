      subroutine csph2llh(sph, llh)
!       sph:  /coord/ structue. Input. spherical coordinate system
!       llh:  /coord/ structure. Output. to contain data in latitude,
!                                       longitude, height.
!    ***  note  ***
!           llh can be the same as xyz.  time component is unchanged
!
      use modAtmosDef
      implicit none

#include  "Zglobalc.h"

#include  "Zcoord.h"
      type(coord)::llh
      type(coord)::sph
      type(coord)::temp
#include  "Zcoordtype.h"

!         ecentricity 0 approximation
#ifdef UNIONMAP
           temp%h = sph%radius - Eradius
           temp%lat = 90.d0 - sph%theta
!
           if(sph%phi .gt. 180.d0) then
              temp%long = sph%phi - 360.d0
           else
              temp%long = sph%phi
           endif
#else
           temp%r(3) = sph%r(3) - Eradius
           temp%r(1) = 90.d0 - sph%r(1)
!
           if(sph%r(2) .gt. 180.d0) then
              temp%r(2) = sph%r(2) - 360.d0
           else
              temp%r(2) = sph%r(2)
           endif
#endif
           temp%sys = coord_types(2)
           llh = temp
       end
