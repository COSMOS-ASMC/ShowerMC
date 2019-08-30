      subroutine cllh2sph(llh, sph)
!       llh:  /coord/ structure. input. lat_long_height
!       sph:  /coord/ structue. output.. spherical coordinate system
!
!    ***  note  ***
!          sph can be the same as llh.  time component is unchanged
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
           temp%radius = llh%h + Eradius
           temp%theta = 90.d0 - llh%lat
!
           if(llh%long .lt. 0.d0) then
              temp%phi = llh%long + 360.d0
           else
              temp%phi = llh%long
           endif
#else
           temp%r(3) = llh%r(3) + Eradius
           temp%r(1) = 90.d0 - llh%r(1)
!
           if(llh%r(2) .lt. 0.d0) then
              temp%r(2) = llh%r(2) + 360.d0
           else
              temp%r(2) = llh%r(2)
           endif
#endif
           temp%sys = coord_types(3)
           sph = temp
       end
