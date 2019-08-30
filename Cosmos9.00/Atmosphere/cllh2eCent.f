      subroutine cllh2eCent(llh, xyz)
!       llh:  /coord/ structure. input. containing data in latitude,
!                                       longitude, height.
!       xyz:  /coord/ structue. output. The coordinate system is
!            such that the origin is at  the center of the earth
!            x-axis is directed to (0, 0) latitude and longitude.
!            y-axis is directed to (0, 90) latitude and longitude.
!            z-axsi is directed to the north pole.
!
!           xyz.r(1) is  x coordinate value in m 
!           xyz.r(2) is  y                    
!           xyz.r(3) is  z                    
!    ***  note  ***
!           xyz can be the same as llh.  time component is unchanged
!
      use modAtmosDef
      implicit none
#include  "Zglobalc.h"

#include  "Zcoord.h"

      type(coord)::llh
      type(coord)::xyz
      type(coord)::temp
!
      real*8 nh
#include  "Zcoordtype.h"
!         ecentricity =0
#ifdef UNIONMAP
      nh = Eradius + llh%h
      temp%x = nh * cos(llh%lat*Torad)* cos(llh%long*Torad)
      temp%y = nh * cos(llh%lat*Torad)* sin(llh%long*Torad)
      xyz%r(3)  = nh * sin(llh%lat * Torad)
      xyz%r(1) = temp%x
      xyz%r(2) = temp%y
#else
      nh = Eradius + llh%r(3)
      temp%r(1) = nh * cos(llh%r(1)*Torad)* cos(llh%r(2)*Torad)
      temp%r(2) = nh * cos(llh%r(1)*Torad)* sin(llh%r(2)*Torad)
      xyz%r(3)  = nh * sin(llh%r(1) * Torad)
      xyz%r(1) = temp%r(1)
      xyz%r(2) = temp%r(2)
#endif
      xyz%sys = coord_types(1)
      end
