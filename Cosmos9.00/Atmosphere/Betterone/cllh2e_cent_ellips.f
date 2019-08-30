!c         test cllh2eCent
c      include '../../Tracking/Zcoord.h'
!      record /coord/ a, b
!      a.lat = 35.
!      a.long =138.
!      a.h =0.
!      do while (a.lat .le. 90.001)
!         read(*, *) a.lat, a.long, a.h
!         call cllh2eCent( a, b)
!         write(*,*) b.x, b.y, b.z
!      enddo   
!      end
!          longitude latitude height to earth center rectangular 
!     coordinate  converstion.
!
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
      implicit none
c----      include '../../Zglobalc.h'
#include  "Zglobalc.h"
c----      include '../../Tracking/Zcoord.h'
#include  "Zcoord.h"
c----      include 'Zearth.h'
#include  "Zearth.h"
      type(coord)::llh,  xyz
      type(coord):: temp
!
      real*8 n
!
      n = Eradius/sqrt(1. -Eecen2 * sin(llh%lat*Torad)**2)
      temp%x = (n + llh%h) *cos(llh%lat*Torad)* cos(llh%long*Torad)
      temp%y = (n + llh%h) *cos(llh%lat*Torad)* sin(llh%long*Torad)
      xyz%r(3)  = (n* Eone_ecen2  + llh%h)*sin(llh%lat * Torad)
      xyz%r(1) = temp%x
      xyz%r(2) = temp%y
      xyz%sys = coord_types(1)
      end
