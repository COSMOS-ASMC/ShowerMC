      subroutine ceCent2llh(xyz, llh)
!       xyz:  /coord/ structue. Input. The coordinate system is
!            such that the origin is at  the center of the earth
!            x-axis is directed to (0, 0) latitude and longitude.
!            y-axis is directed to (0, 90) latitude and longitude.
!            z-axsi is directed to the north pole.
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
      type(coord)::xyz
      type(coord)::temp
!
      real*8 n, cosphi, sinphi, cosl, sinl
#include  "Zcoordtype.h"

!            ecentricity 0 approximation  ??
           n = Eradius
#ifdef UNIONMAP
           temp%h = sqrt( xyz%r(1)**2 + xyz%r(2)**2 +xyz%r(3)**2)
     *                 - n
           cosphi = sqrt( (xyz%r(1)/(n+temp%h))**2 +
     *                      (xyz%r(2)/(n+temp%h))**2 ) 
!           sinphi = xyz.r(3)/(n*Eone_ecen2 + temp.h)
           sinphi = xyz%r(3)/(n + temp%h)
!      
           cosl = xyz%r(1)/(n+temp%h)/cosphi
           sinl = xyz%r(2)/(n+temp%h)/cosphi
           temp%lat = atan2(sinphi, cosphi)*Todeg
           temp%long = atan2(sinl, cosl)*Todeg
           llh%lat = temp%lat
           llh%long = temp%long
           llh%h = temp%h
#else
           temp%r(3) = sqrt( xyz%r(1)**2 + xyz%r(2)**2 +xyz%r(3)**2)
     *                 - n
           cosphi = sqrt( (xyz%r(1)/(n+temp%r(3)))**2 +
     *                      (xyz%r(2)/(n+temp%r(3)))**2 ) 
!           sinphi = xyz.r(3)/(n*Eone_ecen2 + temp.r(3))
           sinphi = xyz%r(3)/(n + temp%r(3))
!      
           cosl = xyz%r(1)/(n+temp%r(3))/cosphi
           sinl = xyz%r(2)/(n+temp%r(3))/cosphi
           temp%r(1) = atan2(sinphi, cosphi)*Todeg
           temp%r(2) = atan2(sinl, cosl)*Todeg
           llh%r(1) = temp%r(1)
           llh%r(2) = temp%r(2)
           llh%r(3) = temp%r(3)
#endif
           llh%sys = coord_types(2)
       end
