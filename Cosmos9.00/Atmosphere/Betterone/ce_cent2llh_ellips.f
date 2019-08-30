c      include 'cllh2eCent.f'
!     program test_ceCent2llh
!
c      include '../Zcoord.h'
!      record /coord/ a, b
!      integer ios
!      real*8 x, y, z 
!      do while (.true.)
!         read(*, *,iostat=ios) x, y, z
!         if(ios .ne.. 0) goto 10
!         call csetCoord('llh', x, y, z, a)
!         call cconv_coord_to('ecen', a, b)
!         call cconv_coord_to('llh',  b, a)
!         write(*,*) a.lat, a.long, a.h
!      enddo   
! 10   continue
!      end
!    earth center rectangular coordinate to longitude latitude height 
!    coordinate  converstion.
!
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
      implicit none
c----      include '../../Zglobalc.h'
#include  "Zglobalc.h"
c----      include '../Zcoord.h'
#include  "Zcoord.h"
c----      include 'Zearth.h'
#include  "Zearth.h"
      type(coord)::llh,  xyz
      type(coord):: temp
!
      real*8 n, cosphi, sinphi, cosl, sinl, hold
      integer i
!        start with ecentricity 0 approximation
         temp%h = 0.
         n = Eradius
         do i = 1, 6
            hold = temp%h
            temp%h = sqrt( xyz%r(1)**2 + xyz%r(2)**2 +(xyz%r(3)/
     *      (1.d0 - n*ecen2/(n+temp%h)))**2 ) - n
             cosphi = sqrt( (xyz%r(1)/(n+temp%h))**2 +
     *                      (xyz%r(2)/(n+temp%h))**2 ) 
             sinphi = xyz%r(3)/(n*one_ecen2 + temp%h)
!      
             cosl = xyz%r(1)/(n+temp%h)/cosphi
             sinl = xyz%r(2)/(n+temp%h)/cosphi
             n = Eradius /sqrt(1.d0 - ecen2*sinphi**2)
            if(abs(temp%h -hold) .lt. 1.d-2) goto 100   ! abs error < 1 mm 
         enddo 
 100     continue
!        write(*, *) ' i=', i  !  <i> = 4 
!
         temp%lat = atan2(sinphi, cosphi)*Todeg
         temp%long = atan2(sinl, cosl)*Todeg
         llh%lat = temp%lat
         llh%long = temp%long
         llh%h = temp%h
         llh%sys = coord_types(2)
       end

