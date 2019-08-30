      implicit none
#include "Zcode.h"
#include "Ztrack.h"

      type(track)::aTrack

      type(coord)::dir
      type(coord)::B
      type(coord)::dr
      real*8 p, sum, Babs, r, len

      logical old
!************ for old ******
      real*8 norm
!*************
      integer i, j

      old = .false.
      read(*,*) old
      aTrack%vec%w%x = 0.
      aTrack%vec%w%y = sqrt(2.d0)/2
      aTrack%vec%w%z = -sqrt(2.d0)/2
      call cmkptc(knuc, 1, -1, aTrack%p)

      p = 1.d0
      aTrack%p%fm%p(1) = aTrack%vec%w%x*p
      aTrack%p%fm%p(2) = aTrack%vec%w%y*p
      aTrack%p%fm%p(3) = aTrack%vec%w%z*p
      sum =0.d0
      do i =1, 3
         sum = sum + aTrack%p%fm%p(i)**2
      enddo
      aTrack%p%fm%p(4) =sqrt( aTrack%p%mass**2 + sum  )
      aTrack%pos%xyz%x = 0.
      aTrack%pos%xyz%y = 0.
      aTrack%pos%xyz%z = 0.
      
      Babs = 0.3d-4
      B%x = 0.
      B%y = sqrt(2.d0)/2 * Babs
      B%z = sqrt(2.d0)/2 * Babs

      call cmagDefR(aTrack, B, r)
!////////////
!      write(*,*) r
!/////////
      len = r/10.
      do i = 1, 1000
         write(*,*)
     *         sngl(aTrack%pos%xyz%x),
     *         sngl(aTrack%pos%xyz%y),
     *         sngl(aTrack%pos%xyz%z)

         call cmagneticDef(aTrack, B, len, dr, dir)
!           for new
         if(.not. old) then
            aTrack%pos%xyz%x = aTrack%pos%xyz%x + dr%x
            aTrack%pos%xyz%y = aTrack%pos%xyz%y + dr%y
            aTrack%pos%xyz%z = aTrack%pos%xyz%z + dr%z
            aTrack%vec%w = dir
         else
!
!
!         for old
            aTrack%pos%xyz%x = aTrack%pos%xyz%x +
     *           len*aTrack%vec%w%x +   dr%x
            aTrack%pos%xyz%y = aTrack%pos%xyz%y +
     *           len*aTrack%vec%w%y +   dr%y
            aTrack%pos%xyz%z = aTrack%pos%xyz%z +
     *           len*aTrack%vec%w%z +   dr%z

            aTrack%vec%w%x = aTrack%vec%w%x + dir%x
            aTrack%vec%w%y = aTrack%vec%w%y + dir%y
            aTrack%vec%w%z = aTrack%vec%w%z + dir%z
!           normalize
            norm = sqrt(aTrack%vec%w%r(1)**2 +
     *           aTrack%vec%w%r(2)**2 +
     *           aTrack%vec%w%r(3)**2)
            do j = 1, 3
               aTrack%vec%w%r(j) = aTrack%vec%w%r(j)/norm
            enddo
         endif
      enddo
      end


         
