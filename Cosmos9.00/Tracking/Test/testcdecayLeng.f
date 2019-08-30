#include  "BlockData/cblkdedx.h"
      implicit none
#include  "Zcondc.h"
#include  "Ztrack.h"
#include  "Zmanagerp.h"

      real*8  length,  E, z
      integer i, imax

      type(track)::aTrack

      write(0,*) ' Enter Energy and heigth, imax'
      z = 10.d3
      E = 5
      imax=100000
      read(*, *) E, z, imax
      call cdedxEleci(200.d-6, .true.)
#if ATMOSPHERE == 1
      AtmosFile = ' '
!          read segmented atmosphere data
      call creadAtmosD
!          manipulate data
      call catmosCnst1
      call catmosCnst2
#endif

      call cmkptc(3, -1, 1, aTrack%p)
      aTrack%p%fm%p(4) = E
      aTrack%p%fm%p(1) = 0.
      aTrack%p%fm%p(2) = 0.
      aTrack%p%fm%p(3) = sqrt(E**2-aTrack%p%mass**2)
      aTrack%pos%height = z

      do i=1, imax
         call cdecayLeng(aTrack, length)
         write(*,*) sngl(length)
      enddo
      end


