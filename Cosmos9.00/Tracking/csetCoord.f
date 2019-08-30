      subroutine csetCoord(sys, x, y, z, a)
!          set coordinate values 
!     x,y,z: real*8. input. 3 components of the coordinate.
!                           meaning depends on sys.
!         a: /coord/. output. coordinate values are set. time is
!                           unchanged.
!       sys: character*4. input.  what coordinate.
!
      use modAtmosDef
      implicit none
#include  "Zcoord.h"
! #include  "Zearth.h"

      real*8 x, y, z
      type(coord)::a
      character*(*) sys
      integer i
!
      character*70 msg
#include  "Zcoordtype.h"
!
      a%r(1) = x
      a%r(2) = y
      a%r(3) = z
      a%sys = sys
      do i = 1, max_coord_types
         if(sys .eq. coord_types(i)) goto 10
      enddo
      write(msg, *) ' coordinate type=', sys, ' invalid',
     *            ' input to csetCoord'
      call cerrorMsg(msg, 0) 
 10   continue
      end
