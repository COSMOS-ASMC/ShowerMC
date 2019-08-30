      subroutine cprintCoord(a)
      implicit none
#include  "Zcoord.h"


      type(coord)::a
!
      character*70 msg

      integer i
#include  "Zcoordtype.h"
!
      if(a%sys .eq. coord_types(1)) then
         write(msg, *) 'x =',a%r(1),' y=', a%r(2), ' z=',a%r(3)
      elseif(a%sys .eq. coord_types(2)) then
!         write(msg, *) ' lat(deg) =',a.lat,' longi(deg)=', a.long,
!     *     ' height=',a.h
         write(msg, *) ' lat(deg) =',a%r(1),' longi(deg)=', a%r(2),
     *     ' height=',a%r(3)
      elseif(a%sys .eq.  coord_types(3)) then
!         write(msg, *) ' theta(deg) =',a.theta,' phi(deg)=', a.phi,
!     *     ' radius=',a.radius
         write(msg, *) ' theta(deg) =',a%r(1),' phi(deg)=', a%r(2),
     *     ' radius=',a%r(3)
      else
         do i = 1, max_coord_types
            if(a%sys .eq. coord_types(i)) then
               write(msg, *) a%r(1), a%r(2), a%r(3), ' sys=', a%sys
               goto 10
            endif
         enddo
         write(msg, *) ' coord. system=', a%sys,
     *    ' unknown to print_coord'
         call cerrorMsg(msg, 0)
      endif
 10   continue
      call cerrorMsg(msg, 1)
      end
