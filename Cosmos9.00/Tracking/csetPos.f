!    csetPos:  use only xyz.
!    csetPos2:  use input height and depth.
!
!       set position information for a given xyz
!       location: input/output.    coord part of location is input.
!     
      subroutine csetPos(location)
!
      use modAtmosDef
      implicit none
#include  "Zobs.h"
#include  "Zcoord.h"
#include  "Zpos.h"
! #include  "Zearth.h"
!
      type(position)::location
      real*8 cvh2thick
      character*70 msg 
!
      if(location%xyz%sys .ne. 'xyz') then
         write(msg, *) ' error coord. sys to csetPos=', location%xyz%sys
         call cerrorMsg(msg, 0)
      endif
      location%radiallen = 
     *      sqrt(location%xyz%r(1)**2 + location%xyz%r(2)**2 +
     *      location%xyz%r(3)**2)      
      location%height = location%radiallen - Eradius
      location%depth = cvh2thick(location%height)
      end
!     --------------------------------------
      subroutine csetPos2(h, d, location)
!     --------------------------------------
!
      use modAtmosDef
      implicit none
!----      include 'Zobs.h'
#include  "Zobs.h"
!----      include 'Zcoord.h'
#include  "Zcoord.h"
#include  "Zpos.h"
!----      include './Atmosphere/Zearth.h'
!#include  "Zearth.h"
!
      real*8 h, d
      type(position)::location
!
      location%height = h
      location%radiallen = h + Eradius
      location%depth = d
      end
