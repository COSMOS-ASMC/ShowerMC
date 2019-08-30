!>brief set Calculated magnetic field to /magfield/ b
!>@param[in] sys which system. 'xyz', 'hva', 'ned' etc data type character*3.
!>@param[in] b1,b2,b3 3 components of mag. in the system 'sys' data type real*8
!>@param[out] b  /magfield/ 
      subroutine csetMagField(sys, b1, b2, b3, b)
!         set magnetic field strength in b
!   sys: character*3. input.  which system. 'xyz', 'hva', 'ned' etc
!  b1,b2,b3:real*8 input.  3 components of mag. in the system 'sys'
!     b: /magfield/ .  output. 
!
      implicit none
#include  "Zmagfield.h"
      character*(*) sys
      real*8 b1, b2, b3
      type(magfield)::b
!
      character*70 msg
       integer i
#include "Zmagfieldtp.h"

       do i=1, max_mag_types
         if(sys .eq. mag_types(i)) then
              b%x = b1
              b%y = b2
              b%z = b3
              b%sys = sys
              goto 10
          endif
       enddo
       write(msg, *) ' error sys=',sys,'to csetMagField'
       call cerrorMsg(msg, 0)
  10   continue
      end
