      subroutine cprintMagF(b)
!         print magnetic field for debug purpose
!       b: /magfield/   input.
!
      implicit none

#include  "Zmagfield.h"
      type(magfield)::b
      character*100 msg
#ifdef UNIONMAP
!
      if(b%sys .eq. 'xyz') then
         write(*, *) 'mag_field: x=',b%x, ' y=',b%y,
     *   ' z=',b%z
      elseif(b%sys .eq. 'hva') then
         write(*,*) ' mag_field: horizontal=',b%h,
     *   ' vertical=', b%v, ' angle=', b%a,
     *   ' deg'
      elseif(b%sys .eq. 'ned') then
         write(*,*) ' mag_field: north=',b%n,
     *   ' east=',b%e, ' down=',b%d
      else
         write(msg,*) ' mag coord type=',b%sys,
     *  ' unkown to cprint_magfield',
     *   ' componen=', b%x, b%y, b%z
         call cerrorMsg(msg, 0)
      endif
#else
      if(b%sys .eq. 'xyz') then
         write(*, *) 'mag_field: x=',b%x, ' y=',b%y,
     *   ' z=',b%z
      elseif(b%sys .eq. 'hva') then
         write(*,*) ' mag_field: horizontal=',b%x,
     *   ' vertical=', b%y, ' angle=', b%z,
     *   ' deg'
      elseif(b%sys .eq. 'ned') then
         write(*,*) ' mag_field: north=',b%x,
     *   ' east=',b%y, ' down=',b%z
      else
         write(msg,*) ' mag coord type=',b%sys,
     *  ' unkown to cprint_magfield',
     *   ' componen=', b%x, b%y, b%z
         call cerrorMsg(msg, 0)
      endif
#endif
      write(*,*) ' Unit is T. Multiply 10**4 to Gausss'
      end
