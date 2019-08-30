      subroutine cqUHookr(i, rv)
      implicit none
#include "Zmanagerp.h"
!           returns i-th userhook real variable value r,
      integer i  ! input.
      real*8  rv  ! output.  r
      if(i .le. 0 .or. i .gt. MAX_USERHOOKR) then
!         call cerrorMsg('out of range for cqUHookr', 0)
         rv = -1.d-60
      else
         rv = UserHookr(i)
      endif
      end
      subroutine cqUHooki(i, iv)
      implicit none
#include "Zmanagerp.h"
!           returns i-th userhook integer variable value iv
      integer i  ! input.
      integer  iv  ! output.  integer value
      if(i .le. 0 .or. i .gt. MAX_USERHOOKI) then
!         call cerrorMsg('out of range for cqUHooki', 0)
         iv = -9999999
      else
         iv = UserHooki(i)
      endif
      end
      subroutine cqUHookc(i, cv)
      implicit none
#include "Zmanagerp.h"
!           returns i-th userhook character variable value cv
      integer i  ! input.
      character*(*) cv  ! output. blank tail of UserHookc is omitted. 

      integer klena
      if(i .le. 0 .or. i .gt. MAX_USERHOOKC) then
!         call cerrorMsg('out of range for cqUHookc', 0)
         cv = ' '
      else
         if( klena(UserHookc(i)) .gt. 0) then
            cv = UserHookc(i)(1:klena(UserHookc(i)))
         else
            cv = ' ' 
         endif
      endif
      end
      

