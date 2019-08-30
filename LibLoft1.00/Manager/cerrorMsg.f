!        printout error message on error out.
!
      subroutine cerrorMsg(msg, needrtn)
      implicit none
#include  "Zmanagerp.h"

      external cblkManager

      character*(*) msg !  input.  message to be printed
      integer  needrtn  !  input. 0  --> stop 9999,else control is returned
      integer klena
      real(8):: dummy

      if(klena(msg) .ge. 1) then
         write(ErrorOut, '(a)') msg(1:klena(msg))
      endif
      if(needrtn .eq. 0)  then
#if defined (PCLinuxIFC64) || (PCLinuxIFC) || (MacIFC)
         call TRACEBACKQQ()
#endif
         call cbackTrace(dummy)
      endif
      end
      subroutine cbackTrace(a, b, c)
!         This is to print back trace information
!       when some error is detected by the user
!       and the user want to see the routine names
!       before the error is detected. (compilation
!       with backtrace mode is needed)
!       the calling program should not contain
!       good "a b c". but may be better to be
!         call cbactTrace(dummy)
!
       implicit none
      real(8),intent(out)::a, b, c(1000000)
      
      c(:) =  0.
      b = 100.
      a = 100.
      end

      
