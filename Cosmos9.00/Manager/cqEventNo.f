!      inquire the event number
      subroutine cqEventNo(num, cumnum)
      implicit none
#include "Zmanagerp.h"
      integer num   ! output. number of events in the current run
      integer cumnum ! output. cummulative number of events so far.
!     may be used  after the initialization of an event, then
!     this gives the number for that event.
!

      num = EventsInTheRun
      cumnum = EventNo
      end
!        added in v7.642; may be used form chookBgEvent
!       when the user resets the primary by reading a file
!       which stores 1ry information,  but the file is 
!       exhausted.  (see Readme in Util/Atmnc3/Interface/)
      subroutine cresetDestEvNo(num)
      implicit none
#include "Zmanagerp.h"
      integer,intent(in):: num
      DestEventNo(2) = num
      end
