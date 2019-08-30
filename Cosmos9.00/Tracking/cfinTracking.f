      subroutine cfinTracking
      implicit none
#include "Ztrack.h"
! #include "Zmagfield.h"            
#include "Ztrackp.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"

      integer i, ib, ie, is
      logical dvlpd

      if( (Trace .gt. 0 .and. Trace .lt. 60) .or. Trace .gt. 100) then
         close(TraceDev)    ! close default trace file  for 1 event.
      elseif(Trace .gt. 60  .and. Trace .lt. 100) then
         call cputCerenkovE     ! cerenkov output trailer
      endif
!                                               for each event
!         get size weighted age
      dvlpd = .false.
      if(ObserveAS) then
         if(Upgoing) then   ! upgoing primary
            ib = NoOfASSites
            ie = 1
            is = -1
         else
            ib = 1
            ie = NoOfASSites
            is = 1
         endif
         do i = ib, ie, is
            if(ASObsSites(i)%esize .gt. 0.) then
               dvlpd = .true.
               ASObsSites(i)%age = ASObsSites(i)%age/ASObsSites(i)%esize
            elseif(dvlpd) then
               ASObsSites(i)%age = 5.0    ! died out.
            else
               ASObsSites(i)%age = 0.0   ! not yet developed
            endif
         enddo
      endif
      end
