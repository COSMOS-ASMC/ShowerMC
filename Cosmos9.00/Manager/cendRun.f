!         to be called when the end of execution status is reached.
      subroutine cendRun
      implicit none
#if defined (KEKA) || defined (KEKB)
#include "mpif.h"
      integer err
#endif          

!           save current status.      
          call csaveStatus 
!            user hook
          call chookEnRun
#if defined (KEKA) || defined (KEKB)
          call mpi_finalize(err)
#endif          
      end
!       ***************************************** stave  cont info
      subroutine csaveStatus
      implicit none
#include "Zmanagerp.h"
      
      integer icon

      PrevEventNo = EventNo    ! current event no.
      call rnd1s(InitRN)  !  current seed.
      if(ContFile .ne. ' ' ) then
         call copenNLfw(TempDev, ContFile, icon)
!      open(TempDev, file=ContFile, form='formatted')
         call cwriteParam(TempDev, 1)
         close(TempDev)
      endif
      end
!       ***************************************** print final status
      subroutine cprintStatus
      implicit none
#include "Zmanagerp.h"
#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryv.h"


      character*80 msg
      integer i

      write(msg, *) ' No of cummulative events =', PrevEventNo, 
     *      ' No of events in this run=', EventsInTheRun
      call cerrorMsg(msg, 1)
      write(msg, *) ' comp.    sampled    accepted'
      call cerrorMsg(msg, 1)
      do i = 1, Prim%no_of_comps
         write(msg, '(i4,2i12)')
     *   i, Prim%NoOfSampComp(i,1), Prim%NoOfSampComp(i,2)
         call cerrorMsg(msg, 1)
      enddo
      call cerrorMsg("###end of run###", 1)
      end
