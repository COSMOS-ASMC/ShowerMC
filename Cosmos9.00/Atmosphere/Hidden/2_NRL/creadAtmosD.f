      subroutine creadAtmosD
      implicit none
#include "Zmanagerp.h"
#include "Zobs.h" 
#include "Zobsp.h" 
      if( AtmosFile == ' ' ) then
         if( any( NRL_period(1:2) <= 0 ) .or.
     *       any( NRL_period(1:2) > 365) .or. 
     *       any( NRL_period(3:4) < 0  ) .or.
     *       any( NRL_period(3:4) > 23 )  ) then
            write(0,*)
     *      ' AtmosModel-2 is being used, i.e'
            write(0,*)
     *      ' you want to use NRL time-place-dependent atmosphere model'
            write(0,*)
     *      ' but "AtmosFile" is blank. In this case, NRL_period(1:2)'
            write(0,*)
     *      ' must be 1~365, and (3:4)  0-23 but not so ',
     *      ' NRL_peirod(:)=',NRL_period(:)
            stop
         endif
         call cNRLGenData(LatitOfSite, LongitOfSite, NRL_period)
      else
         call cNRLdataRead(TempDev, AtmosFile)
         call cNRLdataManip
      endif
      end subroutine creadAtmosD

      subroutine cqAtmosModel(modelno)
      implicit none
      integer,intent(out):: modelno
      modelno = 2
      end subroutine cqAtmosModel
