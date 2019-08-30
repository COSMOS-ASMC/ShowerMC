#include "Zcondc.h"
      subroutine ciniAtmos
      implicit none
#include "Zmanagerp.h"
#include "Zobs.h"
#include "Zobsp.h"
      integer:: atmosmodel
      
#if ATMOSPHERE == 1
!          read segmented atmosphere data
      call creadAtmosD
!          manipulate data
      call catmosCnst1
      call catmosCnst2
#elif ATMOSPHERE == 2
!          read segmented atmosphere data
      call creadAtmosD
!          manipulate data
      call catmosCnst1
#elif ATMOSPHERE == 3
      call creadAtmosD  ! use NRL
!     check Lat,Longit
      if( AtmosFile /=  " " ) then
!     If  Lat,Long differ from those in AtmosFile, stop
         call cNRLLatLongCheck( LatitOfSite, LongitOfSite )
      endif   
#endif
      call cqAtmosModel(atmosmodel)
      write(0,*) ' Atmosphere model # is ', atmosmodel
      end   subroutine ciniAtmos
      
