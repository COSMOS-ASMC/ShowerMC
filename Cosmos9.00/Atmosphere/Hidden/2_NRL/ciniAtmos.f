      subroutine ciniAtmos
      implicit none
#include "Zmanagerp.h"
#include "Zobs.h"
#include "Zobsp.h"

      call creadAtmosD  ! use NRL
!     check Lat,Longit
      if( AtmosFile /=  " " ) then
!     If  Lat,Long differ from those in AtmosFile, stop
         call cNRLLatLongCheck( LatitOfSite, LongitOfSite )
      endif   
      call cqAtmosModel(AtmosModel)
      write(0,*) ' Atmosphere model # is ', AtmosModel
      end   subroutine ciniAtmos
      
