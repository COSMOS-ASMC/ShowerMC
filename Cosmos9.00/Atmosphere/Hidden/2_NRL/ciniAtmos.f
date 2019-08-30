      subroutine ciniAtmos
      implicit none
#include "Zmanagerp.h"
#include "Zobs.h"
#include "Zobsp.h"
      integer:: atmosmodel
      call creadAtmosD  ! use NRL
!     check Lat,Longit
      if( AtmosFile /=  " " ) then
!     If  Lat,Long differ from those in AtmosFile, stop
         call cNRLLatLongCheck( LatitOfSite, LongitOfSite )
      endif   
      call cqAtmosModel(atmosmodel)
      write(0,*) ' Atmosphere model # is ', atmosmodel
      end   subroutine ciniAtmos
      
