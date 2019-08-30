!           test csPrimAng
!
!
!       sPrimAng < ../../Data/Namelist/parameters
!
      implicit none
#include  "Ztrack.h"
#include  "Zincidentp.h"
#include  "Zincidentv.h"

#include "Zglobalc.h"
#include "Zobs.h"      
#include "Zobsp.h"      
#include "Zobsv.h"      
      external cblkObs
      integer i
!
      type(coord)::dir

      call creadParam(5)
      call cinitObs
      write(*, *) ' Latitude=', sngl(LatitOfSite), 
     *  ' deg. Longitude=', sngl(LongitOfSite), 
     *  ' deg. DtGMT=', sngl(DtGMT), ' hours',
     *  ' year for Geomagnetism=', sngl(YearOfGeomag)

      Coszenith=(0.9,1)
      SourceDec= 30.
      Ddelta = 2.5
      Za1ry ='aps'
      call ciniSPrimAng

      write(*, *) ' Za1ry=',Za1ry, ' Obsvhour=',Obsvhour,
     *   ' Sourcedec=',SourceDec
      write(*, *) 'Coszenith=',Coszenith

      do i = 1, 10
         call csPrimAng(dir)
         write(*, *) dir%r(1), dir%r(2), dir%r(3)
      enddo
      end
