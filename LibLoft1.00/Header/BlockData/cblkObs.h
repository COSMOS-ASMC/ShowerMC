!     *******************
      block data cblkObs
!     *******************
      implicit none
#include  "Zobs.h"
#include  "Zobsp.h"
          
       data 
     * DepthList /maxNoOfSites*0./ ,
     * HeightList /maxNoOfSites*0./ ,
     * ASDepthList /maxNoOfASSites*0./ ,
     * ASHeightList /maxNoOfASSites*0./ ,
     * DtGMT /8.0/ ,
     * LatitOfSite /30.11/ ,
     * LongitOfSite /90.53/ ,
     * ObsPlane /1/ ,	
     * XaxisFromSouth /361.0/,
     * YearOfGeomag /2000.5/
      end
