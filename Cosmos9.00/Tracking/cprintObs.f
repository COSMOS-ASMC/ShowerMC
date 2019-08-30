      subroutine cprintObs(io)
      implicit none
#include "Zglobalc.h"
#include "Zcoord.h"
#include "Zpos.h"      
#include "Zmagfield.h"
#include "Zobs.h"      
#include "Zobsp.h"      
#include "Zobsv.h"      

      integer io  ! input. Fortran logical dev. no. for output.
      integer i

      write(io, *) ' Latitude=', sngl(LatitOfSite), 
     *  ' deg. Longitude=', sngl(LongitOfSite), 
     *  ' deg. DtGMT=', sngl(DtGMT), ' hours',
     *  ' year for Geomagnetism=', sngl(YearOfGeomag)
      write(io, *)'-----------------------------'
      write(io, *) '  Position of Obs. sites '
      call cprObsSiteHd(io)
      do i = 1, NoOfSites
         call cprObsSite(io, ObsSites(i))
      enddo
      write(io, *) ' -----------------------------'
      write(io, *) ' As Obs. sites'
      call cprASObsSiteHd(io)
      do i = 1, NoOfASSites
         call cprASObsSite(io, ASObsSites(i))
      enddo
      end
