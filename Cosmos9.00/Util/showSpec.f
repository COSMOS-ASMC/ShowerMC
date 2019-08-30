#include  "BlockData/cblkGene.h"
      implicit none
#include "Zglobalc.h"
#include "Zmanagerp.h"
#include "ZrigCut.h"
#include "Zptcl.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zprimary.h"
#include "Zprimaryc.h"
#include "Zprimaryv.h"
#include "Zincidentp.h"
#ifdef NEXT486
#define IMAG_P dimag
#else
#define IMAG_P imag
#endif


      integer  i, j

      real*8 cosmin, cosmax, rigc
      real*8 flux, E, azmmin, azmmax, dE, zen1, zen2


      call creadParam(5)
      call cbeginRun
      call cprintPrim(ErrorOut)
      cosmax = IMAG_P(CosZenith)
      cosmin = real(CosZenith)
      azmmax = IMAG_P(Azimuth) + XaxisFromSouth
      azmmin = real(Azimuth) + XaxisFromSouth
      dE=10.**0.05
      if(CutOffFile .eq. ' ') then
         rigc = 0.
         call csetCosdeg(.false., .false.)
      else
         rigc = 100.
         if(ZenValue .eq. 'deg') then
            call csetCosdeg(.false., .true.)
         else
            call csetCosdeg(.false., .false.)
         endif
      endif
      zen1 = cosmin
      zen2 = cosmax
      do i = 1, Prim%no_of_comps
         E = Prim%each(i)%energy(1)
         do j = 1, 10000
            call cprimFlux(Prim%each(i), E, rigc, zen1, zen2,
     *      azmmin, azmmax,  flux)
!            if(flux .gt. 0.) then
               write(*,*)
     *         Prim%each(i)%label, sngl(E), sngl(flux)
!            endif
            E = E*dE
            if(E .gt. Prim%each(i)%energy(Prim%each(i)%no_of_seg+1))
     *          goto 100            
         enddo
 100     continue
!         write(*,*)
      enddo
      end
      subroutine chookTrace
      end
      subroutine chookCeren
      end
      subroutine chookCerenS
      end
      subroutine chookCerenE
      end
      subroutine chookBgRun
      end
