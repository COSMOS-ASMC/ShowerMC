!        to compute muon numbers
      subroutine cgetNmu(eth, nmu)
!        this is not yet made.  tentatively nmu = 0 is given.
      implicit none
#include "Ztrack.h"
! #include "Zmagfield.h"      
#include "Zobs.h"
#include "Zobsv.h"
#include "Zelemagp.h"
      real*8 eth  !  input.  Threshold energy of muons. (GeV) 
      real*8 nmu(maxNoOfASSites)  ! output. number of muons E>eth
      integer l
      do   l=1, NoOfASSites
         nmu(l) = 0.
      enddo
      end
