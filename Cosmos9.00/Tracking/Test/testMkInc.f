!     at execution, give namelist data by redirection. e.g.,
!      mkInc < ../../Data/Namelist/parameters
!
      program testMkInc
      implicit none
#include "Ztrack.h"

      type(track)::anIncident

!
      call creadParam(5)
      call cinitObs
      call cmkIncident(anIncident)
      end

