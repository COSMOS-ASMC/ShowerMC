!     at execution, give namelist data by redirection. e.g.,
!       initObs < ../../Data/Namelist/parameters
!
      program testInitObs
      implicit none

!
      call creadParam(5)
      call cinitObs
      call cprintObs
      end
#include "BlockData/cblkElemag.h"
