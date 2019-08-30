      module modcTPXS
! set up TPXS related  consts for Cosmos Air
!  use modXSecMedia
      use modTPXS
      implicit none
!  type(TPXSconst),pointer:: TPXSnega  ! for e-
!  type(TPXSconst),pointer:: TPXSposi  ! for e+
      type(TPXSconst),target,allocatable:: TPXSnega(:)  ! for e-
      type(TPXSconst),target,allocatable:: TPXSposi(:)  ! for e+  
      end module modcTPXS

      subroutine cmkTPXSconst
!  use modXsecMedia
      use modcTPXS
      implicit none
#include "ZmediaLoft.h"  
      integer::i
      allocate( TPXSnega(NoOfMedia))
      allocate( TPXSposi(NoOfMedia))
      do i = 1, NoOfMedia
         call cMWconstForMedia(media(i), -1, TPXSnega(i))
         call cMWconstForMedia(media(i),  1, TPXSposi(i))
      enddo
      end subroutine cmkTPXSconst
