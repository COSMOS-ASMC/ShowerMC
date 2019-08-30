      module modcDCS
! set up DCS related  consts for Cosmos Air 
      use modDCS
      implicit none
!    in the case of pointer, segmentation fault
!  will occure when dealing with making dcs cubic spline
!  coef. 
!  type(DCSconst),pointer:: DCSnega  ! for e-
!    
      type(DCSconst),target,allocatable:: DCSnega(:)  ! for e-
!  type(DCSconst),pointer:: DCSposi  ! for e+
      type(DCSconst),target,allocatable:: DCSposi(:)  ! for e+
      end module modcDCS

      subroutine cmkDCSconst
!  use modXsecMedia
      use modcDCS
      implicit none
#include "ZmediaLoft.h"
      integer::i
      allocate( DCSnega(NoOfMedia) )
      allocate( DCSposi(NoOfMedia) )
      do i = 1, NoOfMedia
         call cDCSconstForMedia(media(i), -1,  DCSnega(i))
         call cDCSconstForMedia(media(i),  1,  DCSposi(i))
      enddo
      end subroutine cmkDCSconst
