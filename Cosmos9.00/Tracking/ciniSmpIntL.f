!     ******************************
      subroutine ciniSmpIntL
!       use modMCScontrol
      implicit none

#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackv.h"

!     clear conter
      call ciniIntInf   !    NumberOfInte = 0      ProcessNo = 0
      MoveStat = ToInteract
      call csetActiveMCS   ! Mol . probably not needed but not harmful
      end
