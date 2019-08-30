!       ********************************************8
!	The total energy loss rate by synchrotron radiation
!       (GeV/m)
!       
	real*8 function csyncTELoss(u)
	implicit none
#include "Zglobalc.h"
	real*8 u  ! input. upsilon
        real*8  cgUpsilon
	csyncTELoss = 0.66666 * SyncConvR * cgUpsilon(u)
	end
