!     **********************************************************
	real*8 function cmBremMFP(erg, upsilon, xc)
	implicit none

#include "Zglobalc.h"

	real*8 erg  ! input electron energy, in GeV.
        real*8 upsilon  ! input. Upsilon value
        real*8 xc   !   input.  fix the cut off, below which the synchroron 
                    !           is treated as energy loss only.
                    ! one of   3.1627 x 10-3, 10-3, ... 10-6 or 0  
	            ! if xc is not one of these, neareset one
                    ! is  chosen.
!              xc should be 0. in v3.0 or later.
!     
!        The 'cross-section' for emitting fractional gamma ray energy x ~ x+dx,
!        in the unit distance  is given by 
!        P(Ee, x, U)dx = root(3)/2/pi (SyncConvR)/Ee* U *
!       [  (1-x)/x K1(2zeta)dx + x K2(2zeta)]dx
!        (/meter, if Ee is in GeV). ( U is   E/m * B/Bc ).
!       where, K1(z) = z Int(z,inf)K5/3(z)dz; known as Brems func.
!       and    K2(z) = zK2/3(z).
!
!       The integral value of the first term and the second term in []
!       above is given cmBremI1 and cmBremI2, resp.
!
!       
        real*8 cmBremI1, cmBremI2
	real*8 const, fz
	parameter ( const =  1.732050808/2/pi * SyncConvR )

	fz = cmBremI1(upsilon, xc) + cmBremI2(upsilon)
!              mfp in meter.
	cmBremMFP = 1.d0/
     * (const / erg * upsilon * fz)
	end

