      real*8 function epmuvmax( E )
      implicit none
#include "Zmuint.h"
#include "Zmass.h"
!         compute max fractional energy of muon brems gamma, or
!         virtual gamma at pair production by muon
!
!
      real*8 E  ! input.  muon total energy

      epmuvmax = 1.d0 - 3.d0/4.d0 *sqrt(exp(1.d0))*(masmu/E) *
     *           Z3
      end
