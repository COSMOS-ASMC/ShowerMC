      real*8 function cmuvmax2(E)
      implicit none
#include "Zcmuint.h"
#include "Zmass.h"
!         compute max fractional energy of muon brems gamma, or
!         virtual gamma at pair production by muon
!      This is the same one as cmuvmax in Util/Mu except for
!      Z3 is  thru media and averaged one if  not sigle atom
!
      real*8 E  ! input.  muon total energy

      cmuvmax2 = 1.d0 - 3.d0/4.d0 *sqrt(exp(1.d0))*(masmu/E) *
     *           Zeff3

      end
