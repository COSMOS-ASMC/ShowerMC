      real*8 function epmuvmax2(media,  E)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
!         compute max fractional energy of muon brems gamma, or
!         virtual gamma at pair production by muon
!      This is the same one as epmuvmax (in Util/Mu; now not) 
!      except for
!      Z3 is  thru media and averaged one if  not sigle atom
!
       type(epmedia):: media  !input.
      real*8 E  ! input.  muon total energy

      epmuvmax2 = 1.d0 - 3.d0/4.d0 *sqrt(exp(1.d0))*(masmu/E) *
     *            media%Z1_3rd
!     *            media.Z**0.3333
!     *            media.Zeff3

      end
