      subroutine cepMagPairL(media, pj, Bfield, rrho)
!        this is to be used when Bfield and rrho are fixed
      implicit none
#include "Zmedia.h"
#incldue "Zptcl.h"
      type(epmedia),intent(in):: media
      type(ptcl),intent(in):: pj
      real(8),intent(in):: Bfield(3) ! \vec B in T  in the current
                                     ! coordinate system.
            
      real(8),intent(in):: rrho !  rho(actual)/media%rho.
                         ! in EPICS this is media%rhoc
!      if(MagPair .eq. 1) then    ! moved to outside
      call epmpairp(pj, Bfield, Xai, pairmfp, dl)
      dx = dl / media%gtocm * rrho  ! once convert to kg/m2
      call csetIntInf( dx, .false., 'mpair')
!      endif
      end
