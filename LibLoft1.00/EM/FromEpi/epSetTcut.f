      subroutine  epSetTcut(media, Recoilmin) ! 9.154 separated from
          ! epSetEmin.f
      implicit none
#include "Zmedia.h"      
      real(8),intent(in):: Recoilmin  ! recoil KE min GeV
      !      for histrical reason, media.sh.tcut and sh.w0 are
!     used for dEdx calculations;
!        in Epics     Det.cmp(cn).recoilE  could be directly used.
      type(epmedia),intent(out)::  media  !  media to be used 
                    !  media.sh.tcut and w0 are set
      media%sh%tcut = Recoilmin  ! GeV
      media%sh%w0  = Recoilmin*1000.d0 ! MeV 
      call epResetUrban( media, media%urb)
      end

