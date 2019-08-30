      subroutine epqEmin(cmpn, minErg)
!        get   Emin for the given component number
!       This can be used after config file and epicsfile
!       have been read.  

      implicit none
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zmass.h"
#include "Zcnfig.h"
      integer,intent(in):: cmpn  ! component #
      real(8),intent(out):: minErg(5) !  kinetic energy in GeV
         ! (1): for gamma (2): electron (3): recoil min.
         ! (4): other charge  (5): neutron
      if(cmpn /= Cn ) then
         call epSetEmin(cmpn)
      endif
      minErg(1) = EminGamma
      minErg(2) = EminElec - masele
      minErg(3) = RecoilKEmin
      minErg(4) = KEmin
      minErg(5) = EminH
      if(cmpn /= Cn ) then
!              reset Emin for safety
         call epSetEmin(Cn)
      endif
      end
