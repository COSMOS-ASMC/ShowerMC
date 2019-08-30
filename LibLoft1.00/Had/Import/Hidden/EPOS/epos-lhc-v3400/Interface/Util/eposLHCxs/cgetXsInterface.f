!           this must be modified depending on the int. model 
!      this is for EPOS
      subroutine cgetXsIni
      call ceposIniAll
      end subroutine cgetXsIni

      subroutine cgetXsInterface(pj, tg, xs)
#include "Zptcl.h"
      type(ptcl):: pj    ! input projectile
      type(ptcl):: tg    ! input target
      real(8),intent(out):: xs ! in mb
      call ceposIniOneEvent(pj, tg, xs)
      end      subroutine cgetXsInterface



