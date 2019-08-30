#include "BlockData/cblkGene.h"  
!           this must be modified depending on the int. model 
!      this is for EPOS
      subroutine cgetXsIni
      use modsibyllXs
      implicit none
      call csibyllinit
      end subroutine cgetXsIni

      subroutine cgetXsInterface(pj, tg, xs)
      use modsibyllXs
      implicit none
#include "Zptcl.h"
      type(ptcl):: pj    ! input projectile
      type(ptcl):: tg    ! input target
      real(8),intent(out):: xs ! in mb
      call csibyllXs(pj, tg, xs)
      end      subroutine cgetXsInterface



