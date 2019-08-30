!           this must be modified depending on the int. model 
!      this is for QGS
      subroutine cgetXsIni
      call ciniQGS
      end subroutine cgetXsIni

      subroutine cgetXsInterface(pj, tg, xs)
#include "Zcode.h"
#include "Zptcl.h"
      type(ptcl):: pj    ! input projectile
      type(ptcl):: tg    ! input target
      real(8),intent(out):: xs ! in mb
      integer::TA
      if(tg.code == kgnuc) then
         TA = tg.subcode         
      else
         TA = 1
      endif
      call  cxsecQGS( pj, TA,   xs )
      end      subroutine cgetXsInterface


