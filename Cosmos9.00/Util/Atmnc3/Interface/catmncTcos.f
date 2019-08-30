      subroutine catmncTcos(katmnc, code, subc, charge)
!      ptcl  code conversion from atmnc3 to cosmos 
!      (vice versa is not yet)
      integer,intent(in):: katmnc  ! atmnc3 code
      integer,intent(out):: code, subc, charge ! cosmos code
#include "Zcode.h"      
      include "../src3/include/atmnc-particle-code.inc"

      select case(katmnc)
      case (kep)
         code = kelec
         subc = antip
         charge = 1
      case (kem)   
         code = kelec
         subc = regptcl
         charge = -1
      case (kgamma)
         code = kphoton
         subc = regptcl
         charge = 0
      case (kkp)   
         code = kkaon
         subc = regptcl
         charge = 1
      case (kkm) 
         code = kkaon
         subc = antip
         charge = -1
      case (kk0l)
         code = kkaon
         subc =k0l
         charge = 0
      case (kk0s)
         code = kkaon
         subc =k0s
         charge = 0
      case (kneut)
         code = knuc
         subc = regptcl
         charge = 0
      case (kneutbar)
         code = knuc
         subc = antip
         charge = 0
      case (kpro)
         code = knuc
         subc = regptcl
         charge = 1
      case (kprobar)
         code = knuc
         subc = antip
         charge = -1
      case (kpip)
         code = kpion
         subc = regptcl
         charge = 1
      case (kpim)
         code = kpion
         subc = antip
         charge = -1
      case (kpi0)
         code = kpion
         subc = 0
         charge = 0
      case (kmup)
         code = kmuon
         subc = antip
         charge = 1
      case (kmum)
         code = kmuon
         subc = regptcl
         charge = -1
      case (knue)
         code = kneue
         subc = regptcl
         charge = 0
      case (knuebar)
         code = kneue
         subc = antip
         charge = 0
      case (knumu) 
         code = kneumu
         subc = regptcl
         charge = 0
      case (knumubar) 
         code = kneumu
         subc = antip
         charge = 0
      case (kalpha)
         code = kgnuc
         subc = 4
         charge = 2
      case default
         code = krare
         subc = 0
         charge = 0
         write(0,*) ' atmnc code=',katmnc, 'invalid'
      end select
      end subroutine  catmncTcos
