!           compute |p| from px, py, pz
       subroutine cpxyzp(po, pabs)
       implicit none
#include  "Zptcl.h"
       type(fmom):: po
       real*8 pabs
       pabs=sqrt(po%p(1)**2+po%p(2)**2+po%p(3)**2)
       end
