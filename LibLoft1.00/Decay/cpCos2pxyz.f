      subroutine cpCos2pxyz(cosa, p, pxyz)
      implicit none
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
      real*8 cosa, p

      type(fmom):: pxyz
      real*8  cs, sn, tmp, sina

      call kcossn(cs, sn)
      tmp=1.- cosa**2
      sina=sqrt(tmp)
      pxyz%p(1) = sina*cs*p
      pxyz%p(2) = sina*sn*p
      pxyz%p(3) = cosa*p
      end
