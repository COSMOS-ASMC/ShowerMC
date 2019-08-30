!           get beta for Time calculation
!
      subroutine  cgetBeta(aPtcl, beta)
      implicit none
!----      include '../Particle/Zptcl.h'
#include  "Zptcl.h"
      type(ptcl)::aPtcl
      real*8 beta
!
      real*8 g, gmin/100./
!
      if(aPtcl%mass .eq. 0.) then
         beta = 1.
      else
         g =aPtcl%fm%p(4)/aPtcl%mass 
         if(g .lt. gmin) then
            beta = sqrt(1. - 1./g/g)
         else
            beta =(-1./g/g/8.0-0.5)/g/g + 1.
         endif
      endif
      end
