      subroutine cprimFlux0(comp, e_or_p, flux)
      implicit none
#include "Zptcl.h"
#include "Zprimary.h"
      type (component):: comp  ! input. primary component
      real*8 e_or_p   ! input.  Energy or momemntum as given in primary file
      real*8 flux     ! output  flux vlaue at e_or_p
                      !      if e_or_p is outside of the table, flux=0


      integer j
      real*8 beta

      call kdwhereis(e_or_p, comp%no_of_seg+1, comp%energy, 1, j)
      if(j .gt. comp%no_of_seg) then
         flux = 0.
      elseif(j .lt. 1) then
         flux  = 0
      else
         beta = comp%beta(j)
         flux = comp%flux(j)*(e_or_p/comp%energy(j))**(-beta)
      endif
      end
