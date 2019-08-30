      subroutine cprimFlux(
     *  comp, e_or_p, rigth, cos1, cos2, fai1, fai2,
     *  flux)
#include  "Zglobalc.h"
#include  "Zptcl.h"
#include  "Zprimary.h"

      real*8 azmmin, azmmax, rig, cosx
      logical degree, cosfactor
      common /ZpirmFlux/ azmmin, azmmax, rig, cosx,
     * degree, cosfactor

      type (component)::comp  ! input primary component
      real*8 e_or_p     ! input.  E or p as given in primary file
      real*8 rigth      ! input.  threshold rigidty below which geomagneic
                        !         effect appears.  make it 0 if
                        !         no rigidity cut.
      real*8 cos1, cos2 ! input.  cos zenith range (cos1 < cos2)
      real*8 fai1, fai2 ! input.  azimuthal angle range (fai1< fai2) deg.
      real*8 flux       ! output. average flux in the above range

      real*8 prob, flux0 
      call cprimAcceptance(
     *   comp, e_or_p, rigth, cos1, cos2, fai1, fai2, prob)
      call cprimFlux0(comp, e_or_p, flux0)

      flux = prob* flux0       

      end
