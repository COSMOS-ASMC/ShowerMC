!        cupsilon:  compute upsilon
!        cy2zeta:  y --> zeta
!        cx2zeta:  x --> zeta
!
!            compute critical value Upsilon
!           = E/m *  H/Hcr
      real*8  function cupsilon(electron, mag)
      implicit none
#include "Zglobalc.h"
#include "Zptcl.h"
#include "Zmagfield.h"

      type(ptcl):: electron ! electron
      type(magfield):: mag  ! magnetic field.
      real*8 bsin, cgetBsin

      bsin = cgetBsin(electron, mag)
!            E/m  * Bsin/Bcr      
      cupsilon = electron%fm%p(4) /electron%mass * bsin/Bcr
      end
!     *****************************************************************
!      compute critical energy (peak energy) of photons for synchrotron
!      radiation.
!        
      real*8 function cegCrit(e, upsilon)
      implicit none
!         y in Erber is Eg/Egcrit

      real*8  e  ! input.  Electron energy in GeV.
      real*8 upsilon   !  Upsilon value  E/m * Bsin/Bcr
!
      cegCrit = e * 3. * upsilon/(2. + 3.*upsilon)
      end
!     *****************************************
!        zeta = y/(2 +  3Upsilon(1-y))
!             = y/(2+3Upsilon)/(1- x)
!             = x (1+x )/ (3Upsilon)          Eq.2.8
      real*8 function cy2zeta(y, upsilon)
      implicit none
      real*8  y  ! input.  Eg/Egcrit
      real*8 upsilon ! input.  

      cy2zeta = y/ (2. + 3.*upsilon*(1. - y))
      end
!     ********************************************
!        This is approx. formula.
!
      real*8 function cx2zeta(x, upsilon)
      implicit none
      real*8  x  ! input.  Eg/Ee
      real*8 upsilon ! input.  

      cx2zeta = x *(1. + x)/(3.*upsilon)
      end
