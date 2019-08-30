      subroutine cprimAcceptance(
     *    comp, e_or_p, rigth, cos1, cos2, fai1, fai2, prob)
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
      real*8 prob       ! output. average prob. that the primary can come

      type (ptcl):: aPtcl
!

      external crigfunc1
      real*8 crigfunc1
      real*8  error, ans

      integer icon
      data eps/1.d-4/



      azmmin = fai1
      azmmax = fai2
      call cconv_prim_e(comp, e_or_p, aPtcl)
!         rigidity
      rig =
     *   sqrt( aPtcl%fm%p(4)**2 - aPtcl%mass**2 )/aPtcl%charge

      if(rig .lt. rigth) then
         call kdexpIntF(crigfunc1, cos1, cos2, eps, ans, error, icon)
         prob= ans/( (fai2-fai1) * (cos2 - cos1 ))   ! fai is in deg. o.k
      else
         prob = 1.0
      endif
      end
!     ******************
      real*8 function crigfunc1(zen)
      implicit none
#include  "Zglobalc.h"
#include  "Zptcl.h"
#include  "Zprimary.h"

      real*8 azmmin, azmmax, rig, cosx
      logical degree, cosfactor
      common /ZpirmFlux/ azmmin, azmmax, rig, cosx, 
     * degree, cosfactor

      real*8 zen

      real*8 ans


      external crigfunc2
      real*8 crigfunc2

!        integrate over phi = azmmin, azmmax
      cosx = zen
      call k16pGaussLeg(crigfunc2, azmmin, azmmax, 16, ans)
      crigfunc1 = ans 
      if(cosfactor) then
         crigfunc1 = crigfunc1 * cosx
      endif
      end
!  ******************
      real*8 function crigfunc2(phi)
      implicit none
#include "Zglobalc.h"
#include "Zptcl.h"
#include "Zprimary.h"
      real*8 phi   ! iput in degree
!         for fixed cos and phi; get  R(theta, fai, rig)

      real*8 azmmin, azmmax, rig, cosx
      logical degree, cosfactor
      common /ZpirmFlux/ azmmin, azmmax, rig, cosx, 
     * degree, cosfactor

      real*8 temp
!
      real*8 prob
!
      if(degree) then
         temp = acos(cosx)/Torad
      else
         temp = cosx
      endif
      call crigCut(phi,  temp,  rig, prob) 
      crigfunc2 = prob
      end
!     ***************************
      subroutine csetCosdeg(cosin,  degin)
      implicit none
      logical cosin, degin
      real*8 azmmin, azmmax, rig, cosx
      logical degree, cosfactor
      common /ZpirmFlux/ azmmin, azmmax, rig, cosx, 
     * degree, cosfactor
      
      cosfactor = cosin
      degree = degin
      end
