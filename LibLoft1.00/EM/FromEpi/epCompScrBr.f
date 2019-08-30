      function epCompScrBrs(media, x) result(ans) 
      implicit none
#include "Zmedia.h"
       type(epmedia)::  media ! input media
!      real(8),intent(in):: Eeme ! Ee/me not needed (no E dep)
      real(8),intent(in):: x    ! Eg/Ee
      real(8):: ans   !  ds/dx in mb for the media

      real(8)::  epCompScrBr
      integer:: i
      ans = 0.
      do i = 1, media%noOfElem
         ans = ans  +
     *        epCompScrBr(media%elem(i)%Z, x)*media%No(i)
      enddo
      end    function epCompScrBrs

      function epCompScrBr(z,  y) result(ans)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"

!         give complete screening brems cross-section
!        by Tsai's formula. (R.M.P,  Eq.3.81)
!          
      real*8 z  ! input.  atomic number of the media (single atom)
      real*8 y  ! input.  Eg/Ee
      real(8)::ans ! ds/dx for the give charge Z atom
!     function value ! output. dsigma/dy in mb
!
      real(8),save:: Lrad,  Lradp
      real(8),save::  zsave=0
      real(8),save::  z2, f

      if(zsave .ne. z) then
         zsave = z
         call epGetLrad(z, Lrad, Lradp, f)
         z2 = z*z
      endif
!
      ans = 4.* ar02 * (
     *          (4./3.- 4./3.*y + y*y) * (
     *         z2 * (Lrad-f) + z*Lradp  )
     *       + (1-y)*(z2 + z)/9.)/y
      end   function epCompScrBr
!     ***************************************
      subroutine epGetLrad(z, Lrad, Lradp, f)
      implicit none          

      real*8 z  ! input.  atomic number of the media
!
      real*8 Lrad  ! output. See Tsai's formula. (R.M.P,  Eq.3.81)
      real*8 Lradp ! ouptut.
      real*8 f     ! output.  Coulomb correction func.

      real*8 z2, z3
      real*8 epCoulombC


      f = epCoulombC((z/137.)**2)         
      z2 = z*z
      if(z .eq. 1.0) then
         Lrad = 5.31
         Lradp = 6.144
      elseif(z .eq. 2.0) then
         Lrad = 4.70
         Lradp = 5.621
      elseif(z .eq. 3.0) then
         Lrad = 4.74
         Lradp = 5.8505
      elseif(z .eq. 4.) then
         Lrad = 4.71
         Lradp = 5.924
      else
         z3=z**(-0.3333333)
         Lrad = log(184.15*z3)
         Lradp = log(1194.0*z3**2)
      endif
      end

