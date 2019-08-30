      function epCompScrPrs(media,  x) result(ans)
      implicit none
#include "Zmedia.h"
       type(epmedia)::  media ! input media
!      real(8),intent(in):: Egme ! Eg/me not needed; E indep.
      real(8),intent(in):: x    ! Ee/Eg
      real(8):: ans   !  ds/dx in mb for the media

      real(8)::epCompScrPr
      integer:: i
      ans = 0.
      do i = 1, media%noOfElem
         ans = ans  +
     *        epCompScrPr(media%elem(i)%Z, x)*media%No(i)
      enddo
!
      end         function epCompScrPrs



      function epCompScrPr(z, y) result(ans)
!         complete screening cross-section by Tsai for pair creation.
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
      
      real*8 z ! input. atomic number of the media. single element
      real*8 y ! input. Ee/Eg
      real(8)::ans  ! ds/dy in mb

      real(8),save:: Lrad, Lradp, z2, f
      real(8),save:: zsave=0.

      if(zsave /= z) then
         zsave = z
         call epGetLrad(z, Lrad, Lradp, f)
         z2 = z*z
      endif

      ans = 4* ar02 * (
     *          (4./3.*y*y- 4./3.*y + 1) * (
     *            z2 * (Lrad-f) +  z*Lradp )
     *       -   y*(1-y)*(z2 + z)/9.)

      end   function epCompScrPr


          

