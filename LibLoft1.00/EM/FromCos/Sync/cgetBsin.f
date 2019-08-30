      real*8  function cgetBsin(aPtcl, mag)
      implicit none
#include "Zglobalc.h"
#include "Zptcl.h"
#include "Zcoord.h"
#include "Zmagfield.h"

      type(ptcl):: aPtcl    ! input. electron or gamma 
      type(magfield):: mag  ! magnetic field.

      type(coord):: p, b, pb      ! pb is  P x B
      real*8 pbsin, pabs
      integer:: i

      call cvecProd(aPtcl%fm%p, mag, pb)
!        get pbsin
      pbsin = sqrt( dot_product(pb%r(1:3), pb%r(1:3) ) )
!
      call cpxyzp(aPtcl%fm, pabs)

      cgetBsin = pbsin/pabs
      
      end
