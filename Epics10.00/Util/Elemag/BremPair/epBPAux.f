!     epPairS and epBremS are interface to epPair and epBrem
!     defeined in epBPfunc1,2,3.  (one of 1,2, or 3 must be selected
!     and used).  
!     Two function here are called from epBrgeneric
!
!    *************************
      function epPairS(media, Egme, x) result(ans)
!    *************************
      use BPPS, only:epPair
      implicit none
#include "Zmedia.h"
#include "Zmass.h"

       type(epmedia)::  media  ! input
      real(8),intent(in):: Egme  ! Eg/me
      real(8),intent(in)::  x    !  Eg/Ee
      real(8):: ans  !    ds/dx in mb
      
      real(8),save:: zsave=0.

      real(8):: epPairLowE, epPairLowNorm 
      integer:: i
      real(8):: temp

      ans = 0.
      do i = 1, media%noOfElem
         if( Egme >= media%cnst%PairNonSc/masele ) then
            temp=  epPair(media%elem(i)%Z, Egme, x)
         else
            temp = epPairLowE(media%elem(i)%Z, Egme, x) *
     *           epPairLowNorm(media, media%elem(i)%Z)
         endif
         ans = ans + temp* media%No(i)
      enddo
      end

!    *************************
      function epBremS(media, Eeme, x) result(ans)
!    *************************
      use BPPS, only:epBrem
      implicit none
#include "Zmedia.h"
       type(epmedia)::  media ! input
      real(8),intent(in):: Eeme ! Ee/me
      real(8),intent(in)::  x   ! Eg/Ee
      real(8):: ans  !  ds/dx in mb for this media


      integer i


      ans = 0.
      
      do i = 1,  media%noOfElem
         ans = ans + epBrem(media%elem(i)%Z, Eeme, x)*
     *               media%No(i)
      enddo
      end
