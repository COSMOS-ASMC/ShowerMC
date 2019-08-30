      real*8 function epPrgenex(x)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

!
!         epPrgenex. in mb
!  This is the version for pair creation corrsponding to
!  epBrgenex
!     
      real*8 x  !  input Ee/Eg. 


      real*8 epPrgene
      real*8 Eg

      Eg= Egme * masele
      epPrgenex = epPrgene(media, Eg, x)

      end function epPrgenex

      subroutine epPrgenePreInte(mediain, Egmein)
!    this must be called before  epBrgenex
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"
       type(epmedia)::  mediain  ! input media
      real(8),intent(in):: Egmein  ! Ee/me

      media = mediain
      Egme = Egmein
      end
