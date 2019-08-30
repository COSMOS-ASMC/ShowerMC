!  old cgetBsin
!     MacGfort
!       bsin =  2.3737987772599083
!
!      real0m0.154s
!      user0m0.144s
!      
!     MacIFC
!       bsin  =  2.37379877725991
!
!      real0m0.178s
!      user0m0.156s
!     sys0m0.011s
!      
!  New cgetBsin
!!     MacGfort
!     2.3737987772599083
!
!      real0m0.145s
!      user0m0.137s
!      sys0m0.005s
!      
!      
!     MacIfC
!    2.37379877725991
!
!      real0m0.073s
!      user0m0.063s
!      sys0m0.006s
      implicit none
#include "Zmagfield.h"
#include "Zptcl.h"      
      type(magfield):: Mag
      type(ptcl):: aptcl
      real(8):: bsin
      
      integer::i
      
      Mag%x = 1.
      Mag%y = -1.
      Mag%z = 3.0

      aptcl%fm%p(1) = 10.
      aptcl%fm%p(2) = -1.
      aptcl%fm%p(3) =  5.


      
      do i = 1, 10000000
         call testspeed(aptcl, Mag, bsin)         
      enddo

      
!      write(0,*) bsin

      write(0,*)  bsin
      end

      subroutine testspeed( aptcl, Mag, bsin)
#include "Zmagfield.h"
#include "Zptcl.h"      
      real(8),external:: cgetBsin
      real(8),external:: cgetBsinx
      real(8):: bsin

!      bsin= cgetBsin(aptcl, Mag)
      bsin= cgetBsinx(aptcl, Mag)      
      end
      
      real*8  function cgetBsinx(aPtcl, mag)
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
!           for safety
!n      do i = 1, 3
!n	      p%r(i) = aPtcl%fm%p(i)
!      enddo

      call cvecProd(aPtcl%fm%p, mag, pb)
!        get pbsin
      pbsin = sqrt( dot_product(pb%r(1:3), pb%r(1:3) ) )
!
      call cpxyzp(aPtcl%fm, pabs)

      cgetBsinx = pbsin/pabs
      
      end
      
      
