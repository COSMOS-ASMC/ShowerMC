      subroutine epPrCSampP(media,  Eg,  prob)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"

       type(epmedia)::  media  ! input
      real*8 Eg  ! input. Gamma energy  in GeV
      real*8 prob  !  output. prob. of Pair cre. /  r.l
      
      prob = 4* ar02 * 
!       ( media.cScrC1 - media.cScrMain/4.d0 
!     *     media.cScrMain/12.d0  )
     *   ( media%cScrC1  - media%cScrMain/6.d0)
     *  * media%mbtoPX0
      end

!
!     ************
      subroutine epPrCSampE(media, Eg, Ee)
!     ************
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
       type(epmedia)::  media  ! input
      real*8 Ee  ! input. Gammaa energy in GeV
      real*8 Eg  !  output. sampled Ee : higher one of pair

      real*8  x, u, u1, u2, u3, term1, term2


      term1  = media%cScrMain/12.d0
      term2  = media%cScrC1 - media%cScrMain/4.d0
      do while (.true.)
         call rndc(u)
         if(u .lt.  term2/(term1+term2)) then
!              uniform
            call rndc(u)
         else
!              x**2
            call rndc(u1)
            call rndc(u2)
            call rndc(u3)
            u = max(u1,u2,u3)
         endif
         x = (u +1.d0)/2.d0 
         if(x  .lt. ( 1.d0 - masele/Eg)) goto 10
      enddo
 10   continue
      Ee = x * Eg
      end

