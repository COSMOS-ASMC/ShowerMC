      subroutine cPrCSampP(Eg, prob)
      implicit none
!        include "Zglobalc.h"
!        include "ZbasicCnst.h"
#include "ZbpCnst.h"
!             Basic constant; excerpt from ZbasicCnst.h
       real*8 r0  !  classical elecron radius. m
       real*8 alpha  ! fine structure const.
       real*8 m2Tomb ! m2 to mb conversion
       real*8 ar02   ! alpha * r0**2 in mb
       parameter(
     * r0 = 2.81794092d-15,
     * alpha = 1./137.0359895d0,
     * m2Tomb = 1./1.0d-31,
     * ar02 = alpha * r0**2 * m2Tomb
     * )


      real*8 Eg  ! input. Gamma energy  in GeV
      real*8 prob  !  output. prob. of Pair cre. /  r%l
      
      prob = 4* ar02 * 
     *   ( cScrC1  - cScrMain/6.d0)
     *  *mbtoPX0

      end

!
!     ************
      subroutine cPrCSampE(Eg, Ee)
!     ************
      implicit none
#include "ZbpCnst.h"
#include "Zmass.h"

      real*8 Ee  ! input. Gammaa energy in GeV
      real*8 Eg  !  output. sampled Ee : higher one of pair

      real*8  x, u, u1, u2, u3, term1, term2


      term1  = cScrMain/12.d0
      term2  = cScrC1 - cScrMain/4.d0
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

