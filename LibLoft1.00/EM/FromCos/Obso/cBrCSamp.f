      subroutine cBrCSampP( Ee, prob)
      implicit none
#include "Zglobalc.h"
#include "Zmass.h"
#include "ZbpCnst.h"
!             Basic constant; excerpt from Zbasiccnst.h
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


 


      real*8 Ee  ! input. Electron energy in GeV
      real*8 prob  !  output. prob. of Brems /  r%l
      
      real*8 vc

      vc = BremEgmin/(Ee-masele)
      prob = 4. * ar02 * ( ( log(1.d0/vc) -
     *   (1.d0-vc)) * cScrMain +
     +   (1.d0-vc)* (1.d0+vc)/2 * cScrC1 )
     *   * mbtoPX0
      end

!
!     ************
      subroutine cBrCSampE(Ee, Eg)
!     ************
      implicit none
#include "Zmass.h"
#include "ZbpCnst.h"

      real*8 Ee  ! input. Electron energy in GeV
      real*8 Eg  !  output. sampled Eg

      real*8 vc, x, u, u1, u2, term1, term2

      vc = BremEgmin/(Ee-masele)


      term1  = cScrMain * (log(1.d0/vc) - (1.d0-vc))
      term2  = cScrC1 * (1.0d0-vc)*(1.0+vc)/2.d0
      call rndc(u)
      if(u .le.  term1/(term1+term2)) then
!         (1/x -1)dx
         do while (.true.)
!              average number of trials is 1.0xx; xx depends on vc
            call rndc(u)
            x = vc**u
            call rndc(u)
            if( u .lt. (1.0-x)) goto 10
         enddo
 10      continue
      else
!          x dx in [vc:1]
         do while (.true.)
            call rndc(u1)
            call rndc(u2)
            x = max(u1,u2)
            if(x .gt. vc ) goto 20
         enddo
 20      continue
      endif
      Eg = (Ee-masele) * x
      end
