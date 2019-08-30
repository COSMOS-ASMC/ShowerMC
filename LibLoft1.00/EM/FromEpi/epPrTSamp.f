      subroutine epPrTSampP(media, Eg, prob)
!
!      ********* This is actually not used now.***********
!
!      sampling of occurence of pair electron production / X0
!      near threshold energy.
!      total cross-section in unit of Z^2 ar02 was computed
!      by using testPairLowE2.f and approximated by 
!      polinomials in two devided regions.
!
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "Zglobalc.h"
#include "ZbasicCnst.h"

       type(epmedia):: media
      real*8 Eg
      real*8 prob  ! output probability of Pair / X0
      real*8 thresh, y, s
      parameter ( thresh = 2.001d0 * masele)
      integer i
      real*8 coef1(5), coef2(5)
!               y< 3; y = (Eg-2me)/me
      data ( coef1(i), i=  1,   5)/
     1    0.77366610E-03, -0.31574438E-01,  0.14388056,
     2   -0.26702841E-01,  0.16673182E-02/
!            3  < y  40
      data ( coef2(i), i=  1,   5)/
     1  -0.40412775,  0.37073868, -0.11235537E-01,
     2   0.21697910E-03, -0.17697363E-05/
!         
 
 
      if(Eg .lt. thresh) then
         prob= 1.d-40
      else
         y = (Eg - thresh)/masele
         if(y .lt. 3.0) then
!             better than 0.1 %
            s  =  0.
            do i = 5,  1, -1
               s =  s * y + coef1(i)
            enddo
         else
            s = 0.
!             better than 0.1 %
            do i = 5, 1, -1
               s =  s * y + coef2(i)
            enddo
         endif
!           media.Z2 is  sum of No*Z^2 of each element
         prob = s * media%Z2 * ar02 * media%mbtoPX0
      endif
      end
!     ************
      subroutine epPrTSampE(Eg, Ee)
!     ************
!          samples higher energy pair electron; 
!       We know that the pair energy spectrum is well approximated
!     by the following function where x is (Ee-me)/(Eg-2me)
!      ( 0<x<1).
!      f(x) = sqrt(x(1-x)) (x**p + (1-x)**p)
!     p is the function of Eg as
!          p = 0.022(Eg/me) + 0.956
!     This approximation is good up to Eg~50me where B.H non screened
!     cross-section can be applied. 
!
      implicit none
#include "Zmass.h"

      real*8 Ee,  Eg

      real*8 u, p, x
      if( Eg > 5e-3 .and. Eg < 10.e-3) then
         p = 0.0291*Eg/masele + 0.911249
      else
         p = 0.022d0* Eg/masele + 0.956d0
      endif
!        the sampling function is x^(1/2+p)(1-x)^(1/2)  dx
!                               + x^(1/2) (1-x)^(p+1/2) dx
!   
      call rndc(u)
      if(u .lt. 0.5) then
!         use  first term; note power is diff. by 1.         
         call kbetar( p+1.5d0, 1.5d0, x)
      else
!         use  second term
         call kbetar( 1.5d0,  p+1.5d0, x)
      endif
      if(x .lt. 0.5) then
         x = 1.-x
      endif
      Ee = x * (Eg-2*masele)  + masele
      end
