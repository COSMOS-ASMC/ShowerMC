      real*8 function epPairLowE(Z, Egme, x)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmass.h"
!          compute differentical pair creation cross-section
!          near threshold. (No screening region).
!          epPairLowE = dsigma/dx of B.H's original cross
!          section in mb
! Note:  B.H's paper (1934), Eq.(21) has wrong sign in the
!          last line. 
      real(8),intent(in):: Z  ! Atom's Z
      real*8 Egme ! input  Eg/me
      real*8 x    ! input  Ee/Eg.   me/Eg =< x <= 1.-me/Eg

      real*8  E0, Ep, p0, pp,  logt, epsp, eps0,  k

      k = Egme * masele

      if(x .ge. 1.d0 - masele/k) then
         epPairLowE =  0.
      elseif(x .le. masele/k) then
         epPairLowE =  0.
      else
         E0 = k*x
         p0 = sqrt(E0**2 - masele**2)
         Ep = k - E0
         pp = sqrt(Ep**2-  masele**2)
         
         logt = 2* log( ( E0*Ep + p0*pp + masele**2)/masele/k)
         epsp = 2 * log( (Ep + pp)/masele )
         eps0 = 2 * log( (E0 + p0)/masele )
!      Psi(x) dE0/k**3... = Psi(x) dx /k**2 ...
         epPairLowE = p0*pp/k**2 * ( -4.d0/3.d0
     *    - 2*E0*Ep * (p0**2 + pp**2)/ (p0*pp)**2 + 

     *    masele**2 *( eps0*Ep/p0**3 + epsp*E0/pp**3 - 
     *    epsp*eps0/p0/pp ) + 
         
     *    ( k**2/(p0*pp)**3 *((E0*Ep)**2 + (p0*pp)**2) -8.d0/3.d0 *
     *    E0*Ep/p0/pp) * logt  -

     *    masele**2 *k/2/p0/pp *
     *   ( (E0*Ep-p0**2)*eps0/p0**3 + (E0*Ep-pp**2)*epsp/pp**3 +
     *   2*k*E0*Ep/(p0*pp)**2)* logt)  

         epPairLowE =  epPairLowE * ar02 * Z**2
       endif
      end   function epPairLowE

      function epPairLowNorm(media, Z) result(ans)
      use BPPS, only: epPair
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
       type(epmedia)::  media  ! input media
      real(8),intent(in):: Z  ! single atom's charge
      real(8):: ans     ! normalzation factor to the epPairLowE
      
      real(8),parameter:: xnorm =0.5
      real(8)::  NonScEme
      real(8)::  epPairLowE

      NonScEme = media%cnst%PairNonSc/masele

      ans = epPair(Z, NonScEme, xnorm)/
     *      epPairLowE(Z, NonScEme, xnorm)
      end function epPairLowNorm


