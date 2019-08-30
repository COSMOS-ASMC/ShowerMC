!      module BPTsai
      module BPPS  ! paritial screening byTsai
      implicit none
!       for given media, these quantities are
!       calculated for each element in the media.
!       The calculation is perfomed by Zpart
!       and is used to get brems and pair cross-sections
!       These are independent of energy (and so x=Eg/Ee)
      private
      real(8),save:: Z= 0.
      real(8),save::  Z13, Z23, ame, apme, delta, cf,
     *        lnZ43,  lnZ83,  al183z
      
      public epBrem, epPair, epBPZpart

      contains
      function epBrem(Zin, Eeme, x) result(ans)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
!          This is by Tsai
!          From table B.4 of  Tsai (Rev. Mod. Phys.vol.46' 74).
!          compute differentical bremstralhng cross-section
!          at low energies where no LPM effect exists.
!        epBrem(media, Egme,x)= dsigma/dx : in  mb
!           alpha * r0**2 = 0.579 mb
!
      real(8),intent(in):: Zin  ! target atom Z
      real(8),intent(in)::Eeme ! input  Ee/me, i.e,  gamma factor 
      real(8),intent(in)::  x    ! input  Eg/Ee.    x <= 1.-me/Ee
      real(8)::ans   ! ds/dx in mb
 

      real(8):: phi1, phi2, psi1, psi2
      real(8):: gamma,  epsil

      call epBPZpart(Zin)
      delta = x/(2*(1.-x))/Eeme
      gamma = 200.* delta/ Z13
      epsil = 200.* delta/ Z23

      call epBPfc( gamma, epsil,  phi1, phi2, psi1, psi2)
      ans = ( (4.d0/3.d0 *(1.-x) + x**2) *
     * ( Z**2*(phi1-4*cf) + Z*psi1 )
     *  + 2.d0/3.d0 *(1-x) *
     * (Z**2*(phi1-phi2) + Z*(psi1-psi2)))/x

      ans= max(ans, 0.d0)
!          To mb, alpha r0**2 = 0.579 must be multiplied
      ans = ans * ar02
      end function epBrem
!     *******************************
      function epPair(Zin, Egme, x) result(ans)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
!          compute differentical pair creation cross-section
!          at low energies where no LPM effect exists.

!        epPair = dsigma/dx in mb
!
      real(8),intent(in):: Zin    ! charge of target
      real(8),intent(in):: Egme ! input  Eg/me
      real(8),intent(in):: x    ! input  Ee/Eg.   me/Eg =< x <= 1.-me/Eg
      real(8):: ans

      real*8 phi1, phi2, psi1, psi2
      real(8):: gamma, epsil
      real*8 epPairLowE

      call epBPZpart(Zin)
      delta = 1./(2*x*(1.-x))/Egme
      gamma = 200.*delta/Z13
      epsil = 200.*delta/Z23
      
         call epBPfc(gamma, epsil,  phi1, phi2, psi1, psi2)
         ans = max(  ( (4.d0/3.d0  *(x -1.)*x + 1) *
     *        (Z**2*(phi1-4*cf) + Z*psi1)
     *        - 2.d0/3.d0 *x *(1-x)* 
     *        (Z**2* (phi1-phi2) + Z*(psi1 -psi2)) ), 0.d0)
!             alpha *r0**2 epPair = 0,579 epPair
         ans = ans * ar02       !  in mb
      end function epPair

      
      subroutine epBPfc( gamma, epsil, phi1, phi2, psi1, psi2)
!
!       compute phi1- 4ln(Z)/3, phi2-4ln(Z)/3
!               psi1- 8ln(Z)/3, psi2-8ln(Z)/3
      implicit none

      real(8),intent(in):: gamma, epsil
      real(8),intent(out):: phi1, phi2, psi1, psi2

 
 
      if( Z <= 2.0 ) then 
         call epBPfc1(gamma,  phi1, phi2, psi1, psi2)
      elseif(Z <= 4.0 ) then
         call epBPfc2(  phi1, phi2, psi1, psi2)
      else
         call epBPfc3(gamma, epsil,  phi1, phi2, psi1, psi2)
      endif
      end subroutine epBPfc
!     ********************
      subroutine epBPfc1(gamma, phi1, phi2, psi1, psi2)
!     ********************
!        for H, He
      implicit none
      real(8),intent(in):: gamma
      real(8):: phi1, phi2, psi1, psi2

      real(8):: eta, c, arctanc, ln1, lnc1, lnc2
      real(8):: alpha, kSpence

      parameter (alpha = 1.d0/137.d0)

      if(Z .eq. 1) then
         eta = 1.
      else
         eta = 1.6875
      endif

      c = gamma*Z13/(400.d0*alpha*eta)
      arctanc = c* atan(1./c)
      ln1 = 4.*log(1./(2.*eta*alpha))
      lnc1 = 2*log(1.+c**2)
      lnc2 = c**2 * log(1.+c**(-2))

      phi1 = ln1 + 13.0/3.d0 - lnc1 - 13.d0/2.*arctanc +
     *      1./(6.*(1+1./c**2))
      phi2 =ln1 + 11.d0/3.d0 - lnc1 + 25.*c**2* 
     *       (1.-arctanc) - 14.*lnc2
      psi1 = ln1 + 23.d0/3d0 - lnc1 - 17.5*arctanc + 
     *  8*lnc2 - 1./(6.*(1+1./c**2))
      psi2 = ln1 + 21.d0/3d0 -lnc1 - 105.0*c**2*(1-arctanc) +
     *      50.0* lnc2 - 24.0*c**2 *
     *      ( -2*log(c)*log(1 + 1/c/c) +
     *        kSpence((1.d0+1/c/c))  - kSpence(1.d0))
      end subroutine epBPfc1
!     ********************
      subroutine epBPfc2( phi1, phi2, psi1, psi2)
!     ********************
!        for Z=3, 4
      implicit none
      real*8 phi1, phi2, psi1, psi2
      real*8 b, bp, lnaz, lnb, arctanb, arctanbp
      real*8 lnb2, lnbp, lnbp2, lnazp
      b = ame * delta
      bp = apme *delta
      lnaz =2* log(ame*Z13)
      lnazp = 2* log(apme*Z23)
      lnb = log(1. + b*b)
      lnb2 = log(1 + 1./b/b)
      lnbp = log(1. + bp*bp)
      lnbp2 = log(1 + 1./bp/bp)
      arctanb = atan(1./b)
      arctanbp = atan(1./bp)

      phi1 = 2.*(1. + lnaz) - 2* lnb - 4*b*arctanb
     *  - lnZ43

      phi2 = 2.*(2.d0/3.d0 + lnaz) - 2* lnb + 
     * 8*b*b*(1.-b*arctanb - 0.75*lnb2) - lnZ43

      psi1 = 2.*(1. + lnazp) - 2 * lnbp - 4*bp*arctanbp -
     *   lnZ83

      psi2 = 2.*(2.d0/3.d0 + lnazp) - 2*lnbp  + 
     * 8*bp*bp*(1.-bp*arctanbp - 0.75*lnbp2)
     * - lnZ83
      end subroutine epBPfc2

!     ********************
      subroutine epBPfc3(g, epsi, phi1, phi2, psi1, psi2)
!     ********************
!        for Z>=5
      implicit none
      real(8),intent(in):: g, epsi
      real(8),intent(out):: phi1, phi2, psi1, psi2

      phi1 = 20.863 - 2*log(1+(0.55846*g)**2) -
     *  4*(1 - 0.6* exp(-0.9*g) - 0.4*exp(-1.5*g))
     *  - lnZ43

      phi2 = phi1 - 2.d0/3.d0/(1+(6.5+ 6*g)*g )

      psi1 = 28.340 - 2*log(1 + (3.621*epsi)**2) -
     * 4*(1-0.7*exp(-8.0*epsi) - 0.3*exp(-29.2*epsi))
     *  - lnZ83
      psi2 = psi1 - 2.d0/3.d0/(1+(40 + 400*epsi)*epsi )

      end subroutine epBPfc3


      subroutine epBPZpart(Zin)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
!        compute Z  part of the pair and brems diff. 
!        corss-section (single Atom)

      real(8),intent(in):: Zin  !  charge of the element atom

!        y = (Ee-me)/(Eg-2me)

                 ! non screening region

      if( Z /= Zin ) then
         Z = Zin
         Z13 = Z**(1.d0/3.d0)
         Z23 =Z13**2
         call epGetTFM(Z)
         lnZ43 = log(Z)*4.d0/3.d0
         lnZ83 = 2* lnZ43
      endif

      end subroutine epBPZpart


      subroutine epGetTFM(Z)
!
!        get Thomas-Fermi-Moliere atomic factor a, a'
! 
!        From table B.4 of  Tsai (Rev. Mod. Phys.vol.46' 74)
!        a*me, and a'*me are obtained 
!
!    
      implicit none
      real(8),intent(in):: Z 

      real*8 sz, epCoulombC

      if( Z == 1.0 )  then
         ame = 122.8
         apme = 282.4
      elseif( Z == 2.0 ) then 
        ame = 90.8/Z13
        apme = 265.8/Z23
      elseif( Z == 3.0 ) then
         ame = 100.0/Z13
         apme = 418.6/Z23
      elseif( Z == 4.0 ) then
         ame = 106.0/Z13
         apme = 571.4/Z23
      else
!         for all  Z>=5
         ame = 111.7 / Z13
         apme = 724.2 / Z23
      endif
!        Coulomb correction fucntion
      sz =( Z/137.d0 ) **2
      cf = epCoulombC(sz)
      end subroutine epGetTFM

!      function epPairLowNorm(media, Z) result(ans)
!      implicit none
!#include "Zmedia.h"
!#include "Zmass.h"
!      record /epmedia/ media  ! input media
!      real(8),intent(in):: Z  ! charge
!      real(8):: ans   !  to be multiplied to epPair
!      real(8),parameter:: ynorm = 0.9d0
!
!      real(8):: xnorm, NonScEme
!
!      NonScEme = media.cnst.PairNonSc/masele
!
!      xnorm = ynorm * (1.-2/NonScEme) +
!     *    1./NonScEme      ! -->0.8489
!      ans = epPair(Z, NonScEme, xnorm)/
!     *      epPairLowE(Z, NonScEme, xnorm)
!      end      function epPairLowNorm
!!      end module BPTsai
      end module BPPS
