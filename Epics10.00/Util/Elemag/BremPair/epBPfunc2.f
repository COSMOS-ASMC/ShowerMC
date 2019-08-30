!      module BPNelson
      module BPPS   ! partial screening by Nelson
!            brems function is  not available for Nelson case
      implicit none
!       For given media, these quantities are
!       calculated for each element in the media.
!       The calculation is perfomed by Zpart
!       and is used to get brems and pair cross-sections
!       These are independent of energy (and so x=Eg/Ee)
      private
      real(8),save:: Z=0
      real(8),save::  Z13, Z23, ame, apme,  cf,
     *        lnZ43,  lnZ83,  l183z
      public epPair, epBrem,  epBPZpart   ! epBrem is dummy

      contains
      function epPair(Zin, Egme,  x) result(ans)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
!         This is by Nelson et al. No brems func included

!          compute differentical pair creation cross-section
!          at low energies where no LPM effect exists.
!        Before calling this, epBPZpart must have been called
!        to give the Z
!        epPair =  dsigma/dx   in mb
!
      real(8),intent(in):: Zin
      real(8),intent(in):: Egme ! input  Eg/me
      real(8),intent(in):: x    ! input  Ee/Eg.   me/Eg =< x <= 1.-me/Eg
      real(8):: ans   ! ds/dx  in mb

      real*8 phi1, phi2, fz, gzai
      real*8 epPairLowE
      real(8):: delta

      call epBPZpart(Zin)
      delta = 136./(x*(1.-x))/Egme/Z13
      if( delta .lt. 1.) then
         phi1 = (0.625 * delta - 3.242)*delta + 20.867
         phi2 = (-0.086 * delta -1.930)*delta + 20.209
      else
         phi1 = 21.12 - 4.184*log(delta + 0.952)
         phi2 = phi1
      endif
      
      if(Egme .lt. 25.0) then         
!         if(Egme .lt. 97.8) then         ! nelson original
!              Eg < 0.05 GeV
         fz = lnZ83
!///////////
!         else
!            fz = lnZ83 + 8*cf
!         endif
      elseif( Egme .lt. 150) then
         fz =
     *       lnZ83 + 8*cf* ((Egme-25.0)/(150.0-25.0))**0.5
      else
         fz = lnZ83 + 8*cf
      endif
!///////////////
         
      gzai = log(1440./Z23)/ (log(183./Z13)-cf)

      ans = Z*(Z+gzai)*
     *        (  (x**2 + (1-x)**2)* (phi1 - fz/2)  +
     *        2.d0/3.d0 * x*(1-x)*(phi2 -fz/2)
     *        )
!            alpha *r0**2 epPair
      ans =max( ar02 * ans, 0.d0) ! in mb
      end function epPair
!     ***********************************
      function epBrem(Zin, Eeme, x) result(ans)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
      real*8 Zin
      real*8 Eeme   
      real*8 x
      real(8):: ans
!        this is dummy;
      ans = 1.
      call cerrorMsg(
     * 'epBrem not avail in Nelson  formula', 1)
      
      end function epBrem
!     ***************************
      subroutine epBPZpart(Zin)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"

!        compute Z  part of the pair and brems diff. 
!        corss-section (single Atom)

      real(8),intent(in):: Zin

!
!       The computed results are put in the moudle
! 

      

      real*8 sz, epCoulombC

      if( Z /= Zin) then
         Z =  Zin
         Z13 = Z**(1.d0/3.d0)
         Z23 = Z13**2
!        Coulomb correction fucntion
         sz =( Z/137.d0 ) **2
         cf = epCoulombC(sz)

         lnZ43 = log(Z)*4.d0/3.d0
         lnZ83 = 2* lnZ43
      endif
      end subroutine epBPZpart

      end module BPPS

