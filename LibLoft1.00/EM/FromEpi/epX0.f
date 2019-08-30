!     ****************************************************************
!     *
!     * epX0: compute radiation length of a mixed matter
!     *
!     ****************************************************************
!
       subroutine epX0(media)
       implicit none
#include "Zmedia.h"

       type(epmedia)::  media
!
       real*8  X0, X0i,  zr
       integer i

       X0i=0.
       do   i=1, media%noOfElem

          zr = media%elem(i)%Z  ! to make real*8
          call epX01(zr, media%elem(i)%A, X0)
          X0i = X0i + media%w(i)/X0

       enddo

       media%X0g = 1./X0i
       media%X0 = media%X0g/media%rho
      end
!     ****************************************************************
!     *                                                              *
!     * epX01:  compute radiation length of a simple matter         *
!     *                                                              *
!     *********************** tested 80.07.11 ************************
!
!   /usage/
!            call epX0l(z, a, x0g)
!
!    z:  charge of the matter
!    a:  mass no.
!  x0g:  //                           g/cm**2
!
!
!
!
      subroutine epX01(z, a, x0g)
      implicit none
!         cnst=
!         4/137* r0**2 * n  where r0 is the classical electron radius
!                           n the avogadro number
!                           r0=2.8176e-13 cm
!                           n=6.0247
!
      real*8  z, a, x0g
      

      real*8  z3, az, Lrad, Lradp, fz
      real*8 t0inv
      real*8  cnst/1.396e-3/
      real*8  epCoulombC 
!     
      az = (z/137.0)**2
!         Coulomb correction func.

      fz = epCoulombC(az)

      if(z .eq. 1.0) then
         Lrad = 5.31
         Lradp = 6.144
      elseif(z .eq. 2.0) then
         Lrad = 4.70
         Lradp = 5.621
      elseif(z .eq. 3.0) then
         Lrad = 4.74
         Lradp = 5.8505
      elseif(z .eq. 4.) then
         Lrad = 4.71
         Lradp = 5.924
      else
         z3=z**(-0.3333333)
         Lrad = log(184.15*z3)
         Lradp = log(1194.0*z3**2)
      endif
!        inverse of r.l in g/cm**2
      t0inv=cnst / a  * ( z*z *(Lrad - fz ) + z*Lradp)
!
      x0g=1. / t0inv
      end
