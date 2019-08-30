      subroutine cmoliere(rhoin, 
     *           z, mass, g1, g2,  leng, teta, cond)
      implicit none
!        Moliere theory of multiple scattering angle.
!
      real*8 rhoin  ! input. average media density in kg/m^3

      integer z  ! input.  charge of the particle is ze
      real*8 mass ! input.  mass of the particle in GeV
      real*8 g1  ! input.  gamma factor at the path head
      real*8 g2  ! input.  gamma factor at the path end
      real*8 leng ! input. length the charged particle travelled in m.      
      real*8 teta ! output. sampled spatial angle in radain.
      integer cond ! output. 0 ok. non-0. Moliere theory not applicable

      real*8  rho, gbeta2, beta2, massratio2
      common /Zcmedia/ rho, gbeta2, beta2,
     *       massratio2

!     *********************
      real*8 xc2, xa2, bp, b, u
      real*8 a0, a1, a2,  sum, ra2, ra2inv

      integer icon
      real*8 rejf1, rejf21, rejf22
      real*8 x
!       
!      rejection function for redueced angle < 1.8
!          x is suqre of reduced angle better than 0.2 %
       rejf1(x)= ((0.1217176d-01*x + 0.3054916d-01)*x -0.2524543d0)*x
     *          + 0.9990352d0          
!
!              at 0 < x = 1/angle^2 < 0.15
       rejf21(x) =(( -162.1568*x + 44.48334)*x + 0.3907116)*x 
     *      + 0.4399338              

!             at  x = 1/angle^2 > 0.15
       rejf22(x) = (( 71.23819*x - 49.61906)*x + 10.77501)*x+ 0.2001145      
!
! ------------------------------------
      rho = rhoin
!        bgeta2   (g-1/g) = g*beet^2 = g*(1-1/g^2)
      gbeta2=(g1 - 1.d0/g1) * (g2 - 1.d0/g2)
      beta2 = 1.d0 - 1.d0/g1/g2
      massratio2= (0.511d-3/mass)**2  ! (me/m)^2
!    ............................

!          get Xc^2
      call ckaic2(z,  leng,  xc2)
!          get Xa^2
      call ckaia2(z,  xa2)

!          b -log(b) = b'
      bp = log(xc2/xa2/1.167)

      if(bp .lt. 3.395) then
!         Moliere theory cannot be appliled; use Gaussian later (almost no scattering)
         cond = 1
      else
         cond = 0
         call cblogb(bp, b, icon)
         a0 = max(1.d0 - 5/b, 0.d0)  ! use single scattering term if b<=5. 
         icon = 1               ! make 0 if no rejection
!                the sampling function decomposition is explained in Test/....tex
         do while (icon .ne. 0)
            a1 = 5.21062/b
            a2 = 0.7128/b
            sum = a0 + a1 + a2
            call rndc(u)
            if(a0/sum  .gt. u) then
!             sample reduced angle from exp(-x) dx where x = reduced
!               angle^2.
               call rndc(u)
               ra2 = -log(u)
               icon =0
            elseif( (a0+a1)/sum .gt. u) then
!            sample reduced angle from exp(-x) dx (same as above but
!                  in the region of ra < 1.8
               call rndc(u)
               ra2 = -log(1.-u/1.04076)
!                rejection function
               call rndc(u)
               if(u .lt. rejf1(ra2)) then
                  icon = 0
               endif
            else
!             sample reduced angle from 2xc2 x^-4dx
               call rndc(u)
               ra2 = 3.24/u
!               rejection function
               call rndc(u)
               ra2inv = 1./ra2
               if(ra2inv .lt. 0.15) then
                  if(u .lt. rejf21(ra2inv)) then
                     icon = 0
                  endif
               elseif(u .lt. rejf22(ra2inv)) then
                  icon = 0
               endif
            endif
         enddo
         teta =sqrt( ra2 * xc2 * b)
      endif
      end
!     ***************************
      subroutine cblogb(c, b, cond)
      implicit none
!        solve  B - log(B) = c
!
      real*8 c ! input.   c>=1.
      real*8 b ! output.  solved b >=1. (b <1 is discarded)
      integer cond ! output. 0 if ok.
                   !         1 if c < 1.
      if(c .lt .1) then
         cond = 1
      else
!         b = 0.7 + 1.32 *c -  0.01451* c*c
          b =(((-0.3733404E-04*c + 0.1902303E-02)*c -0.3841290E-01 )*c
     *         + 1.431932)*c +     0.5200487 
      endif
      end
      subroutine ckaia2(z,   xa2)
      implicit none
#include "Zair.h"
!        compute Xa^2; assume the Xa^2 is weakly
!       dependent on Z, we use average Z=zave for
!       calculation.
!
      integer z  ! input. charge of the charged particle is ze
      real*8 xa2   !  output. Xa^2.

      real*8  rho,  gbeta2, beta2, massratio2
      common /Zcmedia/ rho, gbeta2, beta2, 
     *       massratio2
      
      real*8   alpha, const, pi, large
      parameter (alpha = 1./137., const = (1.13*alpha)**2 )
      parameter (pi = 3.1415, large = (pi/2.)**2)

      if(gbeta2 .le. 0.) then
         xa2 = large
      else
!              since air Z is close N and O's Z, we
!           use simply average of Z here.
         xa2 = const * TargetZ2_3rd * massratio2 *
!     *   (1.13 * beta2 + 3.76*(alpha*z*zave)**2) /gbeta2
     *   (1.13 *beta2 + 3.76*(alpha*z*TargetAtomicN)**2)
     *     /gbeta2
      endif
      end
!
      subroutine ckaic2(z, leng, xc2)
      implicit none
#include "Zair.h"
!
!   note: we neglect atomic electron contribution because it is
!         considered in Moller or Bhabha scattering.
!
!         compute Xc^2 = 4Pi r_0^2 N0 z^2 rho Z^2/A  * integral
!        0 to leng of 1/beta**4/gamma**2  (radian^2)/massratio2
!
      integer z     ! input. charged particle charge is ze.
      real*8 leng  ! input. length traveled by the charge particle in m
      real*8 xc2    ! output. Xc^2 in radian^2

      real*8  rho, gbeta2, beta2, massratio2
      common /Zcmedia/ rho, gbeta2, beta2,
     *       massratio2

      real*8 r0, avoganum, const, large, pi
      parameter (r0=2.817d-15, avoganum=6.022d23, pi= 3.1415)
      parameter (const = 4.*pi* r0**2 * avoganum*1.d3, 
     *    large = (pi/2)**2)
!
!
!      integeral 0 to leng of  1/(beta**4 E**2)
!      is approximated as  leng/(beta1**2 gamma1 beta2**2 gamma2)
!      Note: bata**2 *gamma = gamma - 1/gamma
!

      if(gbeta2 .le. 0.) then
         xc2 = large
      else

         xc2 = const* z* z * massratio2 *
!                    < Z(Z+1)> /<A>
     *   (TargetZ2 + TargetAtomicN)/TargetMassN 
     *   * rho * leng/gbeta2   
      endif
      end
