      subroutine epKnockp(media, aPtcl,  prob, path)
      implicit none
#include "Zmedia.h"
#include "Zptcl.h"
#include "Zmass.h"
      common /KnockC/ beta2, tmax, norm, tcut
      real*8 beta2, tmax, norm, tcut


!      compute prob/r.l for knock-on process of charged particles
!      other than e+/e-
!
       type(epmedia):: media ! input. medeia. media.sh.w0 
       type(ptcl):: aPtcl  ! input. energy and mass are used.
      real*8 prob  ! output prob. of knock-on  per r.l
      real*8 path  ! output sampled path.

      real*8   g, temp,  erg, mass, u

!
      mass = aPtcl%mass
      erg = aPtcl%fm%p(4)

      tcut = media%sh%tcut
      g = erg/mass
      beta2 = 1.0 - 1.d0/g**2
      tmax = 2*masele*(g**2-1.)/
     *       (1.0 + 2*g* (masele/mass) 
     *             + (masele/mass)**2)

      if(tcut .ge. tmax) then
         prob = 1.d-35
      else
         norm = (1.0d0/tcut - 1.0d0/tmax) 
         temp = norm -
     *       beta2 * log(tmax/tcut)/tmax 
!            assume spin 0 particle is only pi, K
         if(mass .lt. 0.11 .or. mass .gt. 0.9) then
!           mu, p, etc. spin = 1/2
            temp = temp + (tmax - tcut)/erg/erg/2
         endif
!                          constm = basearea*2
         prob =
     *      temp * masele * media%basearea * 2 / beta2
     *      * aPtcl%charge**2
      endif
      call rndc(u)
      path = -log(u)/prob
      end
!    ***********
      subroutine epKnockea(aPtcl, erg2, erge, cos1, cosr)
!     actual prodction of e and saving it to Pwork is done in epNEPKnock using
!     this subroutine     
!    **************
      implicit none
#include "Zptcl.h"
#include "Zmass.h"
!     
       type(ptcl):: aPtcl ! input

      real*8 erg2  ! output otal energy of the scattered incident particle  in GeV.
      real*8 erge  ! output ejected electron total energy in GeV
!
      real*8 cos1  ! output. cos angle of the scattered particle
      real*8 cosr  ! output. cos angle of the recoiled electron

      common /KnockC/ beta2, tmax, norm, tcut
      real*8 beta2, tmax, norm, tcut

!       samples delta ray electron energy erge
!       and survival partcle energy and their angles
!
!      this erg is lower than ergsave due to
!      energy loss in the path

      real*8 mass, erg, u, temp, t
      logical accept
      real*8  p02, p12, pe2

      mass = aPtcl%mass
      erg = aPtcl%fm%p(4)

      accept = .false.
      do while(.not. accept)
         call rndc(u)
         t = 1.0d0/(1.0d0/tcut - norm*u)
         call rndc(u)
         temp =1. - beta2*t/tmax
         if(mass .lt. 0.11 .or. mass .gt. 0.9) then
            temp = temp + t*t/erg**2/2
         endif
         if(temp .gt. u) then
            accept = .true.
         endif
      enddo
      erge = t + masele

!     erg2 = olderg - erge + masele
!          above will give larger energy than
!          current erg and energy eventually increase.
!          avoid it by using current erg
! 
      erg2 = max( erg - erge + masele, mass)

!
!      from
!       p1cos1 + p2cos2 = p0
!       p1sin1 = p2sin2
!       E0+m = E1 + E2
!     we obtain
!       cos1 = (p0^2 + p1^2 - p2^2)/(2p0p1)
!       cos2 = (p0^2 + p2^2 - p1^2)/(2p0p2)
!        (p0^2 + p1^2 - p2^2)/2 = E0E1 -M^2 - m(E2-m)
!        (p0^2 + p2^2 - p1^2)/2 = E0E2 - mE1
!
!      p02 =abs(
       p02 =   erg**2 - mass**2
       p12 =   erg2**2 -mass**2
       pe2 =   erge**2 - masele**2
       if(p12 .le. 0. ) then
          cos1 = -1.
          cosr = 1.
       else
          cos1 = (erg*erg2 - mass**2 -masele*t)/
     *     sqrt(p02*p12)      
          if(cos1 .gt. 1.) then
             cos1 = 1.
          elseif(cos1 .lt. -1.) then
             cos1 = -1.
          endif
          cosr = (erg*(t+masele) - masele*erg2)/
     *         sqrt(p02*pe2)      
          if(cosr .gt. 1.) then
             cosr = 1.
          elseif(cosr .lt. -1.) then
             cosr = -1.
          endif
       endif
      end










