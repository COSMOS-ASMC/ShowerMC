      subroutine cKnockp(aPtcl,  prob, path)
      implicit none
#include "Zptcl.h"
#include "Zmass.h"
#include "ZdedxAir.h"


      real*8 cnstm 
      parameter( cnstm = 2.74296d0 * 2.d0)

!      compute prob/r.l for knock-on process of charged particles
!      other than e+/e-
!
      type(ptcl)::aPtcl  ! input. energy, charge  and mass are used.
      real*8 prob  ! output prob. of knock-on  per r%l
      real*8 path  ! output sampled path.

      real*8   g, temp,  erg, mass, u

      if(jdef .eq. 0) then
         call cerrorMsg('cdedxEleci must be called beforehand',0)
      endif
!
      mass = aPtcl%mass
      erg = aPtcl%fm%p(4)

!      tcut = w0
      g = erg/mass
      betasq = 1.0d0 - 1.d0/g**2
      tmax = 2*masele*(g**2-1.)/
     *       (1.0 + 2*g* (masele/mass) 
     *             + (masele/mass)**2)

      if(w0  .ge. tmax) then
         prob = 0.d0
         path = 1.d50
      else
         norm = (1.0d0/w0 - 1.0d0/tmax) 
         temp = norm -
     *        betasq * log(tmax/w0)/tmax 
!            assume spin 0 particle is only pi, K
         if(mass .lt. 0.11d0 .or. mass .gt. 0.9d0) then
!           mu, p, etc. spin = 1/2
            temp = temp + (tmax - w0)/erg/erg/2.d0
         endif
!                          constm = basearea*2
         prob =
     *      temp * masele * cnstm *aPtcl%charge**2 / betasq 
         call rndc(u)
         path = -log(u)/prob
      endif
      end
!    ***********
      subroutine cKnockea(aPtcl, erg2, erge, cos1, cosr)
!    **************
      implicit none
#include "Zptcl.h"
#include "Zmass.h"
#include "ZdedxAir.h"
!     
      type(ptcl)::aPtcl ! input

      real*8 erg2  ! output scattered total energy of the particle  in GeV.
      real*8 erge  ! output recoiled total electron energy in GeV
!
      real*8 cos1  ! output. cos angle of the scattered particle
      real*8 cosr  ! output. cos angle of the recoiled electron



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
         t = 1.0d0/(1.0d0/w0 - norm*u)
         call rndc(u)
         temp =1.d0 - betasq*t/tmax
         if(mass .lt. 0.11d0 .or. mass .gt. 0.9d0) then
            temp = temp + t*t/erg**2/2.d0
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
       if(p12 .le. 0.d0 ) then
          cos1 = -1.d0
          cosr = 1.d0
       else
          cos1 = (erg*erg2 - mass**2 -masele*t)/
     *     sqrt(p02*p12)      
          if(cos1 .gt. 1.d0) then
             cos1 = 1.d0
          elseif(cos1 .lt. -1.d0) then
             cos1 = -1.d0
          endif
          cosr = (erg*(t+masele) - masele*erg2)/
     *         sqrt(p02*pe2)      
          if(cosr .gt. 1.d0) then
             cosr = 1.d0
          elseif(cosr .lt. -1.d0) then
             cosr = -1.d0
          endif
       endif
      end










