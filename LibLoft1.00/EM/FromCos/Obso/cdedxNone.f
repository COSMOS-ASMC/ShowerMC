!     ****************************************************************
!     *                                                              *
!     * epdedxNone: gives -de/dx  (gev/(g/cm2)) of non e+/e-
!     *                                                              *
!     
!
!
      subroutine cdedxNone(aPtcl, rhoin, dedt, dedtF)
      implicit none
#include  "Zcode.h"
#include "Zptcl.h"
#include "Zmass.h"
#include  "ZdedxAir.h"
      real*8 emass, emass2 
!       MeV unit electron mass.
      parameter(emass = masele*1000., emass2=emass**2)

      real*8  rhoin             ! input. density of air in kg/m3
      type(ptcl)::aPtcl        ! input. a particle

      real*8 dedt               ! output. restricedd energy loss rate.   GeV / (g/cm2). 
      real*8 dedtF     ! ouptut. full energy loss rate.
 
      real*8  E,  mass,  Beta2, x, full, restricted, temp
      real*8  wm, wlg, u, delta, atomicEbrem, atomicEbremCut
      real*8 bbbeta2, lindbeta2, truebeta2
      real*8 bbbeta, lindbeta 
      real*8 logbbbeta, loglindbeta
      real*8 a, b, c, xx, gra, gb2, g, integ
      parameter (bbbeta = 0.1d0,  lindbeta = 0.005d0, gra=5.d0/3.d0)
      parameter (bbbeta2 = bbbeta**2, lindbeta2=lindbeta**2)
      parameter (logbbbeta =-2.302585093E+00,
!     *           loglindbeta =-4.605170186E+00)
     *           loglindbeta =-5.298317367E+00)
      parameter (a = (1.+gra)/2.d0/(logbbbeta-loglindbeta),
     *           b = 2*a*loglindbeta + 1.) 

      real*8  sha, shb
!       sha = 0.153536* (media.Z/media.A) = 0.153536*0.49919=0.07664   
!       shb= 2*log(masele*/I) = 2*log(0.511e6/85.7)=                   
      data  sha/0.076643d0/, shb/17.38654d0/
      save   sha, shb

!
!               energy in MeV unit
      E = aPtcl%fm%p(4)*1000.
      mass= aPtcl%mass*1000.
      g = E/mass
      truebeta2 = 1. -(1.0d0/g)**2
      if(truebeta2 .lt. bbbeta2) then
!            fix at beta=0.1
         Beta2 = bbbeta2
         E = mass*(1. + Beta2/2)
         g = E/mass
      else
         Beta2 = truebeta2
      endif
!     
      gb2 = Beta2 * g**2   !  g^2 b^2

!          x=log10(p/mc)  = log10(g*beta)
      x=log10( gb2 ) / 2


!             max kinetic energy of knock-on
      wm = 2* emass * gb2
     *     /(1.0 + 2.0*g*(emass/mass) +(emass/mass)**2)
!           dE/dx fluctuation is not needed in Cosmos
!      Tupper = wm/1000.         !  in GeV; used in Urban
!        wm in unit of Me
      u = wm/emass 
!          first compute full average dE/dx
!           sh.a/Beta2( sh.b +ln(2*g^2b^2wm/m) -2Beta^2 -delta
!                 +spin_term )
! 
      full = shb +  log(2*u*gb2) -2.0*Beta2 
      atomicEbrem = 0.
!            assume spin 0 particle is only pi, K
      if(aPtcl%code .ne. kpion .and. aPtcl%code .ne. kkaon) then
!          spin 1/2 term; (almost negligible)
         full = full +  (wm/E)**2/4.
         if(aPtcl%code .eq. kmuon .and. E .gt. 5000. ) then
!            atomic electron brems term.  at 5GeV, ~0.4 % 100GeV 2%
!            so we neglect below 5GeV
!               sh.a*alpha/(2pi) (log(2g)-1/3 log(2u)log^2(2u)
!                assuming sh.a/Beta2 = sh.a at E>5GeV
!                compute effect without sh.a and add later to full
            temp = log(2*u)
            atomicEbrem = 
     *       0.00116*(log(2*g)-0.3333*temp)*temp*temp
         endif
      endif
!             see if restricted energy is requested
      atomicEbremCut = 0.
      if(wm .gt. w0inMeV) then
!             yes. requested
!                subtract average loss rate from  Ek>w0 region
!                loss for Ek>w0
         integ = log(wm/w0inMeV) - Beta2*(1.0-w0inMeV/wm)
!            assume spin 0 particle is only pi, K
         if(aPtcl%code .ne. kpion .and. aPtcl%code .ne. kkaon) then
!               mu, p, etc. spin = 1/2
            integ = integ + ((wm/E)**2-(w0inMeV/E)**2)/4.
            if(aPtcl%code .eq. kmuon  .and. E .gt. 5000. ) then
!                 Integ(0~wm) =Integ(0~w0) + Integ(w0~wm)
!             so  Integ(w0~wm) = Ineg(0~wm)-Integ(0~w0)
               temp = log(2*w0inMeV/emass)
               atomicEbremCut = atomicEbrem-
     *          0.00116*(log(2*g)-0.3333*temp)*temp*temp
            endif               
         endif
      else
         integ = 0.
      endif
      call cdedxdlt(rhoin, g, delta)
      full = full -delta + atomicEbrem
      restricted = full - integ  - atomicEbremCut

      dedt = sha/Beta2*restricted  
      dedtF = sha/Beta2*full


      if(truebeta2 .lt. bbbeta2) then
         c = log(dedt) + ( a* logbbbeta - b )*logbbbeta
         if( truebeta2 .gt.  lindbeta2) then
            xx = log( truebeta2 )/2.
            dedt =exp( (-a*xx + b)*xx + c)
         else
            dedt = exp( (-a*loglindbeta + b)* loglindbeta + c) *
     *             sqrt(truebeta2)/lindbeta
         endif
         dedtF = dedt
      endif
!         x  Z**2 and to GeV unit
      dedt=dedt * aPtcl%charge**2 * 1.d-3
      dedtF=dedtF * aPtcl%charge**2 * 1.d-3
      end
