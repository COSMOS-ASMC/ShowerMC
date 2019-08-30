!     ****************************************************************
!     *                                                              *
!     * dpdedxe:  gives -de/dx  (gev/(g/cm2) of   e-/e+
!     *                                                              *
!     Restricted energy loss rate dE/dt 
!         Full -  moller or bhabha loss(>RecoilkEmin)
!
      subroutine epdedxe(media,  aPtcl, dedt, dedtfull)
      implicit none
#include "Zmedia.h"
#include "Zptcl.h"
#include "Zmass.h"

	include "ZdEdxSpec.h"

       type(epmedia):: media    !  input. media
       type(ptcl)::  aPtcl  !  input. a Particle (e- or e+) 
      real*8  dedt ! output.  restricted dE/dt GeV/(g/cm^2) 
      real*8  dedtfull ! output  full dE/dt //

      real*8 emass
      real*8 E, gi, Beta2, x, cb
!      real*8 dltx, wm, wlg, g, D, u, F, y
      real*8 dltx, g, D, u, F, y

      real*8 delta

      real*8 bbbeta2, lindbeta2, truebeta2
      real*8 bbbeta, lindbeta
      real*8 logbbbeta, loglindbeta
      real*8 a, b, c, xx, gra, restricted, cF
      real*8 full, tm, vc, B1, B2, B3, B4, g1
      real*8 RKEmin ! w0/emass
      real*8 ln2, tln2, tempsqrt
      parameter (ln2=0.6931471, tln2=2*ln2)

      parameter (bbbeta = 0.1d0,  lindbeta = 0.01d0, gra=4.0d0/3.d0)
      parameter (bbbeta2 = bbbeta**2, lindbeta2=lindbeta**2)
      parameter (logbbbeta =-2.302585093E+00,
     *           loglindbeta = -4.605170186E+00)
      parameter (a = (1.+gra)/2.d0/(logbbbeta-loglindbeta),
     *           b = 2*a*loglindbeta + 1.) 
!

      parameter (emass = masele*1000.d0)   ! in MeV

!       Energy, mass=emass  in MeV unit
      E = aPtcl%fm%p(4)*1000.d0
      g = E/emass
      gi= 1.d0/g
      truebeta2= 1. - gi**2
      if(truebeta2 .lt. bbbeta2) then
         Beta2 = bbbeta2
!         wm = emass*Beta2/2
         g = (1.+Beta2/2)
         E = emass*g
      else
         Beta2 = truebeta2
      endif
      g1 =g + 1.0
      u = g - 1.0  ! incident kinetic energy in Me unit
!       x=log10(p/mc)
!      x=log10( (E/emass)**2 - 1. ) / 2  ! = log10(gbeta) = 0.4343log(gbeta)
      x=log10( g**2 - 1. ) / 2  ! = log10(gbeta) = 0.4343log(gbeta)
!        from v8.80, dltx is not used but usual delta is used.

!       dE/dx = sh.a/Beta2* (B0 +sh.b -delta +other)
!             other=shell correction < 0.75 % effect at samll E
!
!         We get restricted dE/dx; i.e knock-on K.E<sh.w0
!        First get full dE/dx and subtract dE/dx (>sh.w0)
!      (get B0 + sh.b  part)

      RKEmin = media%sh%w0/emass
      if(aPtcl%charge .eq. -1) then
!          electron; max recoil kinetic E
         tm = min(u/2, RKEmin)

         full = media%sh%b + log(u**2*g1/2) +
     *             (1.0+u**2/8. -(u+g)*ln2)/g/g

         if(u/2 .gt. RKEmin) then
!             possible  max recoil > w0; need subtraction 
            vc = RKEmin/u
!               loss in w0~tm; integrate Moller 
            F = log(0.5/vc)
     *       - (u+g)*gi*gi*(0.5-vc) + (1.-vc/(1.0-vc))
     *        +  0.5*(u*gi)**2 * (0.25-vc**2)
     *       - log(2*(1.-vc)) *(1.+(u+g)*gi*gi)
!                next is wrong
!     *       - (u+g)*gi*(0.5-vc) + (1.-vc/(1.0-vc))  
!     *       + (u*gi)**2* (-log(2*(1.-vc)) +  0.5*(0.25-vc**2))
         else
            F =0. 
         endif
      elseif(aPtcl%charge .eq. 1) then
!         positron
         tm = min(u, RKEmin)
         vc = RKEmin/u
         full =  media%sh%b  +   log(u**2*g1/2) + tln2 -
     *   Beta2/12. * (((4./g1 + 10.)/g1 + 14.0)/g1 + 23.)
         if(u .gt. RKEmin) then
            y = 1.d0/g1
            B1 = 2.0-y**2
            B2 = (1.0-2*y)*(3.0+y**2)
            B3 = (1.0-2*y)**2 *( 1. + (1.-2*y))
            B4 = (1.0-2*y)**3
            F = log(1./vc) -Beta2*
     *       (B1*(1.-vc) - B2/2.0*(1.0-vc**2) + B3/3.0*(1.-vc**3)
     *         - B4/4.0*(1.-vc**4) )
         else
            F = 0.
         endif
      else
         write(0,*) ' charge is invalid for epdedxe=', aPtcl%charge
         write(0,*) ' ptcl code =', aPtcl%code, aPtcl%subcode
         stop
      endif

      Tupper = tm*emass/1000.0  ! in GeV; used in Urban
!       get delta
      call epdEdxDenC(media, g, delta)
      full = full - delta
      restricted = full - F
      dedt =media%sh%a/Beta2 *restricted
!  ***********
      full = full * media%sh%a/Beta2
!  ***********

      if(truebeta2 .lt. bbbeta2) then
         cF =  ( a* logbbbeta - b )*logbbbeta
         c = log(dedt) + cF
         if( F > 0.) then
            cF = cF + log(full)
         else
            cF = c
         endif
         if( truebeta2 .gt.  lindbeta2) then
            xx = log( truebeta2 )/2.
            dedt = exp( (-a*xx + b)*xx + c)
            if( F > 0. ) then
               full =exp( (-a*xx + b)*xx + cF)
            else
               full = dedt
            endif
         else
            tempsqrt = sqrt(truebeta2)/lindbeta
            dedt = exp( (-a*loglindbeta + b)* loglindbeta + c) *
     *             tempsqrt
            if( F > 0. ) then
               full = exp( (-a*loglindbeta + b)* loglindbeta + cF) *
     *              tempsqrt
            else
               full = dedt
            endif
         endif
      endif
!          convert it to gev/(g/cm2)
      dedt = dedt *1.d-3
      dedtfull = full * 1.d-3
      end
