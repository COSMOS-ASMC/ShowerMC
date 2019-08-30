!     ****************************************************************
!     *                                                              *
!     * epdedxNone: gives -de/dx  (gev/(g/cm2)) of non e+/e-
!     *                                                              
      module dEdxNone  
!           used by dedexNone and dedxTargetEbrems
      real(8),save:: E, mass, Beta2
      real(8),save:: wm,  u
      real(8),save::truebeta2
      real(8), save:: gb2, g

      end module dEdxNone

#include "Zcondc.h"
      subroutine epdedxNone(media, aPtcl, dedt, dedtfull )
      use dEdxNone
      use modEMcontrol
      implicit none
#include "Zmedia.h"
#include "Zcode.h"
#include "Zptcl.h"
#include "Zmass.h"
!  #include "ZepTrackp.h"
      include "ZdEdxSpec.h"  ! for epUrban. smarter solution will be using
                            !   module dEdxNone
!       MeV unit electron mass.
       type(epmedia):: media         ! input. 
       type(ptcl):: aPtcl        ! input. a particle

      real*8 dedt      ! output. restricted Energy loss rate.   GeV / (g/cm2). 
      real*8 dedtfull  ! output full dE/dx
!       MeV unit electron mass.
      real(8),parameter::emass = masele*1000.d0
      real(8),parameter::emass2 = emass**2

      real(8)::  c, xx, integ, x, full, restricted, temp, tempsqrt
      real(8):: F, cF, delta
      real(8),parameter:: bbbeta = 0.045d0
      real(8),parameter:: lindbeta = 0.0045d0
#if defined MATHLOUSY
      real(8),parameter:: logbbbeta =-3.10109278921181729466d0
      real(8),parameter:: loglindbeta =-5.40367788220586297868d0
#else
      real(8),parameter:: logbbbeta =log(bbbeta)
      real(8),parameter:: loglindbeta =log(lindbeta)
#endif
      real(8),parameter:: gra=5.d0/3.d0
      real(8),parameter:: bbbeta2= bbbeta**2
      real(8),parameter:: lindbeta2=lindbeta**2

      real(8),parameter:: a = (1.+gra)/2.d0/(logbbbeta-loglindbeta)
      real(8),parameter:: b = 2*a*loglindbeta + 1.


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
      Tupper = wm/1000.         !  in GeV; used in Urban
!        wm in unit of Me
      u = wm/emass 
!          first compute full average dE/dx
!           sh.a/Beta2( sh.b +ln(2*g^2b^2wm/m) -2Beta^2 -delta
!                 +spin_term )
! 
      full = media%sh%b +  log(2*u*gb2) -2.0*Beta2
!!!!!!!!!!!!!!!
!      write(0,*) ' sh a, b ', media%sh%a,media%sh%b
!      write(0,*) ' sh W0 ', media%sh%a,media%sh%w0
!!!!!!!!!!!!         
!            assume spin 0 particle is only pi, K
      if(aPtcl%code .ne. kpion .and. aPtcl%code .ne. kkaon) then
!          spin 1/2 term; (almost negligible)
         full = full +  (wm/E)**2/4.
      endif
!         if(aPtcl.code .eq. kmuon .and. E .gt. 5000. ) then ! 
!                   <= v9.131
!            atomic electron brems term.  at 5GeV, ~0.4 % 100GeV 2%
!            so we neglect below 5GeV
!               sh.a*alpha/(2pi) (log(2g)-1/3 log(2u)log^2(2u)
!                assuming sh.a/Beta2 = sh.a at E>5GeV
!                compute effect without sh.a and add later to full
!            pi,K,p, scales as E*sqrt(mu/M)
!             see if restricted energy is requested
      if(wm .gt. media%sh%w0) then
!             yes. requested
!                subtract average loss rate from  Ek>w0 region
!                loss for Ek>w0
         integ = log(wm/media%sh%w0) - Beta2*(1.0-media%sh%w0/wm)
!            assume spin 0 particle is only pi, K
         if(aPtcl%code .ne. kpion .and. aPtcl%code .ne. kkaon) then
!               mu, p, etc. spin = 1/2
            integ = integ + ((wm/E)**2-(media%sh%w0/E)**2)/4
         endif
      else
         integ = 0.
      endif
      call epdEdxDenC(media, g, delta)
      full = full -delta
      F =  integ 
      restricted = full - F

      dedt = media%sh%a/Beta2*restricted
!!!!!!!!!!!
!      write(0,*)' a: dedt=',dedt, ' a  Beta2 rest ',
!     * media%sh%a, Beta2, restricted      
!!!!!!!!!!!!      
      full = full*  media%sh%a/Beta2

      if(truebeta2 .lt. bbbeta2) then
         cF  = ( a* logbbbeta - b )*logbbbeta
         c = log(dedt) + cF
         if(F > 0. ) then
            cF = log(full) + cF
         else
            cF = c
         endif
         if( truebeta2 .gt.  lindbeta2) then
            xx = log( truebeta2 )/2.
            dedt =exp( (-a*xx + b)*xx + c)
            if( F > 0.) then
               full = exp( (-a*xx + b)*xx + cF)
            else
               full = dedt
            endif
               
         else
            tempsqrt = sqrt(truebeta2)/lindbeta

            dedt = exp( (-a*loglindbeta + b)* loglindbeta + c) *
     *             tempsqrt
            if(F> 0.) then
               full = exp( (-a*loglindbeta + b)* loglindbeta + cF) *
     *             tempsqrt
            else
               full = dedt
            endif
         endif
      endif
!         x  Z**2 and to GeV unit
      dedt=dedt * aPtcl%charge**2 * 1.d-3
      dedtfull = full * aPtcl%charge**2 * 1.d-3
      end
      subroutine epdedxTargetBrems(media, aPtcl, dedt, dedtfull )
!           for mu, pi, K, p.  compute dE/dx due to brems by atomic
!      electrons in the medium.
!
      use modEMcontrol
      use dEdxNone
      implicit none
#include "Zmedia.h"
#include "Zcode.h"
#include "Zptcl.h"
#include "Zmass.h"
! #include "ZepTrackp.h"
       type(epmedia):: media         ! input. 
       type(ptcl):: aPtcl        ! input. a particle
      real*8 dedt      ! output. restricted Energy loss rate.   GeV / (g/cm2). 
      real*8 dedtfull  ! output full dE/dx



!       MeV unit electron mass.
      real(8),parameter::emass = masele*1000.d0
      real(8),parameter::emass2 = emass**2
 
      real*8  full, restricted, temp
      real*8  atomicEbrem, atomicEbremCut

      logical:: Incatomicbrems

      if( TargetElecBrems == 0 ) then
         dedt= 0.
         dedtfull = 0.
         return
      endif
!               energy in MeV unit

!          first compute full average dE/dx
      if( aPtcl%code == kmuon ) then
         Incatomicbrems = E > 5000.0 .and. 
     *          btest(TargetElecBrems, 0) 
      elseif(  btest(TargetElecBrems, 1) ) then
         Incatomicbrems = E*sqrt(masmu/aPtcl%mass) > 5000.
      else
         Incatomicbrems = .false.
      endif
      if(Incatomicbrems) then
         temp = log(2*u)
         atomicEbrem = 
     *       0.00116*(log(2*g)-0.3333*temp)*temp*temp  ! 0.001167=alpha/(2pi)
      else
         atomicEbrem = 0.
      endif
!             see if restricted energy is requested
      if(wm .gt. media%sh%w0) then
!             yes. requested
         if( Incatomicbrems ) then
!                 Integ(0~wm) =Integ(0~w0) + Integ(w0~wm)
!             so  Integ(w0~wm) = Ineg(0~wm)-Integ(0~w0)
            temp = log(2*media%sh%w0/emass)
            atomicEbremCut = atomicEbrem-
     *          0.00116*(log(2*g)-0.3333*temp)*temp*temp
         else
            atomicEbremCut = 0.
         endif               
      else
         atomicEbremCut = 0.
      endif
      full =  atomicEbrem
      restricted = full -  atomicEbremCut

      dedt = media%sh%a*restricted  
      full = full*  media%sh%a
!         x  Z**2 and to GeV unit, though Z=1 normally
      dedt=dedt * aPtcl%charge**2 * 1.d-3
      dedtfull = full * aPtcl%charge**2 * 1.d-3
      end
