c     ****************************************************************
c     *                                                              *
c     * dpdedxe:  gives -de/dx  (gev/(g/cm2) of   e-/e+
c     *                                                              *
c
c     Full -  moller or bhabha loss
c
      subroutine epdedxe(sh,  aPtcl, dedt)
      implicit none
#include "Zstern.h"
#include "Zptcl.h"
#include "Zmass.h"

	include "ZdEdxSpec.h"

      record /sternh/sh    !  input. Sternheimer's consts
      record /ptcl/ aPtcl  !  input. a Particle (e- or e+) 
      real*8  dedt ! output.  dE/dt GeV/(g/cm^2)

      real*8 emass
      real*8 E, gi, Beta2, x, cb
c      real*8 dltx, wm, wlg, g, D, u, F, y
      real*8 dltx, g, D, u, F, y


      real*8 bbbeta2, lindbeta2, truebeta2
      real*8 bbbeta, lindbeta
      real*8 logbbbeta, loglindbeta
      real*8 a, b, c, xx, gra
      real*8 full, tm, vc, B1, B2, B3, B4 
      parameter (bbbeta = 0.1d0,  lindbeta = 0.01d0, gra=4.0d0/3.d0)
      parameter (bbbeta2 = bbbeta**2, lindbeta2=lindbeta**2)
      parameter (logbbbeta =-2.302585093E+00,
     *           loglindbeta = -4.605170186E+00)
      parameter (a = (1.+gra)/2.d0/(logbbbeta-loglindbeta),
     *           b = 2*a*loglindbeta + 1.) 
c

      parameter (emass = masele*1000.d0)

c////////////
c      real*8  em
c      integer iii
c      data iii /0/
c      save iii
c      if(iii .lt. 5) then
c         write(*,*) ' sh.sa=',sh.sa, ' c=',sh.c,  ' x1=',sh.x1,
c     *   ' sh.x0=', sh.x0
c
c      sh.sa =0.3755
c      sh.c = -4.435
c      sh.x1 = 2.5
c      sh.x0 = 0.0966
c      em = 2.72
c         write(*,*) ' sh.sa=',sh.sa, ' c=',sh.c,  ' x1=',sh.x1,
c     *   ' sh.x0=', sh.x0
c      iii = 1
c      endif   
c////////////

c       Energy, mass=emass  in MeV unit
      E = aPtcl.fm.p(4)*1000.d0
      g = E/emass
      gi= 1.d0/g
      truebeta2= 1. - gi**2
      if(truebeta2 .lt. bbbeta2) then
         Beta2 = bbbeta2
c         wm = emass*Beta2/2
         g = (1.+Beta2/2)
         E = emass*g
      else
         Beta2 = truebeta2
      endif
      u = g -1.0
c       x=log10(p/mc)
c      x=log10( (E/emass)**2 - 1. ) / 2  ! = log10(gbeta) = 0.4343log(gbeta)
      x=log10( g**2 - 1. ) / 2  ! = log10(gbeta) = 0.4343log(gbeta)
      cb=-sh.c
c
c        we define deltx = -delta + 4.605x = -delta + log(g+1)(g-1)
c                          -delta + log(g^2beta^2)
c
      if(x .lt. sh.x0) then
c        4.605x - dlt  
         dltx=4.605*x
      elseif(x .lt. sh.x1) then
         dltx=cb - (sh.x1-x)**3 * sh.sa
      else
         dltx=cb
      endif
c       dE/dx = sh.a/Beta2* (sh.b +log(2(g+1)) + F -delta)
c     () =  sh.b + log2(g+1) + F + deltx - log(g+1)(g-1) 
c        =  sh.b + log2 + F-log(g-1) + deltx
c        =  sh.b + 0.69315 + F - log(u)
c
      if(aPtcl.charge .eq. -1) then
         tm = min(u*emass/2, sh.w0)
c          electron
         full =  sh.b - 0.69315*(2.*gi - gi**2)
     *           + dltx -Beta2 +1.0 + (1.0-gi)**2/8.0
     *           + log(u*emass)
         if(u/2 .gt. sh.w0/emass) then
           vc = sh.w0/emass/u
           F =
     *      log(0.5/vc)
     *     - (2*g-1.)*gi*(0.5-vc) + (1.-vc/(1.0-vc))
     *     + ((g-1)*gi)**2* (-log(2*(1.-vc)) +  0.5*(0.25-vc**2))
           full = full -F
        endif
      else
c         positron
         tm = min(u*emass, sh.w0)
         vc = sh.w0/emass/u
         full =  sh.b - 0.69315*(2.*gi - gi**2)
     *           + dltx -Beta2 +1.0 + (1.0-gi)**2/8.0
     *           + log(u*emass)
         if(u .gt. sh.w0/emass) then
           y = 1.d0/(g+1.d0)
           B1 = 2.0-y**2
           B2 = (1.0-2*y)*(3.0+y**2)
           B3 = (1.0-2*y)**2 *( 1. + (1.-2*y))
           B4 = (1.0-2*y)**3
           F = log(1./vc) -Beta2*
     *   (     B1*(1.-vc) - B2/2.0*(1.0-vc**2) + B3/3.0*(1.-vc**3)
     *         - B4/4.0*(1.-vc**4) )
           full = full - F
        endif
      endif

      Tupper = tm*emass/1000.0  ! in GeV; used in Urban

      dedt =sh.a/Beta2 *full


      if(truebeta2 .lt. bbbeta2) then
         c = log(dedt) + ( a* logbbbeta - b )*logbbbeta
         if( truebeta2 .gt.  lindbeta2) then
            xx = log( truebeta2 )/2.
            dedt =exp( (-a*xx + b)*xx + c)
         else
            dedt = exp( (-a*loglindbeta + b)* loglindbeta + c) *
     *             sqrt(truebeta2)/lindbeta
         endif
      endif
c///////////////
c      write(23,*) sngl( aPtcl.fm.p(4)  - aPtcl.mass), sngl(dedt),
c     *      sngl(sqrt(truebeta2)), aPtcl.code, aPtcl.charge
c
c///////////////
c          convert it to gev/(g/cm2)
       dedt=dedt *1.d-3
       end





